"""Portable JSON entry point: preflight first, then solve or audit saved IK."""
import argparse
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

DEFAULT_ROTATION=[-1,0,0,0,0,1,0,1,0]


def load_config(path):
    path=Path(path).resolve();c=json.loads(path.read_text(encoding='utf-8-sig'))
    allowed={'schema_version','subject','model','data','c3d','output','geometry_dir','window','weights','marker_map','rotation','exclusions','workers','reuse_existing','overwrite','visualization','analysis_windows'}
    if not isinstance(c,dict) or set(c)-allowed:raise ValueError('Unknown config fields: '+str(set(c)-allowed if isinstance(c,dict) else 'not an object'))
    if type(c.get('schema_version')) is not int or c['schema_version']!=1:raise ValueError('schema_version must be integer 1')
    for key in ('subject','model','output'):
        if not isinstance(c.get(key),str) or not c[key].strip():raise ValueError('Missing/invalid '+key)
    if ('data' in c)==('c3d' in c):raise ValueError('Specify exactly one of data (batch directory) or c3d (single file)')
    for key in ('model','data','c3d','output','geometry_dir'):
        if key in c:
            if not isinstance(c[key],str) or not c[key]:raise ValueError('Invalid path '+key)
            c[key]=str((path.parent/c[key]).resolve())
    c.setdefault('window',{});c.setdefault('weights',{});c.setdefault('marker_map',{});c.setdefault('rotation',DEFAULT_ROTATION);c.setdefault('workers',2);c.setdefault('reuse_existing',False);c.setdefault('overwrite',False);c.setdefault('visualization',{})
    if not isinstance(c['window'],dict) or set(c['window'])-{'start','end'}:raise ValueError('window supports start/end in adapter seconds only')
    for k,v in c['window'].items():
        if v is not None and (isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v)):raise ValueError('Invalid window '+k)
    for name in ('weights','marker_map'):
        if not isinstance(c[name],dict):raise ValueError(name+' must be an object')
    for name,value in c['weights'].items():
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:raise ValueError('Weight must be finite and positive: '+name)
    if not all(isinstance(v,str) and v for v in c['marker_map'].values()):raise ValueError('marker_map values must be labels')
    if not isinstance(c['rotation'],list) or len(c['rotation'])!=9 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in c['rotation']):raise ValueError('rotation must contain 9 finite row-major numbers')
    for key in ('reuse_existing','overwrite'):
        if not isinstance(c[key],bool):raise ValueError(key+' must be boolean')
    if type(c['workers']) is not int or c['workers'] not in (1,2):raise ValueError('workers must be 1 or 2')
    v=c['visualization']
    if not isinstance(v,dict) or set(v)-{'enabled','max_frames_per_trial'}:raise ValueError('Invalid visualization fields')
    v.setdefault('enabled',True);v.setdefault('max_frames_per_trial',3)
    if not isinstance(v['enabled'],bool) or type(v['max_frames_per_trial']) is not int or not 1<=v['max_frames_per_trial']<=20:raise ValueError('Invalid visualization settings')
    # Existing user exclusions persist unless the study manifest is deliberately edited.
    persisted=Path(c['output'])/'excluded_trials.json'
    old=json.loads(persisted.read_text(encoding='utf-8-sig')) if persisted.is_file() else []
    c.setdefault('exclusions',old)
    if not isinstance(c['exclusions'],list) or not all(isinstance(n,str) and n and '/' not in n and '\\' not in n and Path(n).suffix.lower()=='.c3d' for n in c['exclusions']):raise ValueError('exclusions must contain C3D basenames')
    if len({n.casefold() for n in c['exclusions']})!=len(c['exclusions']):raise ValueError('Duplicate exclusions')
    if {n.casefold() for n in old}-{n.casefold() for n in c['exclusions']}:raise ValueError('Config would re-include saved user exclusions; reconcile the study manifest first')
    from review_context import analysis_windows
    c['analysis_windows']=analysis_windows(c.get('analysis_windows',{}))
    return c,path


def preflight(c,action='solve'):
    from blocked_inputs import finish_preflight,table_details
    result=dict(ready=False,action=action,errors=[],global_errors=[],warnings=[],trials=[],environment=dict(python=sys.executable));errors=result['errors'];warnings=result['warnings']
    def global_error(message):errors.append(message);result['global_errors'].append(message)
    # Import OpenSim in scratch space so check-only never creates a solver log in inputs/results.
    old_cwd=Path.cwd()
    with tempfile.TemporaryDirectory(prefix='ik_preflight_') as scratch:
        try:
            os.chdir(scratch)
            try:
                np=importlib.import_module('numpy');o=importlib.import_module('opensim')
                result['environment'].update(opensim=o.GetVersion(),numpy=np.__version__)
                if c['visualization']['enabled']:importlib.import_module('matplotlib')
            except Exception as e:global_error('Runtime dependency: '+str(e));return finish_preflight(result)
            try:
                source=Path(c['model']);tree=ET.parse(source);names=[m.get('name') for m in tree.iter('Marker')]
                if not names or len(names)!=len(set(names)):raise ValueError('Missing/duplicate model marker names')
                for side in ('r','l'):
                    joint=tree.find(f".//CustomJoint[@name='knee_{side}']")
                    if joint is None or [v.get('name') for v in joint.findall('./coordinates/Coordinate')]!=['knee_angle_'+side]:raise ValueError('Original single-rotation knee model required')
                rotation=np.array(c['rotation']).reshape(3,3)
                if not np.allclose(rotation@rotation.T,np.eye(3)) or not np.isclose(np.linalg.det(rotation),1):raise ValueError('rotation must be proper orthogonal')
                if c.get('geometry_dir'):o.ModelVisualizer.addDirToGeometrySearchPaths(c['geometry_dir'])
                o.ModelVisualizer.addDirToGeometrySearchPaths(str(source.parent/'Geometry'))
                from run_ik import assess_observability
                result['observability']=assess_observability(o.Model(str(source)),names)
                result['model_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
                missing=[]
                for mesh in {e.text.strip() for e in tree.iter('mesh_file') if e.text}:
                    options=[source.parent/mesh,source.parent/'Geometry'/mesh]
                    if c.get('geometry_dir'):options.append(Path(c['geometry_dir'])/mesh)
                    if not any(p.is_file() for p in options):missing.append(mesh)
                if missing:warnings.append('Display geometry missing (marker IK can run): '+str(missing))
            except Exception as e:global_error('Model/config: '+str(e));return finish_preflight(result)
            folder=Path(c.get('data',c.get('c3d')))
            files=sorted([p for p in folder.iterdir() if p.is_file() and p.suffix.lower()=='.c3d' and 'static' not in p.stem.lower()],key=lambda p:p.name.lower()) if 'data' in c and folder.is_dir() else ([folder] if 'c3d' in c else [])
            if not files:global_error('No dynamic C3D files found')
            known={f.name.casefold() for f in files};excluded={n.casefold() for n in c['exclusions']}
            if excluded-known:global_error('Exclusions absent from inputs: '+str(excluded-known))
            if set(c.get('analysis_windows',{}))-known:global_error('Analysis windows refer to unknown trials')
            out=Path(c['output'])
            if out in (source.parent,source,folder) or any(out==f for f in files):global_error('Output must be separate from source model/data')
            if len(files) and len(excluded)==len(files):warnings.append('All selected trials are excluded; no IK will run')
            from marker_alignment import plan
            from timing import write_timing
            for f in files:
                row=dict(trial=f.name,file=str(f),excluded=f.name.casefold() in excluded,errors=[]);result['trials'].append(row)
                if row['excluded']:continue
                try:
                    adapter=o.C3DFileAdapter();t=adapter.getMarkersTable(adapter.read(str(f)))
                    alignment=plan(names,list(t.getColumnLabels()),c['marker_map'])
                    if not alignment['ready']:row['errors'].extend(alignment['errors'])
                    if t.getTableMetaDataAsString('Units').strip().lower()!='mm':row['errors'].append('Expected mm marker units')
                    times=np.array(list(t.getIndependentColumn()));valid_time=len(times)>=2 and np.isfinite(times).all() and np.all(np.diff(times)>0)
                    start=c['window'].get('start');end=c['window'].get('end')
                    start=times[0] if start is None and len(times) else start;end=times[-1] if end is None and len(times) else end
                    valid_window=valid_time and start is not None and end is not None and times[0]<=start<end<=times[-1]
                    if not valid_window:row['errors'].append('Invalid time vector/window')
                    mask=(times>=start-1e-9)&(times<=end+1e-9) if valid_window else np.zeros(len(times),dtype=bool)
                    selected=np.flatnonzero(mask)
                    if len(selected)<2:row['errors'].append('Too few selected frames')
                    timing=None
                    try:timing,_=write_timing(scratch,f,times,times)
                    except Exception as e:row['errors'].append('Timing verification: '+str(e))
                    detail=table_details(t,names,alignment,mask,timing);row['input_details']=detail
                    actual=list(alignment['mapping'].values())
                    missing=[m['model_marker'] for m in detail['missing_markers'] if m['selected_missing_frames']]
                    if missing:row['errors'].append('Missing selected marker samples: '+', '.join(missing))
                    unknown=set(c['weights'])-set(names)-set(actual)
                    if unknown:row['errors'].append('Unknown weight names: '+str(unknown))
                    window=c.get('analysis_windows',{}).get(f.name.casefold())
                    if window and len(selected):
                        first=int(selected[0])+1;last=int(selected[-1])+1
                        if window['end_frame']>len(times) or window['start_frame']<first or window['end_frame']>last:row['errors'].append('Analysis window lies outside source trial or computed IK coverage')
                    row.update(frames=len(selected),marker_count=len(actual))
                    if len(selected):row.update(start_s=float(times[selected[0]]),end_s=float(times[selected[-1]]),source_first_frame=int(timing['source_first_frame']+selected[0]) if timing else None,source_last_frame=int(timing['source_first_frame']+selected[-1]) if timing else None)
                    target=out/f.stem if 'data' in c else out
                    if action=='review':
                        if not (target/'metrics.json').is_file():raise ValueError('No saved metrics to review')
                    elif target.exists() and any(target.iterdir()) and not (c['overwrite'] or c['reuse_existing']):raise ValueError('Nonempty output; explicitly choose reuse_existing or authorized overwrite')
                    elif target.exists() and any(target.iterdir()) and c['reuse_existing'] and not c['overwrite']:
                        from result_policy import compatible
                        metric=target/'metrics.json'
                        if not metric.is_file():raise ValueError('Existing output has no reusable metrics')
                        d=json.loads(metric.read_text(encoding='utf-8'))
                        expected=dict(requested_start=c['window'].get('start'),requested_end=c['window'].get('end'),requested_weights=c['weights'],marker_map=c['marker_map'],accuracy=1e-5,locking_policy='zero_sensitivity_profile_v1')
                        if not compatible(target,d,result['model_sha256'],f,o.GetVersion(),c['rotation'],expected):raise ValueError('Existing result implementation/settings/source fingerprint mismatch; audit or explicitly authorize replacement')
                except Exception as e:row['errors'].append(str(e))
                errors.extend(f.name+': '+message for message in row['errors'])
            if action=='solve' and 'c3d' in c and c['reuse_existing']:global_error('reuse_existing is supported for batch only; use review-only to audit a single saved trial')
            return finish_preflight(result)
        finally:
            if 'o' in locals():o.Logger.removeFileSink()
            os.chdir(old_cwd)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);g=p.add_mutually_exclusive_group();g.add_argument('--check-only',action='store_true');g.add_argument('--review-only',action='store_true');p.add_argument('--report-dir',help='Explicit directory for preflight JSON/Markdown; otherwise check-only writes nothing');a=p.parse_args()
    try:c,path=load_config(a.config);r=preflight(c,'review' if a.review_only else 'solve')
    except Exception as e:
        from blocked_inputs import finish_preflight
        c=None;r=finish_preflight(dict(ready=False,errors=[str(e)],global_errors=[str(e)],warnings=[],trials=[]))
    print(json.dumps(r,ensure_ascii=False,indent=2))
    report_dir=Path(a.report_dir).resolve() if a.report_dir else ((Path(c['output'])/'_configuration' if 'data' in c else Path(c['output']).parent/(Path(c['output']).name+'_configuration')) if c and not a.check_only else None)
    if report_dir:
        from blocked_inputs import write_blocked_reports
        report_dir.mkdir(parents=True,exist_ok=True);(report_dir/'preflight.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
        write_blocked_reports(r,report_dir,c['subject'] if c else 'subject')
        text='# IK预检查\n\n全部通过：'+str(r['ready'])+'；存在可运行试次：'+str(r['can_run'])+'；动作：'+r.get('action','unknown')+'。未启动IK求解。计数：'+json.dumps(r['counts'],ensure_ascii=False)+'。\n\n错误：'+json.dumps(r['errors'],ensure_ascii=False)+'\n\n提示：'+json.dumps(r['warnings'],ensure_ascii=False)+'\n\n| 试次 | 状态 | 帧数 | 原始帧范围 | 问题 |\n|---|---|---:|---|---|\n'
        text+='\n'.join(f"| {v['trial']} | {v['status']} | {v.get('frames','—')} | {v.get('source_first_frame','—')}–{v.get('source_last_frame','—')} | {v['errors']} |" for v in r['trials'])
        (report_dir/'预检查报告.md').write_text(text,encoding='utf-8')
    if a.check_only:return 0 if r['ready'] else 2
    if not c or r['global_errors'] or ('c3d' in c and not r['ready']):return 2
    out=Path(c['output']);out.mkdir(parents=True,exist_ok=True);record=out/'_configuration' if 'data' in c else out.parent/(out.name+'_configuration');record.mkdir(parents=True,exist_ok=True)
    (record/'resolved_config.json').write_text(json.dumps(c,ensure_ascii=False,indent=2),encoding='utf-8')
    provenance=dict(config_file=str(path),config_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),action='review_saved_results' if a.review_only else 'solve',python=sys.executable)
    (record/'config_provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
    for key in ('weights','marker_map','exclusions','analysis_windows'):(record/(key+'.json')).write_text(json.dumps(c[key],ensure_ascii=False,indent=2),encoding='utf-8')
    if 'data' in c:(out/'excluded_trials.json').write_text(json.dumps(c['exclusions'],ensure_ascii=False,indent=2),encoding='utf-8')
    if a.review_only:
        from review_ik import review_batch,review_trial
        v=c['visualization']
        if 'data' in c:review_batch(out,record/'exclusions.json',v['enabled'],v['max_frames_per_trial'],c['analysis_windows'],preflight_report=r)
        elif Path(c['c3d']).name.casefold() not in {n.casefold() for n in c['exclusions']}:review_trial(out,visualize=v['enabled'],max_visual_frames=v['max_frames_per_trial'],analysis_window=c['analysis_windows'].get(Path(c['c3d']).name.casefold()))
        print('Reviewed saved results; no solver rerun and no settings retroactively applied.');return 0 if r['ready'] else 2
    if 'c3d' in c and Path(c['c3d']).name.casefold() in {n.casefold() for n in c['exclusions']}:print('User-excluded trial; IK skipped.');return 0
    if 'data' in c:
        from batch_ik import execute_batch
        options=argparse.Namespace(subject=c['subject'],model=c['model'],data=c['data'],output=c['output'],geometry_dir=c.get('geometry_dir'),workers=c['workers'],
                                   reuse_existing=c['reuse_existing'],overwrite=c['overwrite'],exclude_list=str(record/'exclusions.json'),start=c['window'].get('start'),end=c['window'].get('end'),
                                   weights=str(record/'weights.json'),marker_map=str(record/'marker_map.json'),rotation=c['rotation'],no_visuals=not c['visualization']['enabled'],max_visual_frames=c['visualization']['max_frames_per_trial'],analysis_windows=str(record/'analysis_windows.json'))
        return execute_batch(options,r)
    script=Path(__file__).with_name('run_ik.py')
    cmd=[sys.executable,str(script),'--subject',c['subject'],'--model',c['model'],'--output',c['output'],'--rotation']+[str(v) for v in c['rotation']]
    cmd+=['--weights',str(record/'weights.json'),'--marker-map',str(record/'marker_map.json'),'--max-visual-frames',str(c['visualization']['max_frames_per_trial'])]
    cmd+=['--analysis-windows',str(record/'analysis_windows.json')]
    cmd+=['--c3d',c['c3d']]
    if c.get('geometry_dir'):cmd+=['--geometry-dir',c['geometry_dir']]
    for key in ('start','end'):
        if c['window'].get(key) is not None:cmd+=['--'+key,str(c['window'][key])]
    if c['overwrite']:cmd+=['--overwrite']
    if not c['visualization']['enabled']:cmd+=['--no-visuals']
    with tempfile.TemporaryDirectory(prefix='ik_driver_') as scratch:return subprocess.run(cmd,cwd=scratch).returncode


if __name__=='__main__':sys.exit(main())
