"""Apply one fixed IK policy to a subject's dynamic C3D trials."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import opensim as o
from review_ik import exclusions, review_batch


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('model','data','output','subject'):p.add_argument('--'+n,required=True)
    p.add_argument('--geometry-dir');p.add_argument('--workers',type=int,default=2,choices=(1,2))
    p.add_argument('--reuse-existing',action='store_true');p.add_argument('--overwrite',action='store_true')
    p.add_argument('--exclude-list',help='User-confirmed JSON array of discarded trial basenames')
    p.add_argument('--start',type=float);p.add_argument('--end',type=float);p.add_argument('--weights');p.add_argument('--marker-map')
    p.add_argument('--rotation',nargs=9,type=float,default=[-1,0,0,0,0,1,0,1,0])
    p.add_argument('--no-visuals',action='store_true');p.add_argument('--max-visual-frames',type=int,default=3)
    p.add_argument('--analysis-windows',help='JSON per-trial review intervals; does not crop IK')
    args=p.parse_args()
    if not 1<=args.max_visual_frames<=20:p.error('--max-visual-frames must be 1..20')
    return execute_batch(args)


def execute_batch(args,preflight_report=None):
    source=Path(args.model).resolve();data=Path(args.data).resolve();out=Path(args.output).resolve()
    for key in ('weights','marker_map','geometry_dir'):
        if getattr(args,key):setattr(args,key,str(Path(getattr(args,key)).resolve()))
    if out==data or out==source.parent:raise ValueError('Use a separate output folder.')
    out.mkdir(parents=True,exist_ok=True);runner=Path(__file__).with_name('run_ik.py').resolve()
    trials=[dict(trial=f.name,file=str(f),status='pending') for f in sorted(data.glob('*.c3d'),key=lambda p:p.name.lower()) if 'static' not in f.stem.lower()]
    if not trials:raise ValueError('No dynamic C3D trials found.')
    from review_context import load_windows
    windows=load_windows(args.analysis_windows,[r['trial'] for r in trials])
    exclusion_file=Path(args.exclude_list).resolve() if args.exclude_list else out/'excluded_trials.json'
    if args.exclude_list and not exclusion_file.is_file():raise FileNotFoundError(exclusion_file)
    excluded=exclusions(exclusion_file,[r['trial'] for r in trials]) if exclusion_file.is_file() else set()
    for r in trials:
        if r['trial'].casefold() in excluded:r['status']='excluded_by_user'
    from blocked_inputs import write_blocked_reports
    if preflight_report is None:
        from run_config import load_config,preflight
        config=dict(schema_version=1,subject=args.subject,model=str(source),data=str(data),output=str(out),
                    window=dict(start=args.start,end=args.end),weights=json.loads(Path(args.weights).read_text(encoding='utf-8-sig')) if args.weights else {},
                    marker_map=json.loads(Path(args.marker_map).read_text(encoding='utf-8-sig')) if args.marker_map else {},rotation=args.rotation,
                    exclusions=[r['trial'] for r in trials if r['status']=='excluded_by_user'],workers=args.workers,reuse_existing=args.reuse_existing,overwrite=args.overwrite,
                    analysis_windows=windows,visualization=dict(enabled=not args.no_visuals,max_frames_per_trial=args.max_visual_frames))
        if args.geometry_dir:config['geometry_dir']=args.geometry_dir
        with tempfile.TemporaryDirectory(prefix='ik_batch_check_') as scratch:
            path=Path(scratch)/'config.json';path.write_text(json.dumps(config),encoding='utf-8');c,_=load_config(path);preflight_report=preflight(c)
    write_blocked_reports(preflight_report,out,args.subject)
    record=out/'_configuration';record.mkdir(exist_ok=True);(record/'preflight.json').write_text(json.dumps(preflight_report,ensure_ascii=False,indent=2),encoding='utf-8')
    checks={r['trial'].casefold():r for r in preflight_report['trials']}
    for r in trials:
        check=checks.get(r['trial'].casefold())
        if r['status']=='excluded_by_user':continue
        if preflight_report.get('global_errors'):
            r.update(status='blocked_global',error='; '.join(preflight_report['global_errors']))
        elif check is None or check['status']!='runnable':
            r.update(status='blocked_preflight',error='; '.join(check['errors']) if check else 'No valid per-trial preflight record',input_details=check.get('input_details') if check else None)
        if r['status'].startswith('blocked_'):print(r['trial'],r['status'],r['error'],flush=True)
    digest=preflight_report.get('model_sha256')
    def report():
        rows=[]
        for r in trials:
            if r['status']=='excluded_by_user':
                rows.append(f"| {r['trial']} | — | — | — | — | — | 用户确认废弃，排除 |")
                continue
            d=r.get('metrics')
            if d:
                flag='需复核' if d['frames_rms_above_20'] or d['frames_max_above_40'] or d['jump_screen'] else '待可视化验证'
                rows.append(f"| {r['trial']} | {d['frames']} | {d['global_rms_mm']:.2f} | {d['peak_frame_rms_mm']:.2f} | {d['max_marker_error_mm']:.2f} | {d['worst_marker']} | {flag} |")
            else:rows.append(f"| {r['trial']} | — | — | — | — | — | {r.get('error',r['status'])} |")
        text='# '+args.subject+' 原膝模型 IK 批次报告\n\n模型：`'+str(source)+'`。按保存配置使用匹配模型点权重和坐标转换；源模型未改动。各工作副本按无观测坐标策略锁定坐标，具体值见逐试验报告和metrics.json，保留可观测的踝和距下关节自由度。\n\n| 试验 | 帧数 | 整体RMS mm | 逐帧RMS峰值 mm | Max mm | 最大误差点 | 检查状态 |\n|---|---:|---:|---:|---:|---|---|\n'+'\n'.join(rows)+'\n\nRMS/Max为从保存的MOT独立重建标记位置计算的结果。逐帧RMS>20mm、Max>40mm或邻帧角度跳变粗筛会标为需复核；这些是经验/工程提示，不是解剖验证。锁定角是建模假设，不得当作实验测得结果。各试验保留工作模型、配置、TRC、MOT、逐点误差和日志。尚未做人工运动可视化或ID/RRA/SO/CMC/JRA。\n'
        (out/(args.subject+'_IK批次报告.md')).write_text(text,encoding='utf-8')
        (out/'batch_metrics.json').write_text(json.dumps(trials,ensure_ascii=False,indent=2),encoding='utf-8')
    def run(r):
        f=Path(r['file']);target=out/f.stem;metric=target/'metrics.json'
        if args.reuse_existing and metric.is_file():
            d=json.loads(metric.read_text(encoding='utf-8'))
            expected=dict(requested_start=args.start,requested_end=args.end,requested_weights=json.loads(Path(args.weights).read_text(encoding='utf-8-sig')) if args.weights else {},marker_map=json.loads(Path(args.marker_map).read_text(encoding='utf-8-sig')) if args.marker_map else {},accuracy=1e-5,locking_policy='zero_sensitivity_profile_v1')
            from result_policy import compatible
            if compatible(target,d,digest,f,o.GetVersion(),args.rotation,expected):return d
            if not args.overwrite:raise RuntimeError('Existing results lack a matching implementation/settings fingerprint; review them or explicitly authorize overwrite.')
        cmd=[sys.executable,str(runner),'--subject',args.subject,'--model',str(source),'--c3d',str(f),'--output',str(target)]
        cmd+=['--rotation']+[str(v) for v in args.rotation]+['--no-visuals']
        for key in ('start','end','weights','marker_map'):
            if getattr(args,key) is not None:cmd+=['--'+key.replace('_','-'),str(getattr(args,key))]
        if args.geometry_dir:cmd+=['--geometry-dir',args.geometry_dir]
        if args.overwrite:cmd+=['--overwrite']
        # The driver log is kept outside the trial so it does not affect output checks.
        with tempfile.TemporaryDirectory(prefix='ik_trial_') as scratch,(out/(f.stem+'.driver.log')).open('w',encoding='utf-8') as log:
            result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=scratch,timeout=180)
        if result.returncode:raise RuntimeError('Failed; see '+f.stem+'.driver.log')
        return json.loads(metric.read_text(encoding='utf-8'))
    report()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures={executor.submit(run,r):r for r in trials if r['status']=='pending'}
        for future in as_completed(futures):
            r=futures[future]
            try:
                r['metrics']=future.result();r['status']='completed';d=r['metrics'];print(r['trial'],round(d['global_rms_mm'],2),round(d['max_marker_error_mm'],2),flush=True)
            except Exception as e:r['status']='failed';r['error']=str(e);print(r['trial'],'FAILED',str(e),flush=True)
            report()
    review=review_batch(out,str(exclusion_file) if exclusion_file.is_file() else None,not args.no_visuals,args.max_visual_frames,windows)
    print('COMPLETE',sum(r['status']=='completed' for r in trials),'/',len(trials)-len(excluded),'EXCLUDED',len(excluded),'REVIEW',review['counts'],flush=True)
    return 2 if preflight_report.get('global_errors') or any(r['status'] in ('failed','blocked_preflight','blocked_global') for r in trials) else 0


if __name__=='__main__':sys.exit(main())
