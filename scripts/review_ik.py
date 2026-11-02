"""Audit saved IK outputs; localize errors, diagnose label hypotheses, honor exclusions."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import opensim as o
from review_context import add_context, context_markdown, load_windows, position_text, ranges_text


def exclusions(path,known=None):
    if path is None:return set()
    values=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(values,list) or not all(isinstance(v,str) and v and '/' not in v and '\\' not in v and Path(v).suffix.lower()=='.c3d' for v in values):
        raise ValueError('Exclusions must be a JSON array of C3D basenames.')
    names={v.casefold() for v in values}
    if len(names)!=len(values):raise ValueError('Duplicate exclusions.')
    if known is not None and names-{v.casefold() for v in known}:raise ValueError('Exclusion names absent from this batch: '+str(names-{v.casefold() for v in known}))
    return names


def episodes(mask,times,values):
    edges=np.flatnonzero(np.diff(np.r_[False,mask,False]))
    return [dict(start_s=float(times[a]),end_s=float(times[b-1]),frames=int(b-a),peak_mm=float(np.max(values[a:b]))) for a,b in zip(edges[::2],edges[1::2])]


def fit_rigid(reference,observed):
    p=reference-reference.mean(0);q=observed-observed.mean(0)
    u,_,vt=np.linalg.svd(p.T@q);s=np.eye(3);s[2,2]=np.linalg.det(u@vt)
    return float(np.sqrt(np.mean(np.sum((p@u@s@vt-q)**2,axis=1))))


def review_trial(directory,rms_limit=20.,max_limit=40.,visualize=True,max_visual_frames=3,analysis_window=None):
    p=Path(directory).resolve();d=json.loads((p/'metrics.json').read_text(encoding='utf-8'));issues=[]
    def require(condition,message):
        if not condition:issues.append(message)
    with (p/'marker_errors_mm.csv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    times=np.array([float(r['time']) for r in rows]);rms=np.array([float(r['RMS_mm']) for r in rows]);maxima=np.array([float(r['Max_mm']) for r in rows]);names=list(d['weights'])
    require(len(times)==d['frames'],'Independent CSV frame count mismatch')
    require(np.isfinite(times).all() and np.all(np.diff(times)>0),'Nonfinite/non-increasing error times')
    errors=np.array([[float(r[n]) for n in names] for r in rows]);require(np.isfinite(errors).all(),'Nonfinite point errors')
    require(np.allclose(rms,np.sqrt(np.mean(errors**2,axis=1)),atol=.002),'CSV RMS does not match point columns')
    require(np.allclose(maxima,errors.max(axis=1),atol=.002),'CSV Max does not match point columns')
    native=o.TimeSeriesTable(str(p/(Path(d['trial']).stem+'_ik_marker_errors.sto')))
    require(native.getNumRows()==len(times),'Native error frame count mismatch')
    if native.getNumRows()==len(times):
        require(np.allclose(np.array(list(native.getIndependentColumn())),times,atol=1e-8),'Native error times mismatch')
        for label,target in [('marker_error_RMS',rms),('marker_error_max',maxima)]:
            col=native.getDependentColumn(label);actual=np.array([col[i]*1000 for i in range(len(times))])
            require(np.allclose(actual,target,atol=.002),'Native/independent '+label+' mismatch')
    motion=o.TimeSeriesTable(str(p/'ik.mot'));require(motion.getNumRows()==len(times),'Motion frame count mismatch')
    if motion.getNumRows()==len(times):require(np.allclose(np.array(list(motion.getIndependentColumn())),times,atol=1e-8),'Motion times mismatch')
    in_degrees=motion.getTableMetaDataAsString('inDegrees').strip().lower()=='yes'
    coordinate_columns=list(motion.getColumnLabels());work=ET.parse(p/'ik_model.osim');coords={e.get('name'):e for e in work.iter('Coordinate')}
    angular_changes=[]
    for name in coordinate_columns:
        col=motion.getDependentColumn(name);arr=np.array([col[i] for i in range(motion.getNumRows())]);require(np.isfinite(arr).all(),'Nonfinite output '+name)
        if name not in ('pelvis_tx','pelvis_ty','pelvis_tz'):
            degrees=arr if in_degrees else np.rad2deg(arr)
            for i in np.flatnonzero(abs(np.diff(degrees))>20):angular_changes.append(dict(coordinate=name,start_s=float(times[i]),end_s=float(times[i+1]),change_deg=float(degrees[i+1]-degrees[i])))
    for lock in d['locked_coordinates']:
        name=lock['coordinate'];require(name in coords and coords[name].findtext('locked','false').strip()=='true','Working model lock missing: '+name)
        require(name in coordinate_columns,'Locked coordinate missing from MOT: '+name)
        if name in coordinate_columns:
            col=motion.getDependentColumn(name);arr=np.array([col[i] for i in range(motion.getNumRows())]);expected=lock['value_deg'] if in_degrees else lock['value_rad']
            require(np.allclose(arr,expected,atol=1e-8),'Locked coordinate changed: '+name)
    for source,expected in [(d['source_model'],d['source_sha256']),(d['c3d'],d['c3d_sha256'])]:
        if Path(source).is_file():require(hashlib.sha256(Path(source).read_bytes()).hexdigest()==expected,'Source changed since IK: '+source)
        else:issues.append('Source unavailable for provenance verification: '+source)
    timing=None;time_rows=[]
    try:
        if any('Source changed since IK' in v or 'Source unavailable' in v for v in issues):raise ValueError('Cannot regenerate timing from unverified source provenance')
        from timing import write_timing
        adapter=o.C3DFileAdapter();table=adapter.getMarkersTable(adapter.read(d['c3d']))
        timing,time_rows=write_timing(p,d['c3d'],list(table.getIndependentColumn()),times)
    except Exception as e:issues.append('Timing verification: '+str(e))
    segments=[]
    for marker_index,name in enumerate(names):
        for s in episodes(errors[:,marker_index]>max_limit,times,errors[:,marker_index]):segments.append(dict(kind='point_max',marker=name,**s))
    for s in episodes(rms>rms_limit,times,rms):segments.append(dict(kind='frame_rms',marker='ALL',**s))
    for segment in segments:
        for edge in ('start','end'):
            index=int(np.argmin(abs(times-segment[edge+'_s'])))
            segment[edge+'_source_frame']=int(time_rows[index]['source_frame']) if time_rows else None
    hypotheses=[];skipped=[]
    if segments:
        alignment=d['alignment'];identity={alignment['identities'][n]:target for n,target in alignment['mapping'].items()}
        marker_defs={e.get('name'):e for e in work.iter('Marker')};trc=o.TimeSeriesTableVec3(str(p/'markers.trc'))
        for side in ('L','R'):
            canonical=[side+n for n in ('HEE','TOE','SMH','FMH','VMH')];foot=[identity.get(n) for n in canonical]
            if any(n not in names or n not in marker_defs for n in foot):skipped.append(side+': required five foot markers unavailable');continue
            parents={marker_defs[n].findtext('socket_parent_frame') for n in foot}
            if len(parents)!=1:skipped.append(side+': points not on one rigid body; rigid-foot diagnostic inapplicable');continue
            reference=np.array([list(map(float,marker_defs[n].findtext('location').split())) for n in foot])*1000
            data=np.array([[[trc.getDependentColumn(n)[i][k]*1000 for k in range(3)] for n in foot] for i in range(trc.getNumRows())])
            current=float(np.mean([fit_rigid(reference,row) for row in data]));swapped=float(np.mean([fit_rigid(reference,row[[0,1,2,4,3]]) for row in data]))
            if current>15 and swapped<current*.5 and current-swapped>10:
                hypotheses.append(dict(side=side,markers=[foot[3],foot[4]],current_rigid_rms_mm=current,swap_hypothesis_rigid_rms_mm=swapped,confirmed=False))
    result=dict(trial=d['trial'],classification='needs_review' if issues or segments or angular_changes else 'no_numerical_flags',
                verification_passed=not issues,verification_issues=issues,thresholds=dict(frame_rms_mm=rms_limit,point_max_mm=max_limit,adjacent_rotation_deg=20),
                frames=len(times),global_rms_mm=float(np.sqrt(np.mean(errors**2))),max_error_mm=float(maxima.max()),worst_marker=d['worst_marker'],
                abnormal_segments=segments,angular_change_flags=angular_changes,label_hypotheses=hypotheses,diagnostic_skips=skipped,timing=timing)
    add_context(result,times,errors,rms,maxima,time_rows,analysis_window)
    (p/'review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    with (p/'异常时间段.csv').open('w',newline='',encoding='utf-8-sig') as f:
        fields=['kind','marker','start_s','end_s','start_source_frame','end_source_frame','frames','peak_mm','start_ik_row','end_ik_row','ik_total_frames','start_trial_frame','end_trial_frame','trial_total_frames','start_trial_percent','end_trial_percent','analysis_relation','analysis_overlap_frames']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(segments)
    rows_md='\n'.join(f"| {s['marker']} | {s['kind']} | {position_text(s)} | {s['start_source_frame']}–{s['end_source_frame']} | {s['start_s']:.3f}–{s['end_s']:.3f} | {s['frames']} | {s['peak_mm']:.2f} | {s['analysis_relation']} |" for s in segments)
    context=context_markdown(result)
    text=f'''# {d['trial']} IK核查

状态：**{result['classification']}**；输出一致性检查：{'通过' if not issues else '需处理'}。问题：{issues}。

{context}

“第几帧/总帧数”从原始C3D文件首个存储样本起，以1计数，起止均包含；原始C3D帧号另列。百分比为帧序号/整条帧数，不代表步态周期百分比。

| 点 | 类型 | 整条试次中的位置 | 原始C3D帧号 | adapter时间s | 帧数 | 峰值mm | 分析区间关系 |
|---|---|---|---|---|---:|---:|---|
{rows_md or '| — | 无超阈值片段 | — | — | — | 0 | — | — |'}

时间为C3D adapter时间，裁剪窗口不重新归零。原始帧号和所选偏移见time_mapping.csv、timing.json；帧号时间不是墙钟，也不证明外部测力同步。时间信息核查：{'通过' if timing else '失败，见问题列表'}。RMS/Max阈值是筛查参考，不是硬性验收标准；少数异常帧不自动判整段废弃。角度跳变：{angular_changes}。

点名互换诊断假设：{hypotheses}。这是五个同侧刚体点的离线配准比较，不是修正后的IK结果；不能只凭误差降低确认点名。未交换数据、改点、插值、裁剪或修改权重。诊断跳过项：{skipped}。

`no_numerical_flags`仅表示未触发本轮数值筛查，仍需运动可视化和解剖合理性验证。锁定角是模型假设。详情见review.json及异常时间段.csv。
'''
    (p/'核查报告.md').write_text(text,encoding='utf-8')
    if visualize:attach_visuals(p,result,max_visual_frames)
    return result


def attach_visuals(directory,result,max_frames=3):
    p=Path(directory)
    if result['classification']=='needs_review':
        try:
            from visualize_ik import render_overlays, render_timeline
            result['timeline']=render_timeline(p,result)
            result['visualizations']=render_overlays(p,result,max_frames)
            with (p/'核查报告.md').open('a',encoding='utf-8') as f:f.write('\n![整条试次误差位置](异常帧可视化/error_timeline.png)\n')
            note='\n异常帧叠加图见 [图像索引](异常帧可视化/index.md)。固定脚本输出，仅供复核，不自动修改点位。\n' if result['visualizations'] else '\n本次未选出可绘制的异常帧。\n'
        except Exception as e:
            result['visualization_error']=str(e);note='\n异常帧可视化失败：'+str(e)+'。\n'
        with (p/'核查报告.md').open('a',encoding='utf-8') as f:f.write(note)
    (p/'review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')


def review_batch(results,exclude_list=None,visualize=True,max_visual_frames=3,analysis_windows=None,preflight_report=None):
    root=Path(results).resolve();records=json.loads((root/'batch_metrics.json').read_text(encoding='utf-8'))
    selected=Path(exclude_list) if exclude_list else root/'excluded_trials.json';excluded=exclusions(selected,[r['trial'] for r in records]) if selected.is_file() else set()
    if exclude_list and not selected.is_file():raise FileNotFoundError(selected)
    from review_context import analysis_windows as validate_windows
    windows=validate_windows(analysis_windows or {},[r['trial'] for r in records])
    reviews=[];checks={r['trial'].casefold():r for r in preflight_report['trials']} if preflight_report else {}
    for record in records:
        directory=root/Path(record['trial']).stem
        check=checks.get(record['trial'].casefold(),{})
        status=check.get('status') if check.get('status','').startswith('blocked_') else record.get('status')
        blocked=status in ('blocked_preflight','blocked_global');failed=status=='failed'
        if blocked or failed:
            messages=check.get('errors') or ([record['error']] if record.get('error') else [])
            if status=='blocked_global' and preflight_report:messages=preflight_report['global_errors']
            r=dict(trial=record['trial'],classification=status if blocked else 'failed_or_unavailable',verification_passed=False,verification_issues=messages,
                   input_details=check.get('input_details') or record.get('input_details'),previous_results_preserved=(directory/'metrics.json').is_file())
        elif (directory/'metrics.json').is_file():
            try:r=review_trial(directory,visualize=False,analysis_window=windows.get(record['trial'].casefold()))
            except Exception as e:r=dict(trial=record['trial'],classification='needs_review',verification_passed=False,verification_issues=[str(e)])
        else:r=dict(trial=record['trial'],classification='failed_or_unavailable',verification_passed=False)
        if record['trial'].casefold() in excluded:
            r['classification_before_exclusion']=r['classification'];r['classification']='excluded_by_user';r['exclusion_reason']='User-confirmed discarded trial; do not use in subsequent analysis.'
            if (directory/'核查报告.md').is_file():
                with (directory/'核查报告.md').open('a',encoding='utf-8') as f:f.write('\n批次最终选择：**excluded_by_user**。用户确认废弃，保留上述历史诊断，不进入后续分析或新增默认可视化。\n')
        if (directory/'metrics.json').is_file() and not (blocked or failed):
            if visualize:attach_visuals(directory,r,max_visual_frames)
            else:(directory/'review.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
        reviews.append(r)
    counts={s:sum(r['classification']==s for r in reviews) for s in sorted({r['classification'] for r in reviews})}
    output=dict(counts=counts,exclude_list=str(selected) if selected.is_file() else None,trials=reviews)
    (root/'review_status.json').write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')
    (root/'trial_selection.json').write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')
    rows=[]
    for r in reviews:
        analysis=r.get('analysis_context',{}).get('classification','not_evaluated')
        trial=r['trial'];link=Path(trial).stem+'/核查报告.md'
        name=f'[{trial}]({link})' if r['classification'] in ('needs_review','no_numerical_flags') else trial
        detail='；'.join(r.get('verification_issues',[]))
        rows.append(f"| {name} | {r['classification']} | {r.get('global_rms_mm',float('nan')):.2f} | {r.get('max_error_mm',float('nan')):.2f} | {r.get('worst_marker','—')} | {ranges_text(r)} | {analysis} | {detail} |")
    text='# IK核查汇总\n\n计数：'+json.dumps(counts,ensure_ascii=False)+'。\n\n| 试验 | 整条试次状态 | 整体RMS mm | Max mm | 最大误差点 | 异常位置：帧序号/整条总帧数 | 分析区间筛查状态 | 受阻/核查问题 |\n|---|---|---:|---:|---|---|---|---|\n'+'\n'.join(rows)+'\n\n帧序号从C3D文件首个样本起，以1计数，起止均包含；不是原始采集帧号或步态周期百分比。仅评价已计算的IK样本。未配置确认的分析区间时不推断力台接触或自动忽略起止端误差。区间外异常作为背景保留，区间内统计另列，整条试次状态和排除清单不变。\n\nblocked_preflight：单试次输入/输出策略受阻，不启动IK；blocked_global：环境/模型/全局配置受阻，整批不启动。受阻明细见输入受阻明细.md及blocked_inputs.json。旧结果保留，但当前受阻记录不会从旧MOT自动认定完成。failed_or_unavailable：求解失败或输出不可用。no_numerical_flags：未触发数值筛查，仍待可视化验证。needs_review：查看对应试验核查报告及异常时间段。excluded_by_user：只按用户明确确认的清单排除，不删除源文件或既有结果；诊断不会自动废弃试次或确认标签互换。核查不会修改MOT/C3D、点名、模型、权重或裁剪数据。\n'
    (root/'核查汇总.md').write_text(text,encoding='utf-8');return output


def main():
    p=argparse.ArgumentParser(description=__doc__);group=p.add_mutually_exclusive_group(required=True);group.add_argument('--trial');group.add_argument('--results');p.add_argument('--exclude-list');p.add_argument('--no-visuals',action='store_true');p.add_argument('--max-visual-frames',type=int,default=3);p.add_argument('--analysis-windows');a=p.parse_args()
    if not 1<=a.max_visual_frames<=20:p.error('--max-visual-frames must be 1..20')
    windows=load_windows(a.analysis_windows)
    if a.trial:
        trial=json.loads((Path(a.trial)/'metrics.json').read_text(encoding='utf-8'))['trial']
        if set(windows)-{trial.casefold()}:p.error('Analysis windows contain other trials')
        r=review_trial(a.trial,visualize=not a.no_visuals,max_visual_frames=a.max_visual_frames,analysis_window=windows.get(trial.casefold()))
    else:r=review_batch(a.results,a.exclude_list,not a.no_visuals,a.max_visual_frames,windows)
    print(json.dumps(r.get('counts',dict(classification=r.get('classification'))),ensure_ascii=False))


if __name__=='__main__':main()
