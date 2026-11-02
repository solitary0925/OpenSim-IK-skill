"""Portable input-block diagnostics and accounting; never repair or discard data."""
import csv
import json
from pathlib import Path


def finish_preflight(result):
    global_errors=result.setdefault('global_errors',[])
    for row in result['trials']:
        row['status']='excluded_by_user' if row['excluded'] else ('blocked_global' if global_errors else ('blocked_preflight' if row['errors'] else 'runnable'))
    result['counts']={s:sum(r['status']==s for r in result['trials']) for s in ('runnable','blocked_preflight','blocked_global','excluded_by_user')}
    result['counts']['total']=len(result['trials'])
    result['ready']=not result['errors']
    result['can_run']=not global_errors and result['counts']['runnable']>0
    return result


def table_details(table,names,alignment,selected,timing=None):
    import numpy as np
    times=np.asarray(list(table.getIndependentColumn()));n=table.getNumRows();selected=np.asarray(selected,dtype=bool)
    if len(times)!=n or len(selected)!=n:raise ValueError('Input row/time/mask counts differ')
    valid_times=bool(np.isfinite(times).all() and np.all(np.diff(times)>0))
    def ranges(mask):
        edges=np.flatnonzero(np.diff(np.r_[False,mask,False]));rows=[]
        for a,b in zip(edges[::2],edges[1::2]):
            rows.append(dict(start_frame=int(a+1),end_frame=int(b),total_frames=n,frames=int(b-a),
                             start_s=float(times[a]) if valid_times else None,end_s=float(times[b-1]) if valid_times else None,
                             start_source_frame=int(timing['source_first_frame']+a) if timing else None,
                             end_source_frame=int(timing['source_first_frame']+b-1) if timing else None))
        return rows
    complete=np.ones(n,dtype=bool);markers=[]
    for name in names:
        actual=alignment['mapping'].get(name)
        if actual is None:complete[:]=False;continue
        col=table.getDependentColumn(actual);values=np.array([[col[i][k] for k in range(3)] for i in range(n)],dtype=float).reshape(n,3)
        bad=~np.isfinite(values).all(axis=1)|np.all(values==0,axis=1);complete&=~bad
        if bad.any():markers.append(dict(model_marker=name,c3d_marker=actual,missing_frames=int(bad.sum()),selected_missing_frames=int((bad&selected).sum()),missing_ranges=ranges(bad),selected_missing_ranges=ranges(bad&selected)))
    return dict(total_frames=n,selected_frames=int(selected.sum()),selected_ranges=ranges(selected),
                frame_basis='1-based sample ordinal in the full source C3D, inclusive; not gait-cycle percentage',
                source_first_frame=timing['source_first_frame'] if timing else None,
                source_last_frame=timing['source_last_frame'] if timing else None,
                unmatched_model_markers=[r['model_name'] for r in alignment['rows'] if r['vicon_name'] is None],
                unresolved_correspondences=[r for r in alignment['rows'] if r['vicon_name'] is None],
                missing_markers=markers,fully_observed_ranges=ranges(complete),selected_fully_observed_ranges=ranges(complete&selected),
                note='Complete marker ranges are diagnostics only; no automatic window selection, interpolation, anonymous-point matching or additional joint locking.')


def write_blocked_reports(result,directory,subject='subject'):
    p=Path(directory);p.mkdir(parents=True,exist_ok=True)
    blocked=[r for r in result['trials'] if r['status'] in ('blocked_preflight','blocked_global')]
    output=dict(subject=subject,counts=result['counts'],global_errors=result.get('global_errors',[]),trials=blocked)
    (p/'blocked_inputs.json').write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')
    fields=['trial','status','model_marker','c3d_marker','scope','start_frame','end_frame','total_frames','frames','start_source_frame','end_source_frame','start_s','end_s']
    with (p/'输入受阻帧范围.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for row in blocked:
            detail=row.get('input_details',{})
            for marker in detail.get('unmatched_model_markers',[]):
                writer.writerow(dict(trial=row['trial'],status=row['status'],model_marker=marker,scope='unresolved_name',total_frames=detail.get('total_frames')))
            for marker in detail.get('missing_markers',[]):
                for scope,key in [('full_source','missing_ranges'),('selected_window','selected_missing_ranges')]:
                    for segment in marker[key]:writer.writerow(dict(trial=row['trial'],status=row['status'],model_marker=marker['model_marker'],c3d_marker=marker['c3d_marker'],scope=scope,**segment))
    text=f"# {subject} 输入受阻明细\n\n计数：{json.dumps(result['counts'],ensure_ascii=False)}。全局问题：{result.get('global_errors',[])}。\n\n受阻试次不启动IK；正常试次在无全局错误时继续运行。受阻、用户确认排除、求解失败及已求解但待复核分别记录。未删除源文件、插值、改名、移动点位或自动裁剪；没有自动把受阻试次列入废弃清单。完整标记窗口仅供复核，不代表已确认步态周期。\n\n| 试次 | 状态 | 总帧数 | 缺少/歧义点名 | 缺失样本（全段；所选窗口） | 原因 |\n|---|---|---:|---|---|---|\n"
    for row in blocked:
        d=row.get('input_details',{});missing=', '.join(d.get('unmatched_model_markers',[])) or '—'
        gaps='；'.join(f"{r['model_marker']}: {r['missing_frames']}/{d['total_frames']}；窗口内{r['selected_missing_frames']}帧" for r in d.get('missing_markers',[])) or '—'
        reason='；'.join(result.get('global_errors',[]) if row['status']=='blocked_global' else row['errors'])
        text+=f"| {row['trial']} | {row['status']} | {d.get('total_frames','—')} | {missing} | {gaps} | {reason} |\n"
    if not blocked:text+='| — | 无受阻试次 | — | — | — | — |\n'
    text+='\n逐个缺失片段的整条试次帧序号/总帧数、原始C3D帧号及adapter时间见输入受阻帧范围.csv；仅在原始时间映射已验证时填写采集帧号。全段和所选窗口分别记录，两个scope不可相加计数。修复输入后重新预检查；通过检查的试次可重新运行。旧结果保留，但当前受阻状态禁止其被当作本次有效结果使用。\n'
    (p/'输入受阻明细.md').write_text(text,encoding='utf-8')
    return output
