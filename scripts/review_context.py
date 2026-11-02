"""Locate errors in the full acquisition and optional, explicitly configured ROI."""
import json
from pathlib import Path


def analysis_windows(values,known=None):
    if not isinstance(values,dict):raise ValueError('analysis_windows must be an object keyed by C3D basename')
    normalized={}
    for name,window in values.items():
        if not isinstance(name,str) or not name or '/' in name or '\\' in name or Path(name).suffix.lower()!='.c3d':raise ValueError('analysis_windows keys must be C3D basenames')
        if name.casefold() in normalized:raise ValueError('Duplicate analysis window trial: '+name)
        if not isinstance(window,dict) or set(window)-{'start_frame','end_frame','reason'}:raise ValueError('Invalid analysis window: '+name)
        start,end=window.get('start_frame'),window.get('end_frame')
        if type(start) is not int or type(end) is not int or start<1 or end<start:raise ValueError('Analysis window requires 1-based inclusive start_frame <= end_frame: '+name)
        if not isinstance(window.get('reason'),str) or not window['reason'].strip():raise ValueError('Document analysis window reason: '+name)
        normalized[name.casefold()]=dict(window)
    if known is not None and set(normalized)-{n.casefold() for n in known}:raise ValueError('Analysis windows refer to unknown trials: '+str(set(normalized)-{n.casefold() for n in known}))
    return normalized


def load_windows(path,known=None):
    return analysis_windows(json.loads(Path(path).read_text(encoding='utf-8-sig')) if path else {},known)


def add_context(result,times,errors,rms,maxima,time_rows,window=None):
    import numpy as np
    times=np.asarray(times);n=len(times);timing=result.get('timing')
    verified=bool(timing and len(time_rows)==n)
    full=int(timing['source_frames']) if verified else None
    ordinal=np.array([int(r['source_sample_index_0'])+1 for r in time_rows]) if verified else None
    mask=(np.asarray(rms)>result['thresholds']['frame_rms_mm']) | (np.asarray(maxima)>result['thresholds']['point_max_mm'])
    for change in result.get('angular_change_flags',[]):
        # A jump involves both adjacent samples, including a jump across an ROI boundary.
        for edge in ('start','end'):mask[int(np.argmin(abs(times-change[edge+'_s'])))]=True
    inside=None
    if window:
        analysis_windows({result['trial']:window})
        if not verified:raise ValueError('Analysis window requires verified original-frame mapping')
        if window['end_frame']>full or window['start_frame']<ordinal[0] or window['end_frame']>ordinal[-1]:raise ValueError('Analysis window lies outside source trial or computed IK coverage')
        inside=(ordinal>=window['start_frame'])&(ordinal<=window['end_frame'])
    for episode in result['abnormal_segments']+result.get('angular_change_flags',[]):
        a=int(np.argmin(abs(times-episode['start_s'])));b=int(np.argmin(abs(times-episode['end_s'])))
        episode.update(start_ik_row=a+1,end_ik_row=b+1,ik_total_frames=n,
                       start_trial_frame=int(ordinal[a]) if verified else None,end_trial_frame=int(ordinal[b]) if verified else None,trial_total_frames=full,
                       start_trial_percent=float(ordinal[a]/full*100) if verified else None,end_trial_percent=float(ordinal[b]/full*100) if verified else None,
                       analysis_relation='unconfigured')
        if inside is not None:
            count=int(inside[a:b+1].sum());episode['analysis_overlap_frames']=count
            episode['analysis_relation']='outside_analysis_window' if count==0 else ('inside_analysis_window' if count==b-a+1 else 'crosses_analysis_boundary')
    edges=np.flatnonzero(np.diff(np.r_[False,mask,False]))
    abnormal=[dict(start_frame=int(ordinal[a]),end_frame=int(ordinal[b-1]),total_frames=full) for a,b in zip(edges[::2],edges[1::2])] if verified else []
    location=dict(frame_basis='1-based sample ordinal in full source C3D, inclusive; NOT acquisition frame labels',trial_total_frames=full,
                  computed_ik_frames=n,computed_start_trial_frame=int(ordinal[0]) if verified else None,computed_end_trial_frame=int(ordinal[-1]) if verified else None,
                  abnormal_unique_computed_frames=int(mask.sum()),abnormal_ranges=abnormal)
    result['trial_position']=location
    if inside is None:
        result['analysis_context']=dict(status='unconfigured',classification='not_evaluated',note='No confirmed analysis/contact interval. Do not infer force-plate contact or discard errors from trial position alone.')
    else:
        count=int((mask&inside).sum());section=np.asarray(errors)[inside]
        result['analysis_context']=dict(status='configured',window=dict(window),classification='needs_review_in_analysis_window' if count or not result['verification_passed'] else 'no_numerical_flags_in_analysis_window',
                                       analyzed_frames=int(inside.sum()),abnormal_frames_inside=count,abnormal_frames_outside=int((mask&~inside).sum()),
                                       global_rms_mm=float(np.sqrt(np.mean(section**2))),max_marker_error_mm=float(section.max()),
                                       note='Review scope only; source data, MOT/TRC, exclusions and full-trial classification remain unchanged. Interval definition must come from confirmed events or an explicit study decision.')
    return result


def position_text(segment):
    if segment.get('trial_total_frames') is None:return '原始映射不可用'
    return f"第 {segment['start_trial_frame']}–{segment['end_trial_frame']} 帧 / 共 {segment['trial_total_frames']} 帧（{segment['start_trial_percent']:.1f}–{segment['end_trial_percent']:.1f}%）"


def ranges_text(result):
    return '；'.join(f"{r['start_frame']}–{r['end_frame']} / {r['total_frames']}" for r in result.get('trial_position',{}).get('abnormal_ranges',[])) or ('无' if result.get('trial_position',{}).get('trial_total_frames') else '映射不可用')


def context_markdown(result):
    p=result['trial_position'];c=result['analysis_context']
    text=f"整条C3D共 {p['trial_total_frames'] if p['trial_total_frames'] is not None else '未知'} 帧；本次IK覆盖第 {p['computed_start_trial_frame']}–{p['computed_end_trial_frame']} 帧，共 {p['computed_ik_frames']} 帧。异常位置（合并重复点/角度标记）：**{ranges_text(result)}**；异常样本共 {p['abnormal_unique_computed_frames']} 帧。只评价已计算的IK样本。\n\n"
    if c['status']=='unconfigured':return text+'未配置确认的分析保留区间，暂不判断异常是否位于力台接触/有效步态周期以外，也不因此判整条试次废弃。\n'
    w=c['window']
    return text+f"分析保留区间：第 {w['start_frame']}–{w['end_frame']} 帧 / 共 {p['trial_total_frames']} 帧；依据：{w['reason']}。区间内异常 {c['abnormal_frames_inside']} 帧，区间外异常 {c['abnormal_frames_outside']} 帧；区间内整体RMS {c['global_rms_mm']:.2f} mm，Max {c['max_marker_error_mm']:.2f} mm。分析区间筛查状态：**{c['classification']}**。区间外异常保留作背景记录；未自动裁剪或改变整条试次的状态/排除清单。\n"
