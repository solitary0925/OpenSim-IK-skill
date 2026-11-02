"""Marker-based IK for a scaled single-rotation-knee CGM24 model."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import opensim as o
from marker_alignment import plan, rename_tree, save_plan

LOCKABLE={'mtp_angle_r','mtp_angle_l','lumbar_extension','lumbar_bending','lumbar_rotation'}


def positions(model,state,names):
    model.realizePosition(state)
    return np.array([[model.getMarkerSet().get(n).getLocationInGround(state)[k] for k in range(3)] for n in names])


def assess_observability(model,names):
    state=model.initSystem();cs=model.getCoordinateSet();original=[cs.get(i).getValue(state) for i in range(cs.getSize())]
    tests=[]
    for flex in (None,-.3,-.8):
        for i,v in enumerate(original):cs.get(i).setValue(state,v,False)
        if flex is not None:
            for side in ('r','l'):
                cs.get('knee_angle_'+side).setValue(state,flex,False)
                cs.get('hip_flexion_'+side).setValue(state,.2,False)
        columns=[]
        for i in range(cs.getSize()):
            c=cs.get(i);v=c.getValue(state);h=1e-5
            c.setValue(state,v+h,False);plus=positions(model,state,names)
            c.setValue(state,v-h,False);minus=positions(model,state,names)
            c.setValue(state,v,False);columns.append(((plus-minus)/(2*h)).ravel())
        tests.append(np.column_stack(columns))
    norms=np.max(np.array([np.linalg.norm(j,axis=0) for j in tests]),axis=0)
    zeros=[cs.get(i).getName() for i,v in enumerate(norms) if v<1e-8 and not cs.get(i).getDefaultLocked()]
    unexpected=set(zeros)-LOCKABLE
    if unexpected:raise ValueError('Unexpected unobservable coordinates require review: '+str(unexpected))
    active=[i for i in range(cs.getSize()) if cs.get(i).getName() not in zeros and not cs.get(i).getDefaultLocked()]
    ranks=[int(np.linalg.matrix_rank(j[:,active],tol=1e-6)) for j in tests]
    if min(ranks)<len(active):raise ValueError('Tracking markers do not independently constrain remaining coordinates: '+str(ranks))
    return dict(auto_locked=zeros,active_coordinates=len(active),ranks=ranks,coordinate_names=[cs.get(i).getName() for i in range(cs.getSize())],max_marker_sensitivity=norms.tolist(),pose_tests=[None,-.3,-.8])


def write_trc(path,times,names,values,rate,timing):
    with path.open('w',encoding='utf-8') as f:
        f.write(f'PathFileType\t4\t(X/Y/Z)\t{path.name}\nDataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n')
        f.write(f"{rate}\t{rate}\t{len(times)}\t{len(names)}\tm\t{rate}\t{timing['source_first_frame']}\t{timing['source_frames']}\nFrame#\tTime\t")
        f.write('\t\t\t'.join(names)+'\t\t\n\t\t'+'\t'.join(f'X{i}\tY{i}\tZ{i}' for i in range(1,len(names)+1))+'\n\n')
        for i,(t,row) in enumerate(zip(times,values)):
            f.write(f"{timing['selected_first_frame']+i}\t{t:.10g}\t"+'\t'.join(f'{v:.12g}' for v in row.ravel())+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('model','c3d','output'):parser.add_argument('--'+key,required=True)
    parser.add_argument('--subject',default='subject');parser.add_argument('--geometry-dir')
    parser.add_argument('--start',type=float);parser.add_argument('--end',type=float)
    parser.add_argument('--marker-map');parser.add_argument('--weights',help='JSON model/C3D marker name -> positive weight; default 1')
    parser.add_argument('--rotation',nargs=9,type=float,default=[-1,0,0,0,0,1,0,1,0])
    parser.add_argument('--no-visuals',action='store_true');parser.add_argument('--max-visual-frames',type=int,default=3)
    parser.add_argument('--analysis-windows',help='JSON per-trial review intervals; does not crop IK')
    parser.add_argument('--overwrite',action='store_true');args=parser.parse_args()
    if not 1<=args.max_visual_frames<=20:parser.error('--max-visual-frames must be 1..20')
    for key in ('weights','marker_map','geometry_dir','analysis_windows'):
        if getattr(args,key):setattr(args,key,str(Path(getattr(args,key)).resolve()))
    model_path,c3d,out=[Path(getattr(args,n)).resolve() for n in ('model','c3d','output')]
    if out in (model_path.parent,c3d.parent) or out in (model_path,c3d):raise ValueError('Use a separate output directory.')
    if out.exists() and any(out.iterdir()) and not args.overwrite:raise ValueError('Nonempty output needs authorized --overwrite.')
    model_hash=hashlib.sha256(model_path.read_bytes()).hexdigest();tree=ET.parse(model_path)
    for side in ('r','l'):
        joint=tree.find(f".//CustomJoint[@name='knee_{side}']")
        if joint is None or [c.get('name') for c in joint.findall('./coordinates/Coordinate')]!=['knee_angle_'+side]:
            raise ValueError('This workflow requires the original single-rotation knee model.')
    rotation=np.array(args.rotation).reshape(3,3)
    if not np.allclose(rotation@rotation.T,np.eye(3)) or not np.isclose(np.linalg.det(rotation),1):raise ValueError('Rotation must be proper orthogonal.')
    a=o.C3DFileAdapter();t=a.getMarkersTable(a.read(str(c3d)));raw_names=[m.get('name') for m in tree.iter('Marker')]
    mapping=json.loads(Path(args.marker_map).read_text(encoding='utf-8-sig')) if args.marker_map else None
    alignment=plan(raw_names,list(t.getColumnLabels()),mapping);save_plan(alignment,out)
    if not alignment['ready']:raise ValueError('Unresolved names; see 点名称对齐报告.md.')
    rename_tree(tree,alignment['mapping']);names=[alignment['mapping'][n] for n in raw_names]
    if t.getTableMetaDataAsString('Units').strip().lower()!='mm':raise ValueError('Validated converter requires marker units mm.')
    times=np.array(list(t.getIndependentColumn()))
    if len(times)<2 or not np.isfinite(times).all() or not np.all(np.diff(times)>0):raise ValueError('Invalid C3D time vector.')
    start=times[0] if args.start is None else args.start;end=times[-1] if args.end is None else args.end
    if start<times[0] or end>times[-1] or not start<end:raise ValueError('Selected time range is outside data.')
    selected=(times>=start-1e-9)&(times<=end+1e-9)
    if selected.sum()<2:raise ValueError('Too few selected samples.')
    values=np.empty((selected.sum(),len(names),3))
    for j,name in enumerate(names):
        col=t.getDependentColumn(name);arr=np.array([[col[int(i)][k] for k in range(3)] for i in np.flatnonzero(selected)])
        if not np.isfinite(arr).all() or np.any(np.all(arr==0,axis=1)):raise ValueError('Missing marker values: '+name)
        values[:,j]=arr@rotation.T*.001
    from timing import write_timing
    timing,_=write_timing(out,c3d,times,times[selected])
    from review_context import load_windows
    analysis_window=load_windows(args.analysis_windows,[c3d.name]).get(c3d.name.casefold())
    if analysis_window and (analysis_window['start_frame']<int(np.flatnonzero(selected)[0])+1 or analysis_window['end_frame']>int(np.flatnonzero(selected)[-1])+1):raise ValueError('Analysis window lies outside source trial or computed IK coverage')
    times=times[selected];weights={n:1. for n in names}
    if args.weights:
        for name,value in json.loads(Path(args.weights).read_text(encoding='utf-8-sig')).items():
            n=alignment['mapping'].get(name,name)
            if n not in weights or not math.isfinite(value) or value<=0:raise ValueError('Unknown marker/nonpositive weight: '+name)
            weights[n]=float(value)
    # Original model stays untouched; modifications are confined to this work copy.
    out.mkdir(parents=True,exist_ok=True);os.chdir(out);geometry=out/'Geometry';geometry.mkdir(exist_ok=True);missing=[]
    for name in {e.text.strip() for e in tree.iter('mesh_file') if e.text}:
        candidates=[model_path.parent/name,model_path.parent/'Geometry'/name]
        if args.geometry_dir:candidates.append(Path(args.geometry_dir)/name)
        source=next((p for p in candidates if p.is_file()),None)
        if source:shutil.copy2(source,geometry/Path(name).name)
        else:missing.append(name)
    o.ModelVisualizer.addDirToGeometrySearchPaths(str(geometry))
    aligned=out/'aligned_model.osim';tree.write(aligned,encoding='utf-8',xml_declaration=True)
    model=o.Model(str(aligned));observability=assess_observability(model,names);cs=model.updCoordinateSet();locked=[]
    for name in observability['auto_locked']:
        c=cs.get(name);v=c.getDefaultValue();c.setDefaultLocked(True);c.setDefaultSpeedValue(0)
        locked.append(dict(coordinate=name,value_rad=v,value_deg=math.degrees(v),reason='No marker sensitivity in three tested poses; explicit user policy.'))
    working=out/'ik_model.osim';model.printToXML(str(working));working_model=o.Model(str(working));s=working_model.initSystem()
    assert all(working_model.getCoordinateSet().get(r['coordinate']).getLocked(s) for r in locked)
    write_trc(out/'markers.trc',times,names,values,float(t.getTableMetaDataAsString('DataRate')),timing)
    tool=o.InverseKinematicsTool();tool.setName(c3d.stem);tool.set_model_file(str(working));tool.setMarkerDataFileName(str(out/'markers.trc'))
    tool.setStartTime(float(times[0]));tool.setEndTime(float(times[-1]));tool.setResultsDir(str(out));tool.setOutputMotionFileName(str(out/'ik.mot'))
    tool.set_accuracy(1e-5);tool.set_report_errors(True);tool.set_report_marker_locations(True)
    tasks=tool.upd_IKTaskSet()
    for name in names:
        task=o.IKMarkerTask();task.setName(name);task.setApply(True);task.setWeight(weights[name]);tasks.cloneAndAppend(task)
    tool.printToXML(str(out/'ik_setup.xml'))
    log_path=out/'ik_run.log';o.Logger.removeFileSink();log_path.write_text('',encoding='utf-8');o.Logger.addFileSink(str(log_path))
    try:
        if not tool.run():raise RuntimeError('IK returned failure.')
    finally:o.Logger.removeFileSink()
    # Independently reconstruct marker errors from the saved MOT and work model.
    motion=o.TimeSeriesTable(str(out/'ik.mot'));mt=np.array(list(motion.getIndependentColumn()));columns=list(motion.getColumnLabels())
    if len(mt)!=len(times) or not np.allclose(mt,times,atol=1e-8):raise RuntimeError('IK output times differ from selected input.')
    in_degrees=motion.getTableMetaDataAsString('inDegrees').strip().lower()=='yes'
    cs=working_model.getCoordinateSet();angles=np.empty((len(times),cs.getSize()));errors=np.empty((len(times),len(names)))
    for i in range(cs.getSize()):
        c=cs.get(i)
        if c.getName() not in columns:raise RuntimeError('Missing output coordinate '+c.getName())
        col=motion.getDependentColumn(c.getName())
        angles[:,i]=np.array([col[k] for k in range(len(times))])
        if in_degrees and c.getMotionType()==o.Coordinate.Rotational:angles[:,i]=np.deg2rad(angles[:,i])
        if c.getDefaultLocked() and not np.allclose(angles[:,i],c.getDefaultValue(),atol=1e-8):raise RuntimeError('Locked coordinate changed: '+c.getName())
    if not np.isfinite(angles).all():raise RuntimeError('Nonfinite IK coordinates.')
    for frame in range(len(times)):
        for i in range(cs.getSize()):cs.get(i).setValue(s,float(angles[frame,i]),False)
        errors[frame]=np.linalg.norm(positions(working_model,s,names)-values[frame],axis=1)*1000
    rms=np.sqrt(np.mean(errors**2,axis=1));maximum=errors.max(axis=1);worst=np.unravel_index(errors.argmax(),errors.shape)
    with (out/'marker_errors_mm.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['time','RMS_mm','Max_mm','worst_marker']+names)
        for i,time in enumerate(times):writer.writerow([time,rms[i],maximum[i],names[int(errors[i].argmax())]]+errors[i].tolist())
    per_marker=[dict(marker=n,rms_mm=float(np.sqrt(np.mean(errors[:,j]**2))),max_mm=float(errors[:,j].max())) for j,n in enumerate(names)]
    ranges=[];jump_flags=[]
    for i in range(cs.getSize()):
        c=cs.get(i);rotational=c.getMotionType()==o.Coordinate.Rotational;v=np.rad2deg(angles[:,i]) if rotational else angles[:,i]
        ranges.append(dict(coordinate=c.getName(),unit='deg' if rotational else 'm',min=float(v.min()),max=float(v.max()),locked=c.getDefaultLocked()))
        if rotational and not c.getDefaultLocked() and np.max(abs(np.diff(v)))>20:jump_flags.append(c.getName())
    assert hashlib.sha256(model_path.read_bytes()).hexdigest()==model_hash
    metrics=dict(subject=args.subject,trial=c3d.name,source_model=str(model_path),c3d=str(c3d),source_sha256=model_hash,c3d_sha256=hashlib.sha256(c3d.read_bytes()).hexdigest(),frames=len(times),start=float(times[0]),end=float(times[-1]),marker_count=len(names),weights=weights,rotation=rotation.tolist(),locked_coordinates=locked,observability=observability,global_rms_mm=float(np.sqrt(np.mean(errors**2))),mean_frame_rms_mm=float(rms.mean()),peak_frame_rms_mm=float(rms.max()),max_marker_error_mm=float(errors.max()),worst_marker=names[worst[1]],worst_time=float(times[worst[0]]),frames_rms_above_20=int((rms>20).sum()),frames_max_above_40=int((maximum>40).sum()),per_marker=per_marker,coordinate_ranges=ranges,jump_screen=jump_flags,missing_geometry=missing,alignment=alignment,opensim_version=o.GetVersion(),
                 execution_policy=dict(requested_start=args.start,requested_end=args.end,requested_weights=json.loads(Path(args.weights).read_text(encoding='utf-8-sig')) if args.weights else {},marker_map=mapping or {},accuracy=1e-5,locking_policy='zero_sensitivity_profile_v1'),
                 timing=timing,implementation_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(__file__).with_name('marker_alignment.py'),Path(__file__).with_name('timing.py'))})
    (out/'metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2),encoding='utf-8')
    locks='\n'.join(f"- `{r['coordinate']}`：{r['value_deg']:.3f}°。" for r in locked)
    table='\n'.join(f"| {r['marker']} | {r['rms_mm']:.2f} | {r['max_mm']:.2f} |" for r in sorted(per_marker,key=lambda r:r['max_mm'],reverse=True))
    text=f'''# {args.subject} / {c3d.stem} IK 报告

源模型：`{model_path}`；动态C3D：`{c3d}`。工作模型 `ik_model.osim`；保留原膝单旋转定义，未调整标记位置或骨架尺度。源模型SHA256未改变。

## 无标记观测的坐标锁定

按用户要求，在工作副本锁定以下坐标，采用原模型默认值（均是建模假设，不是本次实验测量）：

{locks}

使用三个姿态（默认、膝屈伸坐标-0.3/-0.8 rad并配合髋屈曲0.2rad）的标记位置有限差分核查。剩余{observability['active_coordinates']}个可动坐标的数值秩为{observability['ranks']}。距骨虽无直接标记，下游足部标记可约束踝关节，因此保留踝及距下关节自由度。该数值检查不等于全程或解剖验证。

## 数据和配置

{len(times)}帧，{times[0]:.3f}–{times[-1]:.3f}s，跟踪{len(names)}个模型点；实际权重见metrics.json，未配置的点为1。mm→m，旋转矩阵{rotation.tolist()}。本次不做滤波、插值、模型点移动或自动调整权重。额外实验点及点名对应见 `点名称对齐报告.md`；旧模板中未加入模型的簇标记未参与任务。

原始C3D帧号：{timing['selected_first_frame']}–{timing['selected_last_frame']}；所选窗口距文件首样本偏移{timing['selected_start_offset_s']:.6f}s。逐帧对应见time_mapping.csv，完整时间基准见timing.json。保留adapter时间，裁剪窗口不重新归零；帧号时间不是墙钟，也不能单独证明与外部测力/视频同步。

## 误差结果

- 全部帧/点整体RMS：**{metrics['global_rms_mm']:.2f} mm**；逐帧RMS均值：{metrics['mean_frame_rms_mm']:.2f} mm；逐帧RMS峰值：{metrics['peak_frame_rms_mm']:.2f} mm。
- 全程最大点误差：**{metrics['max_marker_error_mm']:.2f} mm**，{metrics['worst_marker']}，t={metrics['worst_time']:.3f}s。
- RMS>20mm的帧数：{metrics['frames_rms_above_20']}/{len(times)}；Max>40mm的帧数：{metrics['frames_max_above_40']}/{len(times)}。
- 指标由保存的ik.mot重新计算模型点位置，独立与实验坐标比较；详见 `marker_errors_mm.csv`。这些指标不是静态缩放误差。
- 邻帧旋转变化超过20°的粗筛标记：{jump_flags}。20°是工程筛查值，不是生理标准；坐标范围见metrics.json，数值合理性仍需人工核对。

| 点名 | 全程RMS mm | 全程Max mm |
|---|---:|---:|
{table}

[OpenSim IK最佳实践](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim33/pages/53674534)给出的RMS<20mm、Max通常<20–40mm是依模型/动作而变的经验建议，不是模型验收标准。锁定的跖趾和腰椎角不得报告为实验测得结果。尚未进行人工运动可视化，也没有运行ID/RRA/SO/CMC/JRA。

输出：ik.mot、ik_setup.xml、ik_model.osim、markers.trc、OpenSim原生误差/模型点位置STO、独立误差CSV、metrics.json及ik_run.log。显示几何缺失：{missing}。
'''
    (out/'IK报告.md').write_text(text,encoding='utf-8')
    from review_ik import review_trial
    review=review_trial(out,visualize=not args.no_visuals,max_visual_frames=args.max_visual_frames,analysis_window=analysis_window)
    with (out/'IK报告.md').open('a',encoding='utf-8') as f:
        f.write('\n## 自动核查\n\n状态：'+review['classification']+'。见核查报告.md、review.json、异常时间段.csv；该状态不是解剖或动态学验收结论。\n')
    print(json.dumps({k:metrics[k] for k in ('trial','frames','global_rms_mm','max_marker_error_mm','worst_marker')},ensure_ascii=False))


if __name__=='__main__':main()
