"""Deterministic diagnostic marker overlays; no model/data corrections."""
import csv
import json
import math
from pathlib import Path


def render_overlays(directory,review,max_frames=3,include_worst=False):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import opensim as o
    p=Path(directory);d=json.loads((p/'metrics.json').read_text(encoding='utf-8'))
    with (p/'marker_errors_mm.csv').open(encoding='utf-8') as f:errors=list(csv.DictReader(f))
    times=np.array([float(r['time']) for r in errors]);scores=np.array([float(r['Max_mm']) for r in errors])
    candidates=[]
    for segment in review.get('abnormal_segments',[]):
        idx=np.flatnonzero((times>=segment['start_s']-1e-8)&(times<=segment['end_s']+1e-8))
        key='RMS_mm' if segment['kind']=='frame_rms' else segment['marker']
        if len(idx):candidates.append(int(idx[np.argmax([float(errors[j].get(key,errors[j]['Max_mm'])) for j in idx])]))
    for change in review.get('angular_change_flags',[]):candidates.append(int(np.argmin(abs(times-change['end_s']))))
    if include_worst or (review.get('classification')=='needs_review' and not candidates):candidates.append(int(scores.argmax()))
    selected=sorted(set(candidates),key=lambda j:(-scores[j],j))[:max_frames]
    if not selected:return []
    o.ModelVisualizer.addDirToGeometrySearchPaths(str(p/'Geometry'));model=o.Model(str(p/'ik_model.osim'));state=model.initSystem();cs=model.getCoordinateSet()
    motion=o.TimeSeriesTable(str(p/'ik.mot'));trc=o.TimeSeriesTableVec3(str(p/'markers.trc'));names=list(d['weights']);degrees=motion.getTableMetaDataAsString('inDegrees').lower()=='yes'
    mt=np.array(list(motion.getIndependentColumn()));tt=np.array(list(trc.getIndependentColumn()))
    if len(mt)!=len(times) or len(tt)!=len(times) or not np.allclose(mt,times,atol=1e-8) or not np.allclose(tt,times,atol=1e-8):raise ValueError('Overlay time vectors disagree')
    timing=[]
    if review.get('timing') and (p/'time_mapping.csv').is_file():
        with (p/'time_mapping.csv').open(encoding='utf-8') as f:timing=list(csv.DictReader(f))
    dest=p/'异常帧可视化';dest.mkdir(exist_ok=True);outputs=[];md=['# Marker overlay review','', 'Blue circles: experimental; orange crosses: model. Lines show residuals. Top: all tracked points; bottom: local region around largest residual. Axes use model coordinates in mm. No bones/surfaces or automatic corrections are inferred.','']
    projections=[(0,1,'X','Y'),(2,1,'Z','Y'),(0,2,'X','Z')]
    for j in selected:
        for k in range(cs.getSize()):
            c=cs.get(k);value=float(motion.getDependentColumn(c.getName())[j])
            if degrees and c.getMotionType()==o.Coordinate.Rotational:value=math.radians(value)
            c.setValue(state,value,False)
        model.realizePosition(state)
        predicted=np.array([[model.getMarkerSet().get(n).getLocationInGround(state)[k]*1000 for k in range(3)] for n in names])
        observed=np.array([[trc.getDependentColumn(n)[j][k]*1000 for k in range(3)] for n in names])
        residual=np.linalg.norm(predicted-observed,axis=1);worst=int(residual.argmax())
        if not np.allclose(residual,[float(errors[j][n]) for n in names],atol=.002):raise ValueError('Overlay reconstruction disagrees with independent error CSV')
        center=(observed[worst]+predicted[worst])/2
        local=np.flatnonzero(np.linalg.norm(observed-center,axis=1)<300);local=np.unique(np.r_[local,worst])
        frame=int(timing[j]['source_frame']) if timing else None
        ordinal=int(timing[j]['source_sample_index_0'])+1 if timing else None
        full=review['timing']['source_frames'] if timing else None
        fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
        for row,subset in enumerate([np.arange(len(names)),local]):
            for ax,(a,b,la,lb) in zip(axes[row],projections):
                ax.scatter(observed[subset,a],observed[subset,b],s=26,color='#2374AB',label='Experimental')
                ax.scatter(predicted[subset,a],predicted[subset,b],s=34,marker='x',color='#D56A19',label='Model')
                for k in subset:
                    ax.plot([observed[k,a],predicted[k,a]],[observed[k,b],predicted[k,b]],color='#888888',lw=.8,alpha=.8)
                    if row==1 or residual[k]>40:ax.annotate(names[k],observed[k,[a,b]],fontsize=7,xytext=(3,3),textcoords='offset points')
                ax.set(xlabel=la+' (mm)',ylabel=lb+' (mm)');ax.margins(.12);ax.set_aspect('equal',adjustable='datalim');ax.grid(alpha=.18)
        axes[0,0].legend(fontsize=8)
        fig.suptitle(f"{d['trial']} | IK t={times[j]:.3f} s | trial sample {ordinal}/{full} | source frame={frame if frame is not None else 'unknown'}\nRMS={float(errors[j]['RMS_mm']):.2f} mm | Max={residual[worst]:.2f} mm ({names[worst]})",fontsize=12)
        filename=f'frame_{frame if frame is not None else j+1}_t_{times[j]:.3f}.png';fig.savefig(dest/filename,dpi=135);plt.close(fig)
        outputs.append(dict(file=str((dest/filename).relative_to(p)),ik_time_s=float(times[j]),source_frame=frame,trial_frame=ordinal,trial_total_frames=full,worst_marker=names[worst],max_mm=float(residual[worst])))
        md.extend([f"## Trial sample {ordinal}/{full} | t={times[j]:.3f} s | source frame {frame}",'',f'![Marker overlay]({filename})',''])
    if len(set(candidates))>len(selected):md.append(f'Output capped at {max_frames} frames; remaining episodes are listed in ../异常时间段.csv.')
    (dest/'index.md').write_text('\n'.join(md),encoding='utf-8');return outputs


def render_timeline(directory,review):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p=Path(directory)
    with (p/'marker_errors_mm.csv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
    position=review['trial_position'];full=position['trial_total_frames'];mapped=full is not None
    if mapped:
        with (p/'time_mapping.csv').open(encoding='utf-8') as f:mapping=list(csv.DictReader(f))
        x=np.array([int(v['source_sample_index_0'])+1 for v in mapping])
    else:x=np.arange(1,len(rows)+1);full=len(rows)
    fig,ax=plt.subplots(figsize=(11,3.5),layout='constrained')
    ax.plot(x,[float(v['Max_mm']) for v in rows],color='#D56A19',lw=1.1,label='Frame Max (mm)')
    ax.plot(x,[float(v['RMS_mm']) for v in rows],color='#2374AB',lw=1.1,label='Frame RMS (mm)')
    ax.axhline(review['thresholds']['point_max_mm'],color='#D56A19',ls='--',lw=.8)
    ax.axhline(review['thresholds']['frame_rms_mm'],color='#2374AB',ls='--',lw=.8)
    for k,r in enumerate(position['abnormal_ranges']):ax.axvspan(r['start_frame']-.5,r['end_frame']+.5,color='#D56A19',alpha=.12,label='Flagged samples' if k==0 else None)
    context=review['analysis_context']
    if context['status']=='configured':
        w=context['window'];ax.axvspan(w['start_frame']-.5,w['end_frame']+.5,color='#368A55',alpha=.12,label='Configured analysis interval')
    note='Analysis interval not configured; contact phase unknown' if context['status']=='unconfigured' else 'Analysis interval explicitly configured; full-trial errors retained'
    ax.set(xlim=(.5,full+.5),xlabel=f"{'Full source trial sample' if mapped else 'IK row (original mapping unavailable)'} (1-based; total {full})",ylabel='Marker error (mm)',title=review['trial']+'\n'+note)
    ax.grid(alpha=.2);ax.legend(fontsize=8,loc='upper right');dest=p/'异常帧可视化';dest.mkdir(exist_ok=True)
    fig.savefig(dest/'error_timeline.png',dpi=135);plt.close(fig)
    return '异常帧可视化/error_timeline.png'
