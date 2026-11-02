"""Read C3D frame provenance without extra dependencies; never reset selected time."""
import csv
import json
from pathlib import Path
import struct


def c3d_metadata(path):
    with Path(path).open('rb') as f:
        h=f.read(512)
        if len(h)!=512 or h[1]!=80 or not h[0]:raise ValueError('Invalid C3D header')
        f.seek((h[0]-1)*512);prefix=f.read(4)
        if len(prefix)!=4 or prefix[3] not in (84,86):raise ValueError('Timing reader supports Intel/MIPS C3D; DEC requires a supported converter')
        endian='<' if prefix[3]==84 else '>';f.seek((h[0]-1)*512);params=f.read(prefix[2]*512)
    groups={};records=[];pos=4
    while pos+2<=len(params):
        length,group=struct.unpack_from('bb',params,pos)
        if not length or not group:break
        length=abs(length);name=params[pos+2:pos+2+length].decode('ascii');at=pos+2+length
        offset=struct.unpack_from(endian+'h',params,at)[0];body=at+2
        if group<0:groups[-group]=name
        else:
            dtype,ndim=struct.unpack_from('bB',params,body);dims=list(params[body+2:body+2+ndim]);count=1
            for size in dims:count*=size
            begin=body+2+ndim;records.append((group,name,dtype,params[begin:begin+count*abs(dtype)]))
        if not offset:break
        next_pos=at+offset
        if next_pos<=pos or next_pos>len(params):raise ValueError('Malformed C3D parameter offsets')
        pos=next_pos
    values={(groups.get(g),n):(dtype,data) for g,n,dtype,data in records}
    def field(name):
        value=values.get(('TRIAL',name))
        if value is None:return None
        dtype,data=value
        if dtype!=2 or len(data)!=4:raise ValueError('Invalid TRIAL:'+name)
        lo,hi=struct.unpack(endian+'HH',data);return lo+(hi<<16)
    return dict(header_first_frame=struct.unpack_from(endian+'H',h,6)[0],header_last_frame=struct.unpack_from(endian+'H',h,8)[0],
                trial_first_frame=field('ACTUAL_START_FIELD'),trial_last_frame=field('ACTUAL_END_FIELD'),
                point_rate_hz=struct.unpack_from(endian+'f',h,20)[0],analog_samples_per_point_frame=struct.unpack_from(endian+'H',h,18)[0],processor=prefix[3])


def write_timing(directory,c3d,full_times,selected_times):
    import numpy as np
    full=np.asarray(full_times);times=np.asarray(selected_times);meta=c3d_metadata(c3d);rate=meta['point_rate_hz']
    if rate<=0 or len(full)<2 or not np.allclose(np.diff(full),1/rate,rtol=1e-5,atol=1e-8):raise ValueError('C3D header rate/adapter times disagree')
    start,end=meta['trial_first_frame'],meta['trial_last_frame'];basis='TRIAL:ACTUAL_START_FIELD'
    if start is None and end is None:start,end=meta['header_first_frame'],meta['header_last_frame'];basis='C3D header (TRIAL frame parameters absent)'
    if start is None or end is None or end-start+1!=len(full):raise ValueError('Original C3D frame range disagrees with marker row count')
    indices=np.rint((times-full[0])*rate).astype(int)
    if len(times)<2 or np.any(indices<0) or np.any(indices>=len(full)) or not np.allclose(full[indices],times,atol=1e-8):raise ValueError('IK samples do not map to original C3D samples')
    rows=[dict(ik_row=i+1,ik_time_s=float(t),c3d_adapter_time_s=float(full[j]),source_sample_index_0=int(j),source_frame=int(start+j),
               source_elapsed_s=float(j/rate),frame1_based_time_s=float((start+j-1)/rate)) for i,(t,j) in enumerate(zip(times,indices))]
    meta.update(frame_number_basis=basis,source_c3d=str(Path(c3d).resolve()),source_frames=len(full),source_first_frame=start,source_last_frame=end,selected_frames=len(rows),selected_first_frame=rows[0]['source_frame'],selected_last_frame=rows[-1]['source_frame'],
                selected_start_offset_s=rows[0]['source_elapsed_s'],adapter_start_time_s=float(full[0]),
                ik_to_frame1_based_offset_s=float((start-1)/rate-full[0]),
                time_basis_note='frame1_based_time_s = (source_frame-1)/point_rate; a frame-index convention, NOT wall clock or independently verified force/video synchronization. IK retains adapter time, without window re-zeroing. Analog samples per point frame alone do not establish external force latency.')
    p=Path(directory);p.mkdir(parents=True,exist_ok=True)
    with (p/'time_mapping.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (p/'timing.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    return meta,rows
