"""Conservative name-only alignment: model marker names follow C3D labels."""
import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def token(name):
    name=name.rsplit(':',1)[-1]
    name=re.sub(r'^(?:Subject\d+|S\d+)[_ .-]+','',name,flags=re.I)
    return re.sub(r'[^A-Za-z0-9]','',name).upper()


ALIASES={}
for side in ('L','R'):
    for suffix,variants in {
        'ASI':['ASIS'], 'PSI':['PSIS'], 'KNE':['Knee.Lat','Knee.Lateral'],
        'KNM':['Knee.Med','Knee.Medial'], 'ANK':['Ankle.Lat','Ankle.Lateral'],
        'MED':['Ankle.Med','Ankle.Medial'], 'HEE':['Heel'], 'TOE':['Toe.Tip'],
        'FMH':['Met1.Head'], 'SMH':['Met2.Head'], 'VMH':['Met5.Head'],
    }.items():
        canonical=side+suffix
        for name in [canonical]+[side+'.'+v for v in variants]:ALIASES[token(name)]=canonical


def canonical_name(name):
    t=token(name)
    return ALIASES.get(t,t)


def anatomical_identity(source,target=None):
    identity=canonical_name(source)
    if identity in ALIASES.values() or target is None:return identity
    return canonical_name(target)


def derived_candidate(name):
    """Screen only; C3D point channels are not necessarily physical markers."""
    raw=name.rsplit(':',1)[-1]
    if raw.startswith('*') or raw.lower().endswith('_m'):return True
    if re.search(r'angles?|forces?|moments?|powers?',raw,re.I):return True
    if re.search(r'(?:HJC|KJC|AJC)$',raw,re.I):return True
    return bool(re.fullmatch(r'(?:PELVIS|[LR](?:FEMUR|TIBIA|FOOT|TOES))[OALP]',raw,re.I))


def plan(model_names,labels,explicit=None):
    explicit=explicit or {};mapping={};rows=[];errors=[]
    if not isinstance(explicit,dict) or not all(isinstance(k,str) and isinstance(v,str) for k,v in explicit.items()):
        raise ValueError('marker_name_map must be an object of original model name -> exact C3D label.')
    for key in explicit:
        if key not in model_names:errors.append('Explicit mapping source absent from model: '+key)
    if len(set(model_names))!=len(model_names):errors.append('Duplicate model marker names.')
    canonical=[canonical_name(n) for n in model_names]
    if len(set(canonical))!=len(canonical):errors.append('Multiple model markers share the same normalized identity; resolve manually.')
    if len(set(labels))!=len(labels):errors.append('Duplicate C3D labels.')
    for name in model_names:
        if name in explicit:
            choices=[explicit[name]] if explicit[name] in labels else []
            method='explicit'
        elif name in labels:
            choices=[name];method='exact'
        else:
            choices=[n for n in labels if canonical_name(n)==canonical_name(name)];method='normalized/alias'
        if len(choices)!=1:
            errors.append(f'{name}: '+('no C3D counterpart' if not choices else 'ambiguous: '+', '.join(choices)))
            rows.append(dict(model_name=name,vicon_name=None,method=method,candidates=choices));continue
        mapping[name]=choices[0];rows.append(dict(model_name=name,vicon_name=choices[0],method=method,candidates=choices))
    if len(set(mapping.values()))!=len(mapping):errors.append('Several model markers map to the same C3D point; resolve manually.')
    identities={source:anatomical_identity(source,target) for source,target in mapping.items()}
    if len(set(identities.values()))!=len(identities):errors.append('Mapped model markers share one anatomical identity; resolve manually.')
    unmatched=[n for n in labels if n not in mapping.values()]
    return dict(ready=not errors,mapping=mapping,identities=identities,rows=rows,errors=errors,
                vicon_without_model=[n for n in unmatched if not derived_candidate(n)],
                derived_or_unlabelled=[n for n in unmatched if derived_candidate(n)],
                renamed_count=sum(k!=v for k,v in mapping.items()))


def rename_tree(tree,mapping):
    """Preserve every location, body, and scale; update named component references."""
    for e in tree.iter('Marker'):
        old=e.get('name')
        if old in mapping:e.set('name',mapping[old])
    for e in tree.iter():
        if not e.text:continue
        text=e.text
        # Match full marker path components, never perform overlapping substring replacement.
        text=re.sub(r'(/markerset/)([^/\s]+)',lambda m:m[1]+mapping.get(m[2],m[2]),text)
        if e.tag=='members':text=' '.join(mapping.get(n,n) for n in text.split())
        e.text=text
    return tree


def report_text(p):
    escape=lambda s:str(s).replace('|','/').replace('\n',' ')
    lines=['# 模型与 Vicon 点名称对齐','','状态：'+('可继续' if p['ready'] else '需要手动处理，禁止缩放'),'',
           '方向：修改输出 OpenSim 模型的点名，使其与 C3D 标签一致；原模型、点坐标及所属身体分段保持不变。这里只判断名称对应，不能证明解剖位置相同。','',
           '| 原模型点名 | Vicon/C3D 点名 | 依据 |','|---|---|---|']
    lines += [f"| {escape(r['model_name'])} | {escape(r['vicon_name'] or '未确定')} | {r['method']} |" for r in p['rows']]
    if p['errors']:lines+=['','## 需要手动确认的冲突/缺失','']+['- '+escape(e) for e in p['errors']]
    lines+=['','## Vicon 有数据但模型中无对应点','',
            '**以下候选点如需参与 OpenSim 分析，需要你在模型中手动添加，并指定正确身体分段和局部位置。程序不会凭名称虚构位置；无需使用的额外点可以保留在实验数据中。**','']
    lines+=['- '+escape(n) for n in p['vicon_without_model']] or ['无。']
    lines+=['','## 疑似派生量、虚拟点或未标记点','',
            '以下仅按名称进行初步筛查，不自动要求添加为实测标记。请结合 Vicon 定义核实：','',
            ', '.join(escape(n) for n in p['derived_or_unlabelled']) or '无。']
    return '\n'.join(lines)+'\n'


def save_plan(p,directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    (directory/'marker_alignment.json').write_text(json.dumps(p,ensure_ascii=False,indent=2),encoding='utf-8')
    (directory/'点名称对齐报告.md').write_text(report_text(p),encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('model','c3d','output'):parser.add_argument('--'+key,required=True)
    parser.add_argument('--marker-map',help='JSON original-model-name -> exact-C3D-label mapping')
    parser.add_argument('--overwrite',action='store_true');args=parser.parse_args()
    model,c3d,out=Path(args.model).resolve(),Path(args.c3d).resolve(),Path(args.output).resolve()
    if out in (model.parent,c3d.parent) or out==model or out==c3d:raise ValueError('Use a separate output directory.')
    if out.exists() and any(out.iterdir()) and not args.overwrite:raise ValueError('Nonempty output requires --overwrite.')
    import opensim as o
    tree=ET.parse(model);a=o.C3DFileAdapter();labels=list(a.getMarkersTable(a.read(str(c3d))).getColumnLabels())
    explicit=json.loads(Path(args.marker_map).read_text(encoding='utf-8-sig')) if args.marker_map else None
    p=plan([e.get('name') for e in tree.iter('Marker')],labels,explicit);save_plan(p,out)
    if not p['ready']:raise SystemExit(2)
    rename_tree(tree,p['mapping']).write(out/'aligned_model.osim',encoding='utf-8',xml_declaration=True)
    print(f"Renamed {p['renamed_count']} model markers; {len(p['vicon_without_model'])} experimental candidate points need manual model additions if used.")


if __name__=='__main__':main()
