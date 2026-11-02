"""Shared reuse gate for preflight and batch execution."""
import hashlib
from pathlib import Path


def compatible(directory,metrics,model_hash,c3d,version,rotation,policy):
    scripts=Path(__file__).parent
    implementation={name:hashlib.sha256((scripts/name).read_bytes()).hexdigest() for name in ('run_ik.py','marker_alignment.py','timing.py')}
    required=('ik.mot','ik_model.osim','markers.trc','time_mapping.csv','timing.json','marker_errors_mm.csv')
    return (metrics.get('implementation_hashes')==implementation and metrics.get('opensim_version')==version
            and metrics.get('execution_policy')==policy and metrics.get('rotation')==[rotation[i:i+3] for i in (0,3,6)]
            and metrics.get('source_sha256')==model_hash and metrics.get('c3d_sha256')==hashlib.sha256(Path(c3d).read_bytes()).hexdigest()
            and all((Path(directory)/name).is_file() for name in required))
