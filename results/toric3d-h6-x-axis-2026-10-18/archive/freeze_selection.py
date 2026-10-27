"""Freeze the prospective axis scan after verification and excluded timing pilots."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import tomllib

from supervise import JULIA, environment, verify_frozen

ARCHIVE=Path(__file__).resolve().parent
STUDY=ARCHIVE.parent
REFERENCE=STUDY.parent/'toric3d-h6-matched-2026-10-17'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output=STUDY/'data/selection.json'
    if output.exists():
        raise RuntimeError('The selection already exists; do not overwrite a frozen experiment.')
    if list((ARCHIVE/'runs').glob('*/checkpoint.toml')):
        raise RuntimeError('Cannot freeze after production sampling has started.')
    old=json.loads((REFERENCE/'data/selection.json').read_text())
    report=json.loads((REFERENCE/'data/report.json').read_text())
    if not report['complete'] or report['output_sha256']['data/matched-results.csv'] != sha(REFERENCE/'data/matched-results.csv'):
        raise RuntimeError('Diagonal reference is incomplete or its CSV differs from the verified report.')
    verification=json.loads((ARCHIVE/'verification.json').read_text())
    if verification.get('status') != 'passed':
        raise RuntimeError('Independent verification must pass before freezing.')
    for relative,claimed in verification['tested_source_hashes'].items():
        if sha(ARCHIVE/relative) != claimed:
            raise RuntimeError(f'Source changed after verification: {relative}')
    pilot=tomllib.loads((ARCHIVE/'timing-pilot.toml').read_text())
    if not pilot['excluded_from_production'] or len(pilot['results']) != 18:
        raise RuntimeError('The excluded timing pilot is incomplete.')
    channels={ch:{'sizes':old['channels'][ch]['sizes'],'rates':old['channels'][ch]['rates']}
              for ch in ('x_only','z_only')}
    shards={}
    for ch,spec in channels.items():
        for L in spec['sizes']:
            name=f'{ch}-L{L}'
            seed=int(hashlib.sha256(f'physical-x-h6-20261018/{name}'.encode()).hexdigest()[:15],16)
            shards[name]={'channel':ch,'size':L,'rates':spec['rates'],'seed':seed}
    lines=[]
    for name,shard in shards.items():
        lines.extend((f'[shards."{name}"]', f'seed = {shard["seed"]}',
                      f'count = {len(shard["rates"])}',''))
    input_path=ARCHIVE/'selection-input.toml'
    metadata_path=ARCHIVE/'prospective-metadata.toml'
    input_path.write_text('\n'.join(lines))
    code='''
include(ARGS[1])
using .H6AxisStudy
import TOML
s=TOML.parsefile(ARGS[2])
meta=Dict("source_hashes"=>H6AxisStudy.source_hashes(),
    "runtime"=>H6AxisStudy.runtime_metadata(),
    "protocols"=>Dict(c=>H6AxisStudy.protocol_metadata(c) for c in ("x_only","z_only")),
    "point_seeds"=>Dict(name=>H6AxisStudy.hexseed.(H6AxisStudy.point_seeds(shard["seed"],shard["count"]))
        for (name,shard) in s["shards"]))
open(ARGS[3],"w") do io
    TOML.print(io,meta;sorted=true)
end
'''
    subprocess.run([JULIA,'--compiled-modules=existing',f'--project={ARCHIVE/"runtime"}',
                    '--startup-file=no','-e',code,str(ARCHIVE/'run_study.jl'),
                    str(input_path),str(metadata_path)],env=environment(),check=True)
    meta=tomllib.loads(metadata_path.read_text())
    for name,shard in shards.items():
        shard['point_seeds_hex']=meta['point_seeds'][name]
    total_points=sum(len(spec['sizes'])*len(spec['rates']) for spec in channels.values())
    auxiliary=['analyze.py','supervise.py','paper-no-grid.mplstyle','freeze_selection.py',
               'timing_pilot.jl','test_runner.jl','test_supervisor.py','test_analysis.py',
               'test_axis_geometry.jl','test_axis_feedback.jl']
    selection={
        'format_version':1,'stage':'direction_comparison','frozen_before_sampling':True,
        'frozen_at_unix':time.time(),'shots':1000,'batch_size':25,
        'arms':['final_only','one_step','two_slice'],'channels':channels,'shards':shards,
        'source_hashes':meta['source_hashes'],'expected_runtime':meta['runtime'],'protocols':meta['protocols'],
        'auxiliary_source_hashes':{p:sha(ARCHIVE/p) for p in auxiliary},
        'reference':{kind:str(REFERENCE/'data'/file) for kind,file in
                     (('csv','matched-results.csv'),('report','report.json'))},
        'common_z_seed_repair':True,'pilot_excluded':True,
        'total_fault_records':total_points*1000,'total_protocol_outcomes':total_points*3000,
        'verification_sha256':sha(ARCHIVE/'verification.json'),
        'timing_pilot_sha256':sha(ARCHIVE/'timing-pilot.toml'),
        'purpose':'Exploratory physical-x versus diagonal growth comparison; no axis threshold fit is prespecified.',
        'sampling_pairing':'Shared across the three arms within each direction; independent across directions.',
        'geometry_limit':'Same L matches transfer-cell count, not seed geometry, seed exposure, or dZ.',
        'one_step_z_limit':'Physical-x uses a unique forward-edge defect push; diagonal uses matched endpoint routing.'}
    for kind in ('csv','report'):
        selection['reference'][kind+'_sha256']=sha(Path(selection['reference'][kind]))
    verify_frozen(STUDY,selection)
    output.write_text(json.dumps(selection,indent=2)+'\n')
    print(json.dumps({'selection':str(output),'points':total_points,'shared_records':total_points*1000}))


if __name__=='__main__':
    main()
