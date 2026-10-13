"""Independent saved-artifact audit; run from repository root with a stage path."""
import csv, hashlib, itertools, json, math, sys, tomllib
from pathlib import Path
P=('none','one_step','window_complete','full_prefix')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def near(a,b): assert math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-12),(a,b)
def wilson(k,n,z=1.959963984540054):
    p=k/n;den=1+z*z/n;c=(p+z*z/(2*n))/den;r=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0,c-r),min(1,c+r)
b=Path(sys.argv[1]);stem='toric3d-h5-vertex-projector-feedback';raw=b/(stem+'-raw.csv');mp=b/(stem+'-metadata.toml')
m=tomllib.loads(mp.read_text());ck=tomllib.loads((b/'checkpoints/checkpoint.toml').read_text());warm=tomllib.loads((b/'warmup-timing.toml').read_text())
assert m['complete'] and m['invalid_recoveries']==0 and m['raw_csv_sha256']==sha(raw)
assert m['protocols']==list(P) and m['backend_versions']=={'ldpc':'2.4.1','pymatching':'2.4.0'}
assert m['config']==ck['config'] and warm['zero_noise_passed'] and warm['excluded_from_scientific_counts']
assert all(sha(p)==h for p,h in m['source_hashes'].items())
expected={(ch,*dims,float(p)) for ch,g in m['config']['channel_grids'].items() for dims in g['dimensions'] for p in g['rates']}
groups={}
for r in csv.DictReader(raw.open()):
    key=(r['channel'],int(r['T']),int(r['Lu']),int(r['Lv']),float(r['p']));g=groups.setdefault(key,{})
    assert r['metric'] not in g;g[r['metric']]=r
    n=int(r['shots']);c=[int(r[f'joint_mask_{i}']) for i in range(16)]
    assert n==m['shots'] and sum(c)==n and min(c)>=0 and int(r['invalid_recoveries'])==0 and r['complete']=='true'
    near(r['p_eff'],min(1-1e-12,max(1e-12,float(r['p']))))
    for i,pr in enumerate(P):
        k=sum(v for mask,v in enumerate(c) if mask&(1<<i));assert k==int(r[pr+'_failures']);near(r[pr+'_rate'],k/n)
        lo,hi=wilson(k,n);near(r[pr+'_ci_low'],lo);near(r[pr+'_ci_high'],hi)
    for (i,a),(j,z) in itertools.combinations(enumerate(P),2):
        cells=[sum(v for mask,v in enumerate(c) if (bool(mask&(1<<i))+2*bool(mask&(1<<j)))==s) for s in range(4)]
        pr=a+'_vs_'+z
        for suffix,v in zip(('neither','first_only','second_only','both'),cells): assert int(r[pr+'_'+suffix])==v
        d=(cells[2]-cells[1])/n;near(r[pr+'_difference'],d);near(r[pr+'_standard_error'],math.sqrt(max(0,(cells[1]+cells[2])/n-d*d)/(n-1)))
        al,ah=wilson(cells[1],n,2.241402727604947);bl,bh=wilson(cells[2],n,2.241402727604947)
        near(r[pr+'_ci_low'],max(-1,bl-ah));near(r[pr+'_ci_high'],min(1,bh-al))
assert set(groups)==expected and m['completed_points']==m['expected_points']==len(expected)
for g in groups.values():
    assert len(g)==23
    for p in P: assert sum(int(g['pauli_class_'+str(i)][p+'_failures']) for i in range(16))==m['shots']
metrics=['logical_x_1','logical_x_2','logical_x_union','logical_z_1','logical_z_2','logical_z_union','pauli_any']+['pauli_class_'+str(i) for i in range(16)]
ordered=[(ch,*d,float(p)) for ch in ('x_only','z_only') if ch in m['config']['channel_grids'] for d in m['config']['channel_grids'][ch]['dimensions'] for p in m['config']['channel_grids'][ch]['rates']]
assert len(ck['states'])==len(ordered)==len(ck['point_seeds'])
for k,state,seed in zip(ordered,ck['states'],ck['point_seeds']):
    assert state['completed_shots']==m['shots'] and state['invalid_recoveries']==0
    for metric,counts in zip(metrics,state['counts']):
        r=groups[k][metric]
        assert counts==[int(r['joint_mask_'+str(i)]) for i in range(16)]
        assert r['point_seed']==seed
        for name in ('record_digest','signature_digest'): assert r[name]==state[name]
        for name in ('sampling_seconds','replay_seconds','decoder_build_seconds'): near(r[name],state[name])
assert not list(b.rglob('failure.toml'))
primary=[g['logical_x_union' if k[0]=='x_only' else 'logical_z_union'] for k,g in groups.items()]
summary=dict(stage=m['stage'],points=len(groups),paired_records=len(groups)*m['shots'],arm_replays=4*len(groups)*m['shots'],metrics_per_point=23,contrasts_per_metric=6,invalid_recoveries=0,raw_csv_sha256=sha(raw),metadata_sha256=sha(mp),checkpoint_sha256=sha(b/'checkpoints/checkpoint.toml'),source_hashes_verified=True,statistics_verified=True,checkpoint_counts_and_digests_verified=True,warmup_excluded=True,zero_noise_passed=True,phase_seconds_by_channel={ch:sum(float(r[x]) for r in primary if r['channel']==ch for x in ('decoder_build_seconds','sampling_seconds','replay_seconds')) for ch in m['config']['channel_grids']})
(b/'independent-audit.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n');print(json.dumps(summary,indent=2))
