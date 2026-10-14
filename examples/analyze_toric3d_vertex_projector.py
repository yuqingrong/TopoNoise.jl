#!/usr/bin/env python3
"""Frozen scout selection and guarded four-arm direct-projector H+5 analysis."""
import argparse, copy, csv, hashlib, importlib.util, json, math, os, tempfile, tomllib
from pathlib import Path
import numpy as np
PROTOCOLS=('none','one_step','window_complete','full_prefix')
CHANNELS=('x_only','z_only')
_spec=importlib.util.spec_from_file_location('_h5_likelihood_engine',Path(__file__).with_name('fit_toric3d_diagonal_feedback.py'))
engine=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(engine)
# A private module instance reuses the established numerical likelihood, not the
# legacy three-arm parser or public API.
engine.PROTOCOL_BITS={p:1<<i for i,p in enumerate(PROTOCOLS)}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rounded(x):return float(format(x,'.8g'))
def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n');tmp.replace(path)

def _finite_json(value,path='$'):
    if isinstance(value,dict):
        if not all(isinstance(k,str) for k in value):raise ValueError(f'{path}: JSON object keys must be strings')
        for k,v in value.items():_finite_json(v,f'{path}.{k}')
    elif isinstance(value,list):
        for i,v in enumerate(value):_finite_json(v,f'{path}[{i}]')
    elif isinstance(value,float) and not math.isfinite(value):raise ValueError(f'{path}: nonfinite numeric value')
    elif value is not None and not isinstance(value,(str,int,float,bool)):raise ValueError(f'{path}: unsupported JSON value')

def _count(value,label):
    if isinstance(value,bool) or not isinstance(value,int) or value<0:raise ValueError(f'{label} must be a nonnegative integer')
    return value

def _logical_class(value,label):
    if not isinstance(value,list) or len(value)!=2 or not all(isinstance(v,bool) for v in value):raise ValueError(f'{label} must contain two booleans')
    return tuple(value)

def _validate_z_group(group,label):
    required=('syndrome_indices','upper_parity_odd','class_counts','selected_class','baseline_class','mixed','unavoidable_count','baseline_failure_count','control_failure_count')
    if not isinstance(group,dict) or any(k not in group for k in required):raise ValueError(f'{label}: incomplete group')
    syndrome=group['syndrome_indices']
    if not isinstance(syndrome,list) or any(isinstance(i,bool) or not isinstance(i,int) or i<1 for i in syndrome) or syndrome!=sorted(set(syndrome)):raise ValueError(f'{label}: malformed syndrome indices')
    if not isinstance(group['upper_parity_odd'],bool) or not isinstance(group['mixed'],bool):raise ValueError(f'{label}: parity/mixed flags must be boolean')
    counts=group['class_counts']
    if not isinstance(counts,list) or len(counts)!=4:raise ValueError(f'{label}: class_counts must have four entries')
    counts=[_count(v,f'{label}.class_counts') for v in counts]
    if sum(counts)==0:raise ValueError(f'{label}: empty class_counts')
    classes=((False,False),(False,True),(True,False),(True,True));selected=_logical_class(group['selected_class'],f'{label}.selected_class')
    _logical_class(group['baseline_class'],f'{label}.baseline_class')
    expected=classes[counts.index(max(counts))]
    if selected!=expected:raise ValueError(f'{label}: selected class violates deterministic tie break')
    mixed=sum(v>0 for v in counts)>1
    if group['mixed']!=mixed:raise ValueError(f'{label}: mixed flag disagrees with class counts')
    unavoidable=_count(group['unavoidable_count'],f'{label}.unavoidable_count')
    if unavoidable!=sum(counts)-max(counts):raise ValueError(f'{label}: unavoidable count disagrees with class counts')
    baseline=_count(group['baseline_failure_count'],f'{label}.baseline_failure_count')
    control=_count(group['control_failure_count'],f'{label}.control_failure_count')
    if baseline>sum(counts) or control>sum(counts):raise ValueError(f'{label}: failure count exceeds mechanisms')
    return dict(mechanisms=sum(counts),mixed=int(mixed),unavoidable=unavoidable,baseline=baseline,control=control,odd=group['upper_parity_odd'])

def _validate_z_summary(summary,label):
    keys=('mechanism_count','group_count','mixed_group_count','unavoidable_count','baseline_failure_count','control_failure_count')
    if not isinstance(summary,dict) or any(k not in summary for k in keys):raise ValueError(f'{label}: incomplete parity summary')
    return {k:_count(summary[k],f'{label}.{k}') for k in keys}

def load_z_boundary_audit(path):
    try:data=json.loads(Path(path).read_text())
    except (OSError,json.JSONDecodeError) as error:raise ValueError(f'invalid Z-boundary audit: {error}') from error
    _finite_json(data)
    if not isinstance(data,dict) or data.get('format')!='h5-z-boundary-ambiguity-v1':raise ValueError('unknown Z-boundary audit format')
    if data.get('status')!='diagnostic_only' or data.get('protocols')!=list(PROTOCOLS):raise ValueError('invalid Z-boundary audit status/protocols')
    hashes=data.get('source_hashes')
    if not isinstance(hashes,dict) or not hashes or any(not isinstance(k,str) or not isinstance(v,str) or len(v)!=64 or any(c not in '0123456789abcdef' for c in v) for k,v in hashes.items()):raise ValueError('invalid Z-boundary source hashes')
    audits=data.get('audits')
    if not isinstance(audits,list) or not audits:raise ValueError('Z-boundary audit list is empty')
    identities=set()
    total_keys=('mechanism_count','group_count','mixed_group_count','unavoidable_count','baseline_failure_count','control_failure_count')
    for index,audit in enumerate(audits):
        label=f'audits[{index}]'
        if not isinstance(audit,dict):raise ValueError(f'{label}: audit must be an object')
        dimensions=audit.get('dimensions')
        if not isinstance(dimensions,list) or len(dimensions)!=3 or any(isinstance(v,bool) or not isinstance(v,int) or v<1 for v in dimensions) or min(dimensions[1:])<2:raise ValueError(f'{label}: invalid dimensions')
        protocol=audit.get('protocol');final=audit.get('final_layer_only')
        if protocol not in PROTOCOLS or not isinstance(final,bool):raise ValueError(f'{label}: invalid protocol/scope')
        identity=(tuple(dimensions),protocol,final)
        if identity in identities:raise ValueError(f'{label}: duplicate audit identity')
        identities.add(identity)
        totals={k:_count(audit.get(k),f'{label}.{k}') for k in total_keys}
        if totals['mechanism_count']==0:raise ValueError(f'{label}: empty mechanism set')
        active=_count(audit.get('active_output_fault_count'),f'{label}.active_output_fault_count')
        if active!=totals['mechanism_count']:raise ValueError(f'{label}: active-output count mismatch')
        rate=audit.get('unavoidable_rate')
        if isinstance(rate,bool) or not isinstance(rate,(int,float)) or not math.isfinite(rate) or not 0<=rate<=1 or not math.isclose(rate,totals['unavoidable_count']/totals['mechanism_count'],abs_tol=1e-12):raise ValueError(f'{label}: invalid unavoidable rate')
        groups=audit.get('groups')
        if not isinstance(groups,list) or len(groups)!=totals['group_count']:raise ValueError(f'{label}: group count mismatch')
        stats=[_validate_z_group(group,f'{label}.groups[{i}]') for i,group in enumerate(groups)]
        aggregate=dict(mechanism_count=sum(s['mechanisms'] for s in stats),group_count=len(stats),mixed_group_count=sum(s['mixed'] for s in stats),unavoidable_count=sum(s['unavoidable'] for s in stats),baseline_failure_count=sum(s['baseline'] for s in stats),control_failure_count=sum(s['control'] for s in stats))
        if aggregate!=totals:raise ValueError(f'{label}: group totals disagree with audit totals')
        parity=audit.get('parity')
        if not isinstance(parity,dict) or set(parity)!=set(('even','odd')):raise ValueError(f'{label}: invalid parity summaries')
        even=_validate_z_summary(parity['even'],f'{label}.parity.even');odd=_validate_z_summary(parity['odd'],f'{label}.parity.odd')
        for key in total_keys:
            if even[key]+odd[key]!=totals[key]:raise ValueError(f'{label}: parity totals disagree for {key}')
        mixed=audit.get('mixed_groups')
        if not isinstance(mixed,list) or len(mixed)!=totals['mixed_group_count']:raise ValueError(f'{label}: mixed-group list mismatch')
    return data

def z_boundary_guard(audit,selected_sizes,protocol,source_hashes=None):
    if protocol not in PROTOCOLS:raise ValueError('unknown H+5 protocol')
    issues=[]
    if source_hashes is not None and audit.get('source_hashes')!=source_hashes:issues.append('Z boundary audit source hashes do not match production data')
    entries=audit.get('audits',[])
    for size in sorted(set(int(v) for v in selected_sizes)):
        matches=[a for a in entries if a.get('protocol')==protocol and a.get('final_layer_only') is True and isinstance(a.get('dimensions'),list) and a['dimensions'][0]==size]
        if len(matches)!=1:
            issues.append(f'Z boundary audit missing or ambiguous for size={size}, protocol={protocol}')
            continue
        entry=matches[0];geometry=tuple(entry['dimensions'])
        for key,label in (('mixed_group_count','mixed groups'),('unavoidable_count','unavoidable'),('control_failure_count','control failures')):
            count=entry.get(key)
            if count:
                issues.append(f'Z boundary ambiguity at dimensions={geometry}, protocol={protocol}: {label}={count}')
    return issues
def normalize_row(source):
    r=dict(source)
    for k in ('T','Lu','Lv','shots','invalid_recoveries'):r[k]=int(r[k])
    r['steps']=r['T'];r['p']=float(r['p'])
    if r['T']<1 or min(r['Lu'],r['Lv'])<2 or r['shots']<1 or not 0<=r['p']<=1:raise ValueError('invalid point geometry/probability')
    if r['channel'] not in CHANNELS or r['invalid_recoveries'] or str(r['complete']).lower()!='true':raise ValueError('incomplete/invalid point')
    counts=[int(r[f'joint_mask_{m}']) for m in range(16)]
    if min(counts)<0 or sum(counts)!=r['shots']:raise ValueError('invalid sixteen-mask counts')
    r['joint_counts']=counts;r['marginals']={};r['pairs']={}
    for i,p in enumerate(PROTOCOLS):
        k=sum(n for m,n in enumerate(counts) if m&(1<<i))
        if int(r[f'{p}_failures'])!=k:raise ValueError('marginal disagrees with joint masks')
        r[f'{p}_failures']=k;r['marginals'][p]=k
    for i,a in enumerate(PROTOCOLS):
        for j,b in enumerate(PROTOCOLS[i+1:],i+1):
            cells=[0]*4
            for m,n in enumerate(counts):cells[bool(m&(1<<i))+2*bool(m&(1<<j))]+=n
            delta=(cells[2]-cells[1])/r['shots'];key=f'{a}_vs_{b}'
            aci=engine.wilson_interval(cells[1],r['shots'],z=2.241402727604947);bci=engine.wilson_interval(cells[2],r['shots'],z=2.241402727604947)
            r['pairs'][key]={'cells':cells,'difference':delta,'ci':[max(-1.,bci[0]-aci[1]),min(1.,bci[1]-aci[0])]}
            if key+'_difference' in r and not math.isclose(float(r[key+'_difference']),delta,abs_tol=1e-12):raise ValueError('paired contrast disagrees with masks')
            for suffix,val in zip(('neither','first_only','second_only','both'),cells):
                if key+'_'+suffix in r and int(r[key+'_'+suffix])!=val:raise ValueError('pair table disagrees with masks')
    for k in ('decoder_build_seconds','sampling_seconds','replay_seconds'):
        r[k]=float(r[k])
        if not math.isfinite(r[k]) or r[k]<0:raise ValueError('invalid timing')
    return r

def read_points(path):
    primary=[];allrows=[]
    with open(path,newline='') as f:
        for source in csv.DictReader(f):
            r=normalize_row(source);allrows.append(r)
            if r['metric']==('logical_x_union' if r['channel']=='x_only' else 'logical_z_union'):primary.append(r)
    if not primary:raise ValueError('no primary H+5 rows')
    # Every saved Pauli class must partition every protocol's shots.
    groups={}
    for r in allrows:groups.setdefault((r['channel'],r['T'],r['Lu'],r['Lv'],r['p']),{})[r['metric']]=r
    for group in groups.values():
        if len(group)!=23:raise ValueError('missing logical/Pauli metrics')
        for p in PROTOCOLS:
            if sum(group[f'pauli_class_{i}'][f'{p}_failures'] for i in range(16))!=next(iter(group.values()))['shots']:raise ValueError('Pauli classes do not partition shots')
    return primary

def key(r):return r['channel'],r['T'],r['Lu'],r['Lv'],r['p']
def merge_points(shards):
    merged={}
    for rows in shards:
        for r in rows:
            k=key(r)
            if k in merged and merged[k]!=r:raise ValueError('conflicting duplicate scout point')
            merged[k]=r
    return [merged[k] for k in sorted(merged)]

def validate_ancestry(inputs,metas,previous=None,parent_sha=None):
    if previous is None:
        if any(i['stage']!='scout' or i['extension_round']!=0 for i in inputs):raise ValueError('extension needs previous selection')
        return 0
    old=previous['scout_inputs']
    if inputs[:len(old)]!=old or len(inputs)!=len(old)+1:raise ValueError('reselection requires unchanged cumulative inputs plus one extension')
    last=inputs[-1];meta=metas[-1];round_=previous['extension_round']+1
    if last['stage']!='extend-scout' or last['extension_round']!=round_ or last['parent_selection_sha256']!=parent_sha:raise ValueError('extension parent/round mismatch')
    if meta.get('parent_scout_inputs')!=old or meta.get('parent_selection_sha256')!=parent_sha:raise ValueError('extension ancestry differs')
    return round_

def load_scout(csvs,metafiles,previous_path=None):
    if len(csvs)!=len(metafiles) or not csvs:raise ValueError('paired CSV/metadata inputs required')
    inputs=[];metas=[];shards=[];baseline=None
    for cp,mp in zip(csvs,metafiles):
        with open(mp,'rb') as f:meta=tomllib.load(f)
        if meta['stage'] not in ('scout','extend-scout') or not meta['complete'] or meta['invalid_recoveries']:raise ValueError('only complete valid scout shards select grids')
        if meta['raw_csv_sha256']!=sha(cp):raise ValueError('raw CSV hash mismatch')
        identity={k:meta[k] for k in ('seed','shots','decoder_settings','source_hashes','protocols','format_version','backend_versions')}
        if identity['protocols']!=list(PROTOCOLS):raise ValueError('protocol tuple mismatch')
        if baseline is not None and identity!=baseline:raise ValueError('incompatible scout shard provenance')
        baseline=identity;rows=read_points(cp)
        if any(r['shots']!=meta['shots'] for r in rows):raise ValueError('shot count mismatch')
        expected={(ch,*dims,float(p)) for ch,g in meta['config']['channel_grids'].items() for dims in g['dimensions'] for p in g['rates']}
        if {key(r) for r in rows}!=expected or len(rows)!=len(expected):raise ValueError('scout point grid incomplete/duplicated')
        inputs.append(dict(csv_sha256=sha(cp),metadata_sha256=sha(mp),path_ids=[str(cp),str(mp)],stage=meta['stage'],extension_round=meta['extension_round'],parent_selection_sha256=meta.get('parent_selection_sha256','')))
        metas.append(meta);shards.append(rows)
    previous=json.loads(Path(previous_path).read_text()) if previous_path else None
    round_=validate_ancestry(inputs,metas,previous,sha(previous_path) if previous_path else None)
    if previous:
        expected={(ch,*d,float(p)) for ch,g in previous['extension_grids'].items() for d in g['dimensions_T_Lu_Lv'] for p in g['rates']}
        if {key(r) for r in shards[-1]}!=expected:raise ValueError('extension does not equal declared grid')
    return merge_points(shards),metas[0],inputs,round_

def crossings(rows,protocol,unsaturated=False):
    sizes=sorted({r['T'] for r in rows});out=[]
    for small,large in zip(sizes,sizes[1:]):
        curves=[{r['p']:r[f'{protocol}_failures']/r['shots'] for r in rows if r['T']==L} for L in (small,large)]
        grid=sorted(set(curves[0])&set(curves[1]));diff=[curves[1][p]-curves[0][p] for p in grid]
        valid=[not unsaturated or all(.02<=curve[p]<=.98 for curve in curves) for p in grid]
        for i in range(len(grid)-1):
            if valid[i] and valid[i+1] and diff[i]*diff[i+1]<0:
                p=grid[i]-diff[i]*(grid[i+1]-grid[i])/(diff[i+1]-diff[i])
                if p>0:out.append({'sizes':[small,large],'p':p})
        for i in range(1,len(grid)-1):
            if valid[i] and valid[i-1] and valid[i+1] and diff[i]==0 and diff[i-1]*diff[i+1]<0 and grid[i]>0:out.append({'sizes':[small,large],'p':grid[i]})
    return out

def clusters(candidates):
    grouped=[]
    for item in sorted(candidates,key=lambda i:i['p']):
        if not grouped or item['p']/grouped[-1][0]['p']>1.5:grouped.append([])
        grouped[-1].append(item)
    def rank(g):return (-len({tuple(i['sizes']) for i in g}),-len(g),math.log(g[-1]['p']/g[0]['p']),math.exp(np.mean(np.log([i['p'] for i in g]))))
    return sorted(grouped,key=rank)

def timing_selection(rows,nrates):
    measured=[]
    for L in (3,5,7):
        rr=[r for r in rows if r['T']==L]
        if not rr:raise ValueError('timing selection requires scout sizes 3,5,7')
        measured.append([L,float(np.median([(r['sampling_seconds']+r['replay_seconds'])/r['shots'] for r in rr])),float(np.median([r['decoder_build_seconds'] for r in rr]))])
    x=np.log([r[0] for r in measured]);A=np.column_stack((np.ones(3),x))
    shot=np.linalg.lstsq(A,np.log(np.maximum([r[1] for r in measured],1e-12)),rcond=None)[0]
    build=np.linalg.lstsq(A,np.log(np.maximum([r[2] for r in measured],1e-12)),rcond=None)[0]
    prefixes=[list(range(3,end+1,2)) for end in (9,11,13)]
    project=lambda ls:nrates*sum(math.exp(build[0]+build[1]*math.log(L))+1000*math.exp(shot[0]+shot[1]*math.log(L)) for L in ls)
    viable=[ls for ls in prefixes if project(ls)<=43200];chosen=viable[-1] if viable else prefixes[0]
    return chosen,dict(shot=shot.tolist(),build=build.tolist(),measured_inputs=measured,projections=[dict(sizes=ls,seconds=project(ls)) for ls in prefixes],runtime_budget_exceeded=not bool(viable)),project(chosen)

def select_grids(points,meta,inputs,extension_round=0):
    if meta['stage']!='scout':raise ValueError('selection derives only from scout data')
    if any((r['T'],r['Lu'],r['Lv'])!=(r['T'],)*3 for r in points):raise ValueError('scout selection requires cubic sizes')
    # Extension bounds constrain scouts; exact five-rate production neighborhoods
    # may extend above .49 within the physical Bernoulli domain [0,1].
    selected=dict(format_version=1,rule_version='h5-scout-selection-v1',status='ready',extension_round=extension_round,scout_inputs=inputs,seed=meta['seed'],decoder_settings=meta['decoder_settings'],source_hashes=meta['source_hashes'],production_shots=1000,extension_grids={},channels={})
    for ch in CHANNELS:
        rows=[r for r in points if r['channel']==ch];sizes=sorted({r['T'] for r in rows});rates=sorted({r['p'] for r in rows});positive=[p for p in rates if p>0]
        if len(sizes)<3 or not positive:raise ValueError('insufficient scout grid')
        v=dict(fit_windows={},candidate_centers={},extension_states={},flags={},saturation_bounds=[1e-8,.49]);production=set();extend=set()
        for pr in PROTOCOLS:
            grouped=clusters(crossings(rows,pr,True));flags={'no_scout_crossing':not bool(grouped),'extension_exhausted':False,'threshold_ineligible':False,'clusters':grouped}
            centers=[rounded(math.exp(np.mean(np.log([i['p'] for i in g])))) for g in grouped];windows=[]
            if not centers:
                last=[{r['p']:r[f'{pr}_failures']/r['shots'] for r in rows if r['T']==L} for L in sizes[-2:]]
                common=sorted(set(last[0])&set(last[1]));unsat=[p for p in common if min(rates)<p<max(rates) and p>0 and all(.02<=a[p]<=.98 for a in last)]
                if unsat:centers=[min(unsat,key=lambda p:(abs(last[1][p]-last[0][p]),p))]
                else:
                    median=float(np.median([np.mean([a[p] for a in last]) for p in common]));new=[]
                    if median>=.5:new.extend((max(1e-8,min(positive)/10),max(1e-8,min(positive)/3)))
                    if median<=.5:new.extend((min(.49,1.5*max(rates)),min(.49,2*max(rates))))
                    new=sorted({rounded(p) for p in new}-set(rates));flags['saturation_median']=median
                    if new and extension_round<3:
                        v['extension_states'][pr]='requesting';extend.update(new)
                    else:
                        v['extension_states'][pr]='exhausted_diagnostic';flags.update(extension_exhausted=True,threshold_ineligible=True)
                        center=min(positive) if median>.5 else max(positive)
                        selected_rates=sorted(sorted(positive,key=lambda p:abs(math.log(p/center)))[:5]);production.update(selected_rates)
                        centers=[center];windows=[[min(selected_rates),max(selected_rates)]]
            if centers and not windows:
                v['extension_states'][pr]='resolved_candidates'
                for center in centers:
                    around=sorted({rounded(center*f) for f in (.8,.9,1,1.1,1.2)});production.update(around);windows.append([around[0],around[-1]])
            v['fit_windows'][pr]=windows;v['candidate_centers'][pr]=centers;v['flags'][pr]=flags
        ls,model,seconds=timing_selection(rows,len(production))
        v.update(dimensions_T_Lu_Lv=[[L]*3 for L in ls],rates=sorted(production),timing_model=model,projected_seconds=seconds)
        selected['channels'][ch]=v
        if extend:selected['extension_grids'][ch]=dict(dimensions_T_Lu_Lv=[[L]*3 for L in (3,5,7)],rates=sorted(extend))
    if selected['extension_grids']:selected['status']='needs_extended_scout'
    return selected

def fit_window(points,channel,protocol,window,bootstrap=1000,seed=1234,degree=2,min_size=0,correction=False,ineligible=False,z_guard_issues=()):
    lo,hi=window;rows=[r for r in points if r['channel']==channel and lo<=r['p']<=hi and r['T']>=min_size]
    sizes=sorted({r['T'] for r in rows});rates=sorted({r['p'] for r in rows});issues=[]
    fit=dict(channel=channel,protocol=protocol,fit_window=window,selected_sizes=sizes,selected_probabilities=rates,degree=degree,correction=correction,status='diagnostic_candidate',pc=None,nu=None,candidate=None,diagnostic=issues,bootstrap=dict(replicates=bootstrap,accepted=0,samples=[],pc_interval=None,nu_interval=None))
    if channel=='z_only':issues.extend([z_guard_issues] if isinstance(z_guard_issues,str) else z_guard_issues)
    if len(sizes)<4:issues.append('need at least four sizes')
    if len(rates)<5:issues.append('need at least five probabilities')
    if ineligible:issues.append('scout extension exhausted: threshold ineligible')
    if not 0<lo<hi<=1:issues.append('invalid fit window')
    if issues:return fit
    cross=crossings(rows,protocol);fit['crossings']=cross;expected=set(zip(sizes,sizes[1:]))
    if not any({tuple(x['sizes']) for x in g}==expected for g in clusters(cross)):issues.append('no single crossing cluster covers every adjacent size pair')
    try:
        c=engine._candidate_fit(rows,channel,protocol,p_min=lo,p_max=hi,degree=degree,nu_bounds=(.3,5.),goodness_alpha=.01,correction=correction)
    except (ArithmeticError,RuntimeError,ValueError) as e:issues.append(str(e));return fit
    fit['candidate']=c;issues.extend(engine._candidate_guard_issues(c,p_min=lo,p_max=hi,goodness_alpha=.01))
    if issues:return fit
    # Stable channel/protocol/window stream: iteration order never changes a CI.
    salt=int(hashlib.sha256(f'{channel}/{protocol}/{window}'.encode()).hexdigest()[:16],16)
    rng=np.random.default_rng(np.random.SeedSequence([seed,salt]));samples=[]
    for i in range(bootstrap):
        sample=[dict(r,**{f'{protocol}_failures':int(rng.binomial(r['shots'],p))}) for r,p in zip(rows,c['fitted_probabilities'])]
        try:
            b=engine._candidate_fit(sample,channel,protocol,p_min=lo,p_max=hi,degree=degree,nu_bounds=(.3,5.),goodness_alpha=.01,correction=correction,warm_start=c['parameters'])
            problems=engine._candidate_guard_issues(b,p_min=lo,p_max=hi,goodness_alpha=.01,include_goodness=False)
            samples.append(dict(replicate=i,accepted=not problems,pc=b['pc'],nu=b['nu']))
        except (ArithmeticError,RuntimeError,ValueError):samples.append(dict(replicate=i,accepted=False,pc=None,nu=None))
    accepted=[b for b in samples if b['accepted']];fit['bootstrap'].update(accepted=len(accepted),samples=samples,method='parametric marginal binomial; intervals conditional on frozen window and scaling model')
    if bootstrap<2 or len(accepted)<max(2,math.ceil(.8*bootstrap)):
        issues.append('fewer than 80% usable bootstrap fits or fewer than two requested');return fit
    fit['bootstrap'].update(pc_interval=np.quantile([b['pc'] for b in accepted],[.025,.975]).tolist(),nu_interval=np.quantile([b['nu'] for b in accepted],[.025,.975]).tolist())
    fit.update(status='ok',pc=c['pc'],nu=c['nu']);return fit

def analyze(points,selection,bootstrap=1000,seed=1234,z_boundary_audit=None):
    fits=[];sensitivity=[]
    for ch in CHANNELS:
        config=selection['channels'][ch]
        for pr in PROTOCOLS:
            guard_issues=()
            if ch=='z_only' and z_boundary_audit is not None:
                guard_issues=z_boundary_guard(z_boundary_audit,
                    [dimensions[0] for dimensions in config['dimensions_T_Lu_Lv']],
                    pr,selection.get('source_hashes'))
            candidates=[]
            for index,window in enumerate(config['fit_windows'][pr]):
                f=fit_window(points,ch,pr,window,bootstrap,seed,
                    ineligible=config['flags'][pr].get('threshold_ineligible',False),
                    z_guard_issues=guard_issues);f['candidate_index']=index;candidates.append(f)
                sizes=f['selected_sizes'];rates=f['selected_probabilities']
                variants=[('cubic',window,0,3,False),('drop_smallest_size',window,sizes[1] if len(sizes)>1 else 9999,2,False),('fixed_L_inverse_correction',window,0,2,True)]
                if len(rates)>1:variants.extend([('trim_lowest_p',[rates[1],window[1]],0,2,False),('trim_highest_p',[window[0],rates[-2]],0,2,False)])
                for name,w,minimum,degree,corr in variants:
                    s=fit_window(points,ch,pr,w,0,seed,degree,minimum,corr,
                        z_guard_issues=guard_issues);sensitivity.append(dict(channel=ch,protocol=pr,candidate_index=index,name=name,status=s['status'],candidate_pc=(s['candidate'] or {}).get('pc'),diagnostic=s['diagnostic']))
            accepted=[f for f in candidates if f['status']=='ok']
            claim=accepted[0] if len(accepted)==1 else None
            fits.append(dict(channel=ch,protocol=pr,status='threshold' if claim else 'diagnostic_candidates',pc=claim['pc'] if claim else None,nu=claim['nu'] if claim else None,pc_interval=claim['bootstrap']['pc_interval'] if claim else None,candidates=candidates))
    return dict(fits=fits,sensitivity=sensitivity,selection_rule=selection['rule_version'])

def write_csv(path,rows,fields):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def write_report(report,points,output,hashes):
    output=Path(output);output.mkdir(parents=True,exist_ok=True);report=dict(report,hashes=hashes)
    atomic_json(output/'fits.json',report)
    write_csv(output/'sensitivity.csv',report['sensitivity'],['channel','protocol','candidate_index','name','status','candidate_pc','diagnostic'])
    samples=[dict(channel=f['channel'],protocol=f['protocol'],candidate_index=c['candidate_index'],**b) for f in report['fits'] for c in f['candidates'] for b in c['bootstrap']['samples']]
    write_csv(output/'bootstrap.csv',samples,['channel','protocol','candidate_index','replicate','accepted','pc','nu'])
    pairs=[dict(channel=r['channel'],T=r['T'],Lu=r['Lu'],Lv=r['Lv'],p=r['p'],contrast=k,difference=v['difference'],ci_low=v['ci'][0],ci_high=v['ci'][1]) for r in points for k,v in r['pairs'].items()]
    write_csv(output/'paired-differences.csv',pairs,['channel','T','Lu','Lv','p','contrast','difference','ci_low','ci_high'])
    lines=['# Direct-projector H+5 guarded analysis','', 'Primary observables: either-pair logical X for X-only faults; state-independent either-pair logical Z for Z-only faults. Z is not |00_L> infidelity.','', '| Channel | Protocol | Status | Threshold | 95% conditional interval |','|---|---|---|---|---|']
    for f in report['fits']:lines.append(f"| {f['channel']} | {f['protocol']} | {f['status']} | {f['pc']} | {f['pc_interval']} |")
    lines+=['','Statistical bootstrap intervals are conditional on the frozen window and scaling model. Model, size, and window sensitivity is listed separately in sensitivity.csv. A diagnostic candidate is not a threshold.','']
    for f in report['fits']:
        for c in f['candidates']:lines.append(f"- {f['channel']}/{f['protocol']} candidate {c['candidate_index']}: window {c['fit_window']}; candidate pc={(c['candidate'] or {}).get('pc')}; guards: {c['diagnostic']}; bootstrap {c['bootstrap']['accepted']}/{c['bootstrap']['replicates']}.")
    (output/'report.md').write_text('\n'.join(lines)+'\n')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for ch in CHANNELS:
        fig,axes=plt.subplots(2,4,figsize=(16,8),constrained_layout=True)
        for j,pr in enumerate(PROTOCOLS):
            ff=next(f for f in report['fits'] if f['channel']==ch and f['protocol']==pr)
            for L in sorted({r['T'] for r in points if r['channel']==ch}):
                rr=sorted([r for r in points if r['channel']==ch and r['T']==L],key=lambda r:r['p'])
                y=[r[f'{pr}_failures']/r['shots'] for r in rr];ci=[engine.wilson_interval(r[f'{pr}_failures'],r['shots']) for r in rr]
                axes[0,j].errorbar([r['p'] for r in rr],y,yerr=[[a-b[0] for a,b in zip(y,ci)],[b[1]-a for a,b in zip(y,ci)]],marker='o',markersize=3,label=f'T=Lu=Lv={L}')
                for c in ff['candidates']:
                    if c['candidate']:
                        cc=c['candidate'];grid=np.linspace(*c['fit_window'],100);axes[1,j].plot(grid,[engine._model_probability(cc,p,L) for p in grid],label=f'L={L}, candidate {c["candidate_index"]}')
            axes[0,j].set_title(pr+'\n'+ff['status']);axes[1,j].set_title('Scaling candidate (diagnostic unless guarded)');axes[0,j].set_ylim(-.02,1.02)
            for ax in axes[:,j]:ax.set_xlabel('post-gate p');ax.set_ylabel('logical failure probability')
        axes[0,0].legend(fontsize=8);fig.suptitle('H+5 '+ch+(' — |00_L> failure' if ch=='x_only' else ' — state-independent logical Z'))
        for ext in ('pdf','svg','png'):fig.savefig(output/f'fits-{ch}.{ext}',dpi=180)
        plt.close(fig)

def main(args=None):
    parser=argparse.ArgumentParser(description=__doc__);subs=parser.add_subparsers(dest='command',required=True)
    select=subs.add_parser('select');select.add_argument('--input-csv',action='append',required=True);select.add_argument('--input-metadata',action='append',required=True);select.add_argument('--previous-selection');select.add_argument('--output',required=True)
    fit=subs.add_parser('fit');fit.add_argument('--selection',required=True);fit.add_argument('--csv',required=True);fit.add_argument('--metadata',required=True);fit.add_argument('--bootstrap',type=int,default=1000);fit.add_argument('--seed',type=int,default=1234);fit.add_argument('--z-boundary-audit');fit.add_argument('--output',required=True)
    report=subs.add_parser('report');report.add_argument('--fits',required=True);report.add_argument('--csv',required=True);report.add_argument('--output',required=True)
    a=parser.parse_args(args)
    if a.command=='select':
        rows,meta,inputs,round_=load_scout(a.input_csv,a.input_metadata,a.previous_selection);selection=select_grids(rows,meta,inputs,round_);atomic_json(a.output,selection);print(selection['status']);return
    if a.command=='report':
        data=json.loads(Path(a.fits).read_text())
        if data.get('hashes',{}).get('raw_csv')!=sha(a.csv):raise ValueError('report CSV hash mismatch')
        write_report(data,read_points(a.csv),a.output,data['hashes']);return
    selection=json.loads(Path(a.selection).read_text());meta=tomllib.loads(Path(a.metadata).read_text())
    if selection['status']!='ready' or meta['stage']!='production' or not meta['complete'] or meta['invalid_recoveries']:raise ValueError('fit requires complete production from a ready selection')
    if meta.get('selection_sha256')!=sha(a.selection) or meta['raw_csv_sha256']!=sha(a.csv):raise ValueError('production ancestry/data hash mismatch')
    for k in ('source_hashes','decoder_settings','seed'):
        if meta[k]!=selection[k]:raise ValueError('production provenance differs from selection')
    if meta['shots']!=selection['production_shots'] or meta['protocols']!=list(PROTOCOLS):raise ValueError('production shots/protocols mismatch')
    rows=read_points(a.csv);expected={(ch,*d,p) for ch,v in selection['channels'].items() for d in v['dimensions_T_Lu_Lv'] for p in v['rates']}
    if {key(r) for r in rows}!=expected or len(rows)!=len(expected):raise ValueError('production grid differs from frozen selection')
    z_audit=load_z_boundary_audit(a.z_boundary_audit) if a.z_boundary_audit else None
    report=analyze(rows,selection,a.bootstrap,a.seed,z_audit)
    hashes=dict(selection=sha(a.selection),raw_csv=sha(a.csv),metadata=sha(a.metadata),analysis=sha(__file__),likelihood_engine=sha(Path(__file__).with_name('fit_toric3d_diagonal_feedback.py')))
    if a.z_boundary_audit:hashes['z_boundary_audit']=sha(a.z_boundary_audit)
    write_report(report,rows,a.output,hashes)
if __name__=='__main__':main()
