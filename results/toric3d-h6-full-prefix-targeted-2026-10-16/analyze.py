#!/usr/bin/env python3
"""Frozen targeted H+6 analysis. Fresh data primary; same-cell pooling secondary.

Report schema v1: fresh/pooled_secondary -> x/z -> points, fits, profiles,
failures. Fresh X additionally has primary_crossing, secondary_crossing,
common_model, independent_size_model and bootstrap. Every result carries its
channel, dataset and fit identifier; provenance binds input, plan and sources.
"""
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
import tomllib

import numpy as np
import scipy

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PREVIOUS=HERE.parent/'toric3d-h6-full-prefix-reanalysis-2026-10-16'
OLD=HERE.parent/'toric3d-h6-full-prefix-2026-10-15'
OLD_SHARDS=('production-x-L11-a','production-x-L11-b','production-x-L11-c',
            'production-x-L9','production-x-small','production-z')
INTEGER_FIELDS=('steps','Ly','Lz','shots','target_shots','seed','final_only_logical_1_failures',
 'final_only_logical_2_failures','final_only_failures','full_prefix_logical_1_failures',
 'full_prefix_logical_2_failures','full_prefix_failures')


def import_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


X=import_module(PREVIOUS/'x_analysis.py','targeted_unchanged_x_helpers')
Z=import_module(PREVIOUS/'z_models.py','targeted_unchanged_z_helpers')


def require(condition,message):
    if not condition: raise ValueError(message)


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cell_key(row): return row['channel'],int(row['steps']),round(float(row['p']),12)


def ordered(rows): return sorted(rows,key=cell_key)


def read_plan(study):
    study=Path(study);path=study/'selection.json';digest_path=study/'selection.sha256'
    require(digest_path.is_file(),'missing frozen selection digest')
    expected=digest_path.read_text().split()[0]
    require(sha(path)==expected,'frozen selection hash mismatch')
    return json.loads(path.read_text())


def channel_plan(plan,channel):
    require(channel in ('x_only','z_only'),'unknown channel')
    return plan['analysis']['x' if channel=='x_only' else 'z']


def validate_counts(row):
    n=row['shots']; k=row['full_prefix_failures']
    require(type(n) is int and n>0 and type(k) is int and 0<=k<=n,'invalid binomial counts')
    require(row['complete'] is True,'incomplete production point')
    require(row['steps']==row['Ly']==row['Lz'],'geometry is not T=L=Ly=Lz')
    require(math.isfinite(row['p']) and 0<row['p']<1,'invalid physical probability')
    require(len(row['point_seed'])==16 and all(c in '0123456789abcdef' for c in row['point_seed']),'invalid stream seed')
    j=row['joint']
    require(len(j)==4 and all(type(v) is int and v>=0 for v in j) and sum(j)==n,'joint histogram count mismatch')
    require(j[1]+j[3]==row['final_only_failures'] and j[2]+j[3]==k,'joint marginal count mismatch')
    for arm in ('final_only','full_prefix'):
        a,b,u=[row[arm+suffix] for suffix in ('_logical_1_failures','_logical_2_failures','_failures')]
        require(all(type(v) is int for v in (a,b,u)) and 0<=max(a,b)<=u<=min(n,a+b),'logical marginal count mismatch')
    require(row.get('failures',k)==k and row.get('none_failures',k)==k,'analysis count alias mismatch')


def assert_independent(fresh,old,excluded_seeds=()):
    fresh_seeds=[r['point_seed'] for r in fresh]; old_seeds=[r['point_seed'] for r in old]
    require(len(set(fresh_seeds))==len(fresh_seeds),'duplicate fresh random stream')
    require(len(set(old_seeds))==len(old_seeds),'duplicate old production random stream')
    require(not set(fresh_seeds)&(set(old_seeds)|set(excluded_seeds)),'fresh random stream collides with old/excluded stream')
    return {'fresh_streams':len(fresh_seeds),'old_production_streams':len(old_seeds),
            'excluded_unique_streams':len(set(excluded_seeds)),'point_seeds_disjoint':True,
            'interpretation':'Distinct recorded Xoshiro point seeds verified; this is a provenance/stream-reuse check, not a statistical proof of RNG independence.'}


def validate_fresh_points(points,plan):
    expected={}
    seeds=[]
    for job in plan['jobs']:
        seeds.append(job['seed'])
        for size in job['sizes']:
            for p in job['rates']:
                key=(job['channel'],size,round(p,12))
                require(key not in expected,'duplicate planned grid cell')
                expected[key]=job
    require(len(set(seeds))==len(seeds),'duplicate planned base seed')
    keys=[cell_key(r) for r in points]
    require(len(set(keys))==len(keys),'duplicate fresh grid cell')
    require(set(keys)==set(expected),'incomplete or unexpected fresh grid')
    for row in points:
        validate_counts(row);job=expected[cell_key(row)]
        require(row['shots']==row['target_shots']==job['shots'],'production shot count differs from plan')
        require(row['seed']==job['seed'],'production base seed differs from plan')
        require(row['source_job']==job['name'],'production shard differs from plan')
    assert_independent(points,[])
    by_channel={channel:sum(r['shots'] for r in points if r['channel']==channel) for channel in ('x_only','z_only')}
    require(by_channel==plan['new_fault_records'],'channel record totals differ from plan')
    require(sum(by_channel.values())==plan['total_new_fault_records'],'total record count differs from plan')
    for channel in by_channel:
        cp=channel_plan(plan,channel)
        want={(channel,L,round(p,12)) for L in cp['sizes'] for p in cp['rates']}
        require({k for k in keys if k[0]==channel}==want,'analysis grid differs from production')
    return {'records':sum(by_channel.values()),'by_channel':by_channel,'cells':len(points),'complete':True}


def read_shard(directory,job,source_hashes=None,protocol=None,earliest_start=None):
    directory=Path(directory)
    files={name:directory/name for name in ('raw.csv','metadata.toml','checkpoint.toml')}
    for path in files.values(): require(path.is_file(),f'missing production file: {path}')
    metadata=tomllib.loads(files['metadata.toml'].read_text())
    checkpoint=tomllib.loads(files['checkpoint.toml'].read_text())
    require(metadata['complete'] is True,f'incomplete shard: {directory.name}')
    require(sha(files['raw.csv'])==metadata['raw_sha256'],'raw checksum mismatch')
    require(sha(files['checkpoint.toml'])==metadata['checkpoint_sha256'],'checkpoint checksum mismatch')
    snapshot=directory/metadata['current_snapshot']
    require(snapshot.resolve().is_relative_to(directory.resolve()),'snapshot outside shard')
    require(snapshot.is_file() and sha(snapshot)==sha(files['checkpoint.toml']),'immutable snapshot mismatch')
    for key in ('config','runtime','source_hashes','protocol','point_seeds_hex','generation','started_at_unix','updated_at_unix'):
        require(metadata[key]==checkpoint[key],f'checkpoint metadata mismatch: {key}')
    for key in ('channel','sizes','rates','shots','seed','batch_size'):
        require(metadata['config'][key]==job[key],f'shard configuration mismatch: {key}')
    require(metadata['config']['dimensions']==[[L,L,L] for L in job['sizes']],'shard dimensions mismatch')
    if source_hashes is not None: require(metadata['source_hashes']==source_hashes,'simulation source manifest mismatch')
    if protocol is not None: require(metadata['protocol']==protocol,'physical protocol/decoder mismatch')
    if earliest_start is not None: require(metadata['started_at_unix']>=earliest_start,'production predates frozen plan')
    with files['raw.csv'].open(newline='') as stream: raw=list(csv.DictReader(stream))
    expected=[(L,p) for L in job['sizes'] for p in job['rates']]
    require(len(raw)==len(expected)==len(checkpoint['states'])==len(metadata['point_seeds_hex']),'shard grid size mismatch')
    points=[]
    for index,(r,(L,p),state) in enumerate(zip(raw,expected,checkpoint['states'])):
        row=dict(r)
        for field in INTEGER_FIELDS: row[field]=int(row[field])
        row['p']=float(row['p']);row['elapsed_seconds']=float(row['elapsed_seconds'])
        row['complete']=row['complete']=='true';row['joint']=[int(row.pop(f'joint_{i:02d}')) for i in range(4)]
        row.update(source_job=job['name'],failures=row['full_prefix_failures'],none_failures=row['full_prefix_failures'])
        require(row['channel']==job['channel'] and row['steps']==L and row['p']==p,'raw point order/grid mismatch')
        require(row['seed']==job['seed'] and row['point_seed']==metadata['point_seeds_hex'][index],'raw seed mismatch')
        require(row['shots']==row['target_shots']==job['shots']==state['completed_shots'],'incomplete/wrong shard shot count')
        for rawkey,statekey in (('joint','joint_counts'),): require(row[rawkey]==state[statekey],'checkpoint joint mismatch')
        for statekey,suffix in (('failures','_failures'),('logical_1_counts','_logical_1_failures'),('logical_2_counts','_logical_2_failures')):
            require([row[arm+suffix] for arm in ('final_only','full_prefix')]==state[statekey],'checkpoint counts mismatch')
        validate_counts(row);points.append(row)
    require(metadata['expected_points']==metadata['completed_points']==len(points),'metadata completed point mismatch')
    require(metadata['completed_shots']==sum(r['shots'] for r in points),'metadata completed shot mismatch')
    return points,{'directory':str(directory),'files_sha256':{name:sha(path) for name,path in files.items()},
                   'snapshot':metadata['current_snapshot'],'snapshot_sha256':sha(snapshot),
                   'point_seeds':metadata['point_seeds_hex'],'base_seed':job['seed'],
                   'protocol':metadata['protocol'],'runtime':metadata['runtime'],
                   'source_hashes':metadata['source_hashes'],'started_at_unix':metadata['started_at_unix']}


def load_old_production(plan,old_dir=OLD):
    old_dir=Path(old_dir);manifest_path=old_dir/'ARTIFACTS.sha256.json'
    require(sha(manifest_path)==plan['old_manifest_sha256'],'old artifact manifest changed')
    manifest=json.loads(manifest_path.read_text());points=[];provenance=[]
    for name in OLD_SHARDS:
        directory=old_dir/name
        for filename in ('raw.csv','metadata.toml','checkpoint.toml'):
            require(sha(directory/filename)==manifest[f'{name}/{filename}'],'old production artifact hash mismatch')
        metadata=tomllib.loads((directory/'metadata.toml').read_text())
        job=dict(metadata['config'],name=name)
        rows,record=read_shard(directory,job,plan['simulation_source_sha256'])
        points.extend(rows);provenance.append(record)
    require(len(points)==70 and sum(r['shots'] for r in points)==70000,'old production grid/count mismatch')
    assert_independent([],points)
    return ordered(points),{'manifest_sha256':sha(manifest_path),'shards':provenance}


def collect_excluded_streams(study_dir,old_dir,plan):
    dirs=[Path(study_dir)/name for name in plan['excluded_timing_runs']]
    dirs.extend(sorted(path for path in Path(old_dir).glob('pilot*') if path.is_dir()))
    seeds=[];records=[]
    for directory in dirs:
        path=directory/'metadata.toml'
        if not path.is_file(): path=directory/'checkpoint.toml'
        require(path.is_file(),f'missing excluded stream metadata: {directory}')
        data=tomllib.loads(path.read_text());seeds.extend(data['point_seeds_hex'])
        records.append({'path':str(path),'sha256':sha(path),'point_seeds':data['point_seeds_hex']})
    return seeds,records


def pool_points(fresh,old):
    assert_independent(fresh,old)
    output={cell_key(r):copy.deepcopy(r) for r in fresh};added=[]
    for row in output.values():
        row['fresh_shots']=row['shots'];row['old_shots']=0
        row['stream_components']=[{'point_seed':row['point_seed'],'source_job':row['source_job'],'shots':row['shots']}]
    for row in old:
        key=cell_key(row)
        if key not in output: continue
        target=output[key];target['old_shots']+=row['shots'];target['shots']+=row['shots']
        for field in ('failures','none_failures','full_prefix_failures','final_only_failures',
                      'full_prefix_logical_1_failures','full_prefix_logical_2_failures',
                      'final_only_logical_1_failures','final_only_logical_2_failures'):
            target[field]+=row[field]
        target['joint']=[a+b for a,b in zip(target['joint'],row['joint'])]
        target['target_shots']=target['shots']
        target['stream_components'].append({'point_seed':row['point_seed'],'source_job':row['source_job'],'shots':row['shots']})
        added.append(row)
    for row in output.values():
        row.pop('point_seed',None);row.pop('seed',None);row['source_job']='independent_same_cell_pool'
    return ordered(output.values()),{'old_cells_added':len(added),'old_records_added':sum(r['shots'] for r in added),
        'pooled_records':sum(r['shots'] for r in output.values()),
        'policy':'Old production only, restricted to exactly the fresh channel/size/probability cells; no new rates/sizes from old data.'}


def select_data(points,plan,channel,min_size,window):
    low,high=channel_plan(plan,channel)['windows'][window]
    rows=ordered([r for r in points if r['channel']==channel and r['steps']>=min_size and low<=r['p']<=high])
    require(rows,'empty fit subset')
    sizes=np.array([r['steps'] for r in rows],float)
    reference=float(np.exp(np.mean(np.log(np.unique(sizes)))))
    return {'sizes':sizes,'p':np.array([r['p'] for r in rows]),
            'shots':np.array([r['shots'] for r in rows],float),
            'failures':np.array([r['failures'] for r in rows],float),
            'L_reference':reference,'relative_log':np.log(sizes/reference),
            'p_min':low,'p_max':high,'p_scale':high-low,'min_size':min_size,'window':window}


def fit_schedule(plan,channel):
    cp=channel_plan(plan,channel)
    return [(minimum,window,model) for minimum in cp['min_sizes'] for window in cp['windows'] for model in cp['models']]


def profile_schedule(plan):
    cp=channel_plan(plan,'z_only')
    return [(minimum,window,model) for minimum in cp['min_sizes'] for window in cp['windows'] for model in cp['profile_models']]


def anchor_assessment(points,plan):
    spec=plan['analysis']['x']['primary_crossing'];small,large=spec['size_pair'];rates=spec['rates'];alpha=spec['alpha']
    lookup={cell_key(r):r for r in points};cells=[];by_key={}
    for p in rates:
        for L in (small,large):
            row=lookup['x_only',L,round(p,12)];k,n=row['failures'],row['shots']
            lo,hi=X.clopper_pearson(k,n,alpha/4)
            cell={'L':L,'p':p,'shots':n,'failures':k,'rate':k/n,'ci_low':lo,'ci_high':hi}
            cells.append(cell);by_key[L,p]=cell
    contrasts=[]
    for p in rates:
        a,b=by_key[small,p],by_key[large,p]
        contrasts.append({'L_small':small,'L_large':large,'p':p,'difference':b['rate']-a['rate'],
                          'ci_low':b['ci_low']-a['ci_high'],'ci_high':b['ci_high']-a['ci_low']})
    return {'method':spec['method'],'family_alpha':alpha,'cell_alpha':alpha/4,'number_cells':4,
            'cells':cells,'contrasts':contrasts,
            'success':contrasts[0]['ci_high']<0 and contrasts[1]['ci_low']>0,
            'interpretation':'Prespecified simultaneous opposite-sign anchor contrast; neither proves a unique common crossing nor an asymptotic threshold.'}


def secondary_crossing(points):
    result=X.crossing_assessment(points)
    result['method']=f'equal-tail exact Clopper-Pearson; Bonferroni over all {len(points)} binomial cells'
    result['secondary_note']='All size-pair/adjacent compatibility regions were prespecified as secondary diagnostics; the primary check uses four anchor cells.'
    return result


def fit_channel(points,plan,channel,dataset):
    report={'channel':channel,'dataset':dataset,'points':[r for r in points if r['channel']==channel],
            'fits':[],'profiles':[],'failures':[],'fit_schedule':[list(v) for v in fit_schedule(plan,channel)]}
    for minimum,window,model in fit_schedule(plan,channel):
        identity=f'Lmin{minimum}-{window}-{model}'
        try:
            fit=Z.fit_model(select_data(points,plan,channel,minimum,window),model)
            fit.update(fit_id=identity,channel=channel,dataset=dataset)
            report['fits'].append(fit)
            print(dataset,channel,identity,'pc',fit['pc'],'GOF',fit['deviance_p_value'],flush=True)
        except (ArithmeticError,RuntimeError,ValueError) as error:
            report['failures'].append({'stage':'fit','fit_id':identity,'error':str(error)})
    for minimum in channel_plan(plan,channel)['min_sizes']:
        for window in channel_plan(plan,channel)['windows']:
            subset=[f for f in report['fits'] if f['min_size']==minimum and f['window']==window]
            finite=[f['aicc'] for f in subset if f['aicc'] is not None]
            for f in subset: f['delta_aicc_same_subset']=f['aicc']-min(finite) if finite and f['aicc'] is not None else None
    report['fit_count']=len(report['fits'])
    report['expected_fit_count']=len(fit_schedule(plan,channel))
    return report


def profile_for_fit(data,model,fit,grid_count=301):
    try:
        result=Z.profile_pc(data,model,fit,grid_count)
    except TypeError as error:
        # The unchanged helper can reach its final minimum-agreement calculation
        # with nll=None after a failed free-pc nuisance solve. Preserve this as a
        # critical helper exception; do not manufacture an interval or success.
        result={'fit_id':fit['fit_id'],'model':model,'min_size':data['min_size'],
            'confidence_level':.95,'pc_bounds':[data['p_min'],data['p_max']],
            'grid_count_requested':grid_count,'free_pc':fit['pc'],'free_nll':fit['nll'],
            'minimum_agreement':False,'regions':[],'evaluations':[],
            'critical_helper_exception':True,
            'numerical_failures':[{'stage':'profile_helper','exception_type':type(error).__name__,
                                   'message':str(error),'critical':True}],
            'interpretation':'Profile unavailable after a critical helper exception. No interval is reported; subsequent scheduled profiles are still attempted.',
            'evaluation_retention_note':'The unchanged helper raised before returning its local evaluation cache; those intermediate evaluations are unavailable.'}
    # Historical helper hardcodes the label full, while data control the mathematics.
    result.update(window=data['window'],channel='z_only',dataset='fresh')
    return result


def run_bootstrap(points,plan):
    # A private module keeps bootstrap settings isolated from other analyses.
    engine=import_module(PREVIOUS/'x_analysis.py','targeted_private_bootstrap')
    engine.B=plan['analysis']['x']['bootstrap_replicates']
    engine.SEED=plan['analysis']['x']['bootstrap_seed']
    common=engine.fit_common(points);require(engine.common_usable(common),'primary common fit failed')
    independent=engine.fit_independent(points,common)
    require(independent['optimizer_success'],'primary independent-size fit failed')
    require(independent['deviance']<=common['deviance']+1e-6,'primary nested likelihood ordering failed')
    return {'common_model':common,'independent_size_model':independent,
            'observed_likelihood_ratio':max(0.,common['deviance']-independent['deviance']),
            'bootstrap':engine.bootstrap(points,common,independent)}



def summarize_numerics(report,plan):
    """Separate missing/failed calculations from statistical/model diagnostics."""
    def finite(row,keys):
        return all(isinstance(row.get(k),(int,float)) and math.isfinite(row[k]) for k in keys)
    fit_summary={'expected_count':0,'returned_count':0,'missing_count':0,'unexpected_count':0,
                 'duplicate_count':0,'optimizer_failure_count':0,'nonfinite_fit_count':0,'problems':[]}
    profiles={'expected_count':len(profile_schedule(plan)),'returned_count':0,'missing_count':0,
              'unexpected_count':0,'duplicate_count':0,'invalid_profile_count':0,
              'numerical_failure_count':0,'minimum_disagreement_count':0,
              'failed_evaluation_count':0,'numerical_gap_region_count':0,'critical_helper_exception_count':0,'problems':[]}
    flags={'fit_bound_hit_count':0,'fit_nonmonotone_count':0,'fit_gof_p_lt_0p01_count':0,
           'fit_gof_p_lt_0p05_count':0,'fit_gof_unavailable_count':0,
           'profile_nuisance_bound_evaluation_count':0,'profile_bound_truncated_region_count':0,
           'bootstrap_bound_hit_count':0,'bootstrap_nonmonotone_count':0}
    exceptions=0
    def count_ids(items,wanted,summary,prefix):
        ids=[item.get('fit_id') for item in items];unique=set(ids)
        summary['returned_count']+=len(items)
        summary['missing_count']+=len(wanted-unique)
        summary['unexpected_count']+=len(unique-wanted)
        summary['duplicate_count']+=len(ids)-len(unique)
        if wanted-unique or unique-wanted or len(ids)!=len(unique):
            summary['problems'].append(dict(prefix,missing=sorted(wanted-unique),
                unexpected=sorted(str(v) for v in unique-wanted),duplicate_count=len(ids)-len(unique)))
    for dataset in ('fresh','pooled_secondary'):
        for channel,label in (('x_only','x'),('z_only','z')):
            current=report.get(dataset,{}).get(label,{})
            exceptions+=len(current.get('failures',[]))
            wanted={f'Lmin{m}-{w}-{model}' for m,w,model in fit_schedule(plan,channel)}
            fit_summary['expected_count']+=len(wanted)
            fits=current.get('fits',[])
            count_ids(fits,wanted,fit_summary,{'dataset':dataset,'channel':channel})
            for fit in fits:
                optimizer_failure=fit.get('optimizer_success') is not True
                nonfinite=not finite(fit,('pc','nu','deviance','nll'))
                fit_summary['optimizer_failure_count']+=int(optimizer_failure)
                fit_summary['nonfinite_fit_count']+=int(nonfinite)
                if optimizer_failure or nonfinite:
                    fit_summary['problems'].append({'dataset':dataset,'channel':channel,'fit_id':fit.get('fit_id'),
                        'optimizer_failure':optimizer_failure,'nonfinite':nonfinite})
                flags['fit_bound_hit_count']+=bool(fit.get('bound_hits'))
                flags['fit_nonmonotone_count']+=fit.get('monotone_increasing') is False
                p=fit.get('deviance_p_value')
                if isinstance(p,(int,float)) and math.isfinite(p):
                    flags['fit_gof_p_lt_0p01_count']+=p<.01
                    flags['fit_gof_p_lt_0p05_count']+=p<.05
                else:flags['fit_gof_unavailable_count']+=1
    saved_profiles=report.get('fresh',{}).get('z',{}).get('profiles',[])
    wanted={f'Lmin{m}-{w}-{model}' for m,w,model in profile_schedule(plan)}
    count_ids(saved_profiles,wanted,profiles,{'dataset':'fresh','channel':'z_only'})
    for profile in saved_profiles:
        failures=len(profile.get('numerical_failures',[]))
        disagreement=profile.get('minimum_agreement') is not True
        evaluations=profile.get('evaluations',[])
        bad_evaluations=sum(e.get('optimizer_success') is not True or not finite(e,('nll','lr')) for e in evaluations)
        gaps=sum(bool(r.get('lower_numerical_gap') or r.get('upper_numerical_gap') or r.get('refinement_failures')) for r in profile.get('regions',[]))
        critical=bool(profile.get('critical_helper_exception'))
        invalid=bool(failures or disagreement or bad_evaluations or gaps or critical or not evaluations)
        profiles['invalid_profile_count']+=invalid
        profiles['numerical_failure_count']+=failures
        profiles['minimum_disagreement_count']+=disagreement
        profiles['failed_evaluation_count']+=bad_evaluations
        profiles['numerical_gap_region_count']+=gaps
        profiles['critical_helper_exception_count']+=critical
        if invalid: profiles['problems'].append({'fit_id':profile.get('fit_id'),'numerical_failure_count':failures,
            'minimum_disagreement':disagreement,'failed_evaluation_count':bad_evaluations,
            'numerical_gap_region_count':gaps,'critical_helper_exception':critical,'empty_evaluations':not evaluations})
        flags['profile_nuisance_bound_evaluation_count']+=int(profile.get('nuisance_bound_evaluation_count',0))
        flags['profile_bound_truncated_region_count']+=sum(bool(r.get('lower_bound_truncated') or r.get('upper_bound_truncated')) for r in profile.get('regions',[]))
    x=report.get('fresh',{}).get('x',{});boot=x.get('bootstrap',{});records=boot.get('records',[])
    expected=plan['analysis']['x']['bootstrap_replicates'];ids=[r.get('replicate') for r in records]
    expected_ids=set(range(expected));actual_ids=set(ids)
    bootstrap={'expected_replicates':expected,'returned_replicates':len(records),
        'missing_replicate_count':len(expected_ids-actual_ids),'unexpected_replicate_count':len(actual_ids-expected_ids),
        'duplicate_replicate_count':len(ids)-len(actual_ids),
        'declared_replicate_count_mismatch':boot.get('replicates')!=expected,
        'invalid_common_count':sum(r.get('common_valid') is not True for r in records),
        'invalid_alternative_count':sum(r.get('alternative_valid') is not True for r in records),
        'nonfinite_replicate_count':sum(not finite(r,('pc','nu','common_deviance','alternative_deviance','likelihood_ratio')) for r in records),
        'observed_model_failure_count':sum(x.get(name,{}).get('optimizer_success') is not True or not finite(x.get(name,{}),keys)
            for name,keys in [('common_model',('pc','nu','deviance')),('independent_size_model',('deviance',))])}
    flags['bootstrap_bound_hit_count']=sum(bool(r.get('bound_hits')) for r in records)
    flags['bootstrap_nonmonotone_count']=sum(r.get('monotone') is False for r in records)
    failed=bool(exceptions or any(fit_summary[k] for k in ('missing_count','unexpected_count','duplicate_count','optimizer_failure_count','nonfinite_fit_count'))
        or any(profiles[k] for k in ('missing_count','unexpected_count','duplicate_count','invalid_profile_count'))
        or any(v for k,v in bootstrap.items() if k not in ('expected_replicates','returned_replicates')))
    return {'status':'complete_with_numerical_failures' if failed else 'complete',
            'has_numerical_failures':failed,'exception_count':exceptions,'fits':fit_summary,
            'profiles':profiles,'bootstrap':bootstrap,'flags':flags,
            'interpretation':'Completion describes numerical validity of every planned calculation, not threshold acceptance. Bound hits, nonmonotonicity and small model goodness-of-fit p-values are reported separately and alone do not imply numerical failure.'}


def verify_numerical_sources(plan):
    reports={name:json.loads((PREVIOUS/name).read_text()) for name in ('x-report.json','z-report.json')}
    for name,digest in plan['prior_analysis_report_sha256'].items(): require(sha(PREVIOUS/name)==digest,'prior report hash mismatch')
    require(sha(PREVIOUS/'x_analysis.py')==reports['x-report.json']['provenance']['analysis_sha256'],'X helper source changed')
    require(sha(PREVIOUS/'z_models.py')==reports['z-report.json']['source_sha256'],'Z helper source changed')
    bound=reports['x-report.json']['provenance']['numerical_source_sha256']
    for name,digest in bound.items(): require(sha(ROOT/name)==digest,'likelihood engine source changed')
    paths=[Path(__file__),PREVIOUS/'x_analysis.py',PREVIOUS/'z_models.py']+[ROOT/name for name in bound]
    return {str(path.relative_to(ROOT)):sha(path) for path in paths}


def write_json(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def export_tables(output,report):
    def write(name,rows):
        if not rows:return
        fields=list(dict.fromkeys(key for row in rows for key in row))
        with (output/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
            writer.writerows({key:json.dumps(value) if isinstance(value,(dict,list)) else value for key,value in row.items()} for row in rows)
    for dataset in ('fresh','pooled_secondary'):
        for label,channel in report.get(dataset,{}).items():
            prefix=f'{dataset}-{label}'
            fields=('fit_id','channel','dataset','model','min_size','window','pc','nu','deviance','degrees_of_freedom',
                    'deviance_p_value','aicc','delta_aicc_same_subset','optimizer_success','bound_hits','monotone_increasing','diagnostic_three_sizes_only')
            write(prefix+'-fits.csv',[{key:f.get(key) for key in fields} for f in channel['fits']])
            write(prefix+'-predictions.csv',[dict(p,fit_id=f['fit_id']) for f in channel['fits'] for p in f['predictions']])
            for profile in channel['profiles']:write(prefix+'-profile-'+profile['fit_id']+'.csv',profile['evaluations'])
    x=report.get('fresh',{}).get('x',{})
    if 'bootstrap' in x:write('fresh-x-bootstrap.csv',x['bootstrap']['records'])
    if 'primary_crossing' in x:write('fresh-x-primary-contrasts.csv',x['primary_crossing']['contrasts'])
    if 'secondary_crossing' in x:
        write('fresh-x-cells.csv',x['secondary_crossing']['cells']);write('fresh-x-pair-bands.csv',x['secondary_crossing']['bands'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-dir',type=Path,default=HERE)
    parser.add_argument('--old-dir',type=Path,default=OLD)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--profiles-grid',type=int,default=301)
    parser.add_argument('--validate-only',action='store_true')
    args=parser.parse_args();require(args.profiles_grid>=31,'profile grid must have at least 31 points')
    start=time.time();study=args.study_dir;plan_path=study/'selection.json';plan=read_plan(study)
    require(plan['frozen_before_production_at_unix']<start,'plan not frozen before analysis')
    numerical_sources=verify_numerical_sources(plan)
    for name,digest in plan['simulation_source_sha256'].items(): require(sha(ROOT/name)==digest,'simulation source differs from frozen selection')
    old,old_provenance=load_old_production(plan,args.old_dir)
    protocols={r['protocol']['channel']:r['protocol'] for r in old_provenance['shards']}
    fresh=[];shards=[]
    for job in plan['jobs']:
        rows,record=read_shard(study/job['name'],job,plan['simulation_source_sha256'],protocols[job['channel']],plan['frozen_before_production_at_unix'])
        fresh.extend(rows);shards.append(record)
    fresh=ordered(fresh);validation=validate_fresh_points(fresh,plan)
    excluded,excluded_provenance=collect_excluded_streams(study,args.old_dir,plan)
    validation['independence']=assert_independent(fresh,old,excluded)
    pooled,pool_info=pool_points(fresh,old);validation['pooling']=pool_info
    output=args.output_dir or study/'analysis';output.mkdir(parents=True,exist_ok=True)
    snapshot={'schema_version':1,'fresh_points':fresh,'pooled_points':pooled,'validation':validation,
              'selection_sha256':sha(plan_path),'fresh_shards':shards,'old_production':old_provenance,'excluded_streams':excluded_provenance}
    write_json(output/'input.json',snapshot)
    if args.validate_only:
        print(json.dumps(validation,indent=2));return
    report={'schema_version':1,'created_at_unix':start,'status':'running','validation':validation,'fresh':{},'pooled_secondary':{},
            'provenance':{'selection_sha256':sha(plan_path),'input_sha256':sha(output/'input.json'),
                          'numerical_sources_sha256':numerical_sources,'simulation_sources_sha256':plan['simulation_source_sha256'],
                          'python':platform.python_version(),'python_executable':sys.executable,'numpy':np.__version__,'scipy':scipy.__version__,
                          'thread_environment':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}},
            'interpretation':'Fresh fixed-design estimates are primary; old-production same-cell pooling is secondary. Conditional model intervals omit finite-size model uncertainty. No automatic asymptotic-threshold acceptance.',
            'aicc_interpretation':Z.AICC_NOTE,
            'correction_normalization_note':'Each subset uses its geometric-mean L_reference. Correction amplitudes across subsets cannot be compared directly.'}
    def save(): write_json(output/'report.json',report);export_tables(output,report)
    for dataset,points in (('fresh',fresh),('pooled_secondary',pooled)):
        for channel,label in (('x_only','x'),('z_only','z')):
            report[dataset][label]=fit_channel(points,plan,channel,dataset);save()
    x=report['fresh']['x'];xp=x['points']
    x['primary_crossing']=anchor_assessment(fresh,plan);x['secondary_crossing']=secondary_crossing(xp);save()
    try:
        x.update(run_bootstrap(xp,plan))
    except (ArithmeticError,RuntimeError,ValueError) as error:
        x['failures'].append({'stage':'bootstrap','error':str(error)})
    save()
    z=report['fresh']['z'];z['profile_schedule']=[list(v) for v in profile_schedule(plan)]
    for minimum,window,model in profile_schedule(plan):
        identity=f'Lmin{minimum}-{window}-{model}'
        try:
            fit=next((f for f in z['fits'] if f['fit_id']==identity),None)
            require(fit is not None and fit['optimizer_success'],'profile free fit unavailable or unconverged')
            profile=profile_for_fit(select_data(fresh,plan,'z_only',minimum,window),model,fit,args.profiles_grid)
            z['profiles'].append(profile)
            print('PROFILE',identity,profile['regions'],flush=True)
        except (ArithmeticError,RuntimeError,ValueError) as error:
            z['failures'].append({'stage':'profile','fit_id':identity,'error':str(error)})
        save()
    z['profile_count']=len(z['profiles']);z['expected_profile_count']=len(profile_schedule(plan))
    report['numerical_summary']=summarize_numerics(report,plan)
    report['status']=report['numerical_summary']['status']
    report['completed_at_unix']=time.time();report['elapsed_seconds']=time.time()-start;save()
    print(report['status'],flush=True)


if __name__=='__main__': main()
