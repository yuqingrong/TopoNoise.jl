#!/usr/bin/env python3
"""Exploratory X uncertainty analysis of a frozen, validated data snapshot."""
import csv
import hashlib
import importlib.metadata
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from scipy.special import expit,xlogy
from scipy.stats import beta

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ENGINE_PATH=ROOT/'examples/fit_toric3d_h6_threshold.py'
SPEC=importlib.util.spec_from_file_location('unchanged_h6_x_engine',ENGINE_PATH)
H6=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(H6)
B=1000
SEED=2026101601
ALPHA=.05
GLOBAL_CHECK_STRIDE=25
GLOBAL_DEVIANCE_TOLERANCE=1e-5


def require(condition,message):
    if not condition:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clopper_pearson(k,n,alpha=.05):
    require(type(k) is int and type(n) is int and n>0 and 0<=k<=n,'invalid binomial counts')
    require(math.isfinite(alpha) and 0<alpha<1,'invalid alpha')
    return (0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),
            1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k)))


def merge_components(components):
    result=[]
    for low,high in sorted(components):
        require(low<=high,'reversed interval')
        if result and low<=result[-1][1]+16*np.finfo(float).eps*max(1.,abs(low)):
            result[-1][1]=max(result[-1][1],high)
        else:result.append([float(low),float(high)])
    return result


def compatible_components(rates,lower,upper):
    """All p with linearly interpolated lower(p)<=0<=upper(p), including points."""
    require(len(rates)==len(lower)==len(upper)>=2,'invalid band lengths')
    require(all(math.isfinite(v) for values in (rates,lower,upper) for v in values),'nonfinite band')
    require(all(a<b for a,b in zip(rates,rates[1:])),'rates must increase')
    require(all(a<=b for a,b in zip(lower,upper)),'reversed band')
    components=[]
    for i,(left,right) in enumerate(zip(rates,rates[1:])):
        lo,hi=0.,1.
        for a,b in ((lower[i],lower[i+1]),(-upper[i],-upper[i+1])):
            slope=b-a
            if slope==0:
                if a>0:lo,hi=1.,0.;break
            elif slope>0:hi=min(hi,-a/slope)
            else:lo=max(lo,-a/slope)
        if lo<=hi:
            components.append([left+(right-left)*lo,left+(right-left)*hi])
    return merge_components(components)


def intersect_components(first,second):
    return merge_components([[max(a,c),min(b,d)] for a,b in first for c,d in second
                             if max(a,c)<=min(b,d)])


def region_description(components,rates):
    return dict(components=components,empty=not components,
                touches_grid_lower=bool(components and components[0][0]<=rates[0]),
                touches_grid_upper=bool(components and components[-1][1]>=rates[-1]),
                covers_full_grid=components==[[rates[0],rates[-1]]])


def crossing_assessment(points):
    sizes=sorted({r['steps'] for r in points});rates=sorted({r['p'] for r in points})
    require(len(points)==len(sizes)*len(rates),'complete rectangular grid required')
    cells=[];lookup={}
    for r in points:
        lo,hi=clopper_pearson(r['failures'],r['shots'],ALPHA/len(points))
        cell=dict(L=r['steps'],p=r['p'],shots=r['shots'],failures=r['failures'],
                  rate=r['failures']/r['shots'],ci_low=lo,ci_high=hi)
        require((cell['L'],cell['p']) not in lookup,'duplicate cell')
        cells.append(cell);lookup[cell['L'],cell['p']]=cell
    pair_results=[];bands=[]
    adjacent=set(zip(sizes[:-1],sizes[1:]))
    for small,large in itertools.combinations(sizes,2):
        pair_rows=[]
        for p in rates:
            a,b=lookup[small,p],lookup[large,p]
            difference=b['rate']-a['rate']
            se=math.sqrt(a['rate']*(1-a['rate'])/a['shots']+b['rate']*(1-b['rate'])/b['shots'])
            entry=dict(L_small=small,L_large=large,adjacent=(small,large) in adjacent,p=p,
                difference=difference,standard_error=se,z_score=difference/se if se else None,
                rate_small=a['rate'],rate_large=b['rate'],
                ci_low=b['ci_low']-a['ci_high'],ci_high=b['ci_high']-a['ci_low'])
            pair_rows.append(entry);bands.append(entry)
        components=compatible_components(rates,[r['ci_low'] for r in pair_rows],[r['ci_high'] for r in pair_rows])
        pair_results.append(dict(sizes=[small,large],adjacent=(small,large) in adjacent,
                                 **region_description(components,rates)))
    common_adjacent=[[rates[0],rates[-1]]];common_all=[[rates[0],rates[-1]]]
    for pair in pair_results:
        common_all=intersect_components(common_all,pair['components'])
        if pair['adjacent']:common_adjacent=intersect_components(common_adjacent,pair['components'])
    return dict(method='equal-tail exact Clopper-Pearson; Bonferroni over all 35 binomial cells',
        family_alpha=ALPHA,cell_alpha=ALPHA/len(points),number_cells=len(points),cells=cells,
        bands=bands,pairs=pair_results,primary_adjacent=region_description(common_adjacent,rates),
        secondary_all_pairs=region_description(common_all,rates),
        secondary_note='All ten size pairs use the same simultaneous cell coverage event; secondary analysis requested after seeing wide adjacent-only bands, not a replacement for the frozen primary result.',
        interpretation='Zero-compatible regions are conditional on piecewise-linear interpolation of probability bounds. Adjacent-pair compatibility does not guarantee a common ordinate for all sizes. The all-pair secondary intersection is equivalent to max_L lower_L(p)<=min_L upper_L(p), hence all-size common-crossing compatibility. Neither establishes existence, uniqueness, or an asymptotic threshold; edge-touching regions are truncated by the sampled grid.')


def logistic_objective(coefficients,design,counts,shots):
    eta=design@coefficients;q=expit(eta);total=float(np.sum(shots))
    loss=float(np.sum(shots*np.logaddexp(0,eta)-counts*eta)/total)
    gradient=design.T@(shots*q-counts)/total
    hessian=design.T@((shots*q*(1-q))[:,None]*design)/total
    return loss,gradient,hessian


def binomial_deviance(counts,shots,probabilities):
    observed=counts/shots;safe=np.clip(probabilities,1e-14,1-1e-14)
    value=float(2*np.sum(xlogy(counts,observed/safe)+xlogy(shots-counts,(1-observed)/(1-safe))))
    require(value>=-1e-6 and math.isfinite(value),'invalid deviance')
    return max(0.,value)


def mapped_coefficients(common,L):
    scale=(L/common['L_reference'])**(1/common['nu'])
    center=(common['pc']-common['p_min'])/common['p_scale']
    a0,a1,a2=common['coefficients'];A=scale/2;C=scale*(.5-center)
    return np.array([a0+a1*C+a2*C*C,A*(a1+2*a2*C),a2*A*A])


def mapped_common_probabilities(points,common):
    answer=[]
    for r in points:
        t=2*(r['p']-common['p_min'])/common['p_scale']-1
        answer.append(float(expit(mapped_coefficients(common,r['steps'])@np.array([1,t,t*t]))))
    return answer


def fit_independent(points,common=None):
    """Unconstrained convex quadratic logits by size, initialized in the null model."""
    rates=sorted({r['p'] for r in points});pmin,pmax=rates[0],rates[-1]
    all_prob=np.empty(len(points));curves=[]
    for L in sorted({r['steps'] for r in points}):
        indices=[i for i,r in enumerate(points) if r['steps']==L]
        rows=[points[i] for i in indices]
        t=np.array([2*(r['p']-pmin)/(pmax-pmin)-1 for r in rows])
        design=np.column_stack([np.ones(len(t)),t,t*t])
        require(np.linalg.matrix_rank(design)==3,'independent quadratic needs three distinct rates')
        counts=np.array([r['failures'] for r in rows],float);shots=np.array([r['shots'] for r in rows],float)
        empirical=float(np.clip(counts.sum()/shots.sum(),1e-6,1-1e-6))
        coefficients=mapped_coefficients(common,L) if common else np.array([math.log(empirical/(1-empirical)),0.,0.])
        converged=False;reason='maximum iterations'
        for iteration in range(100):
            loss,gradient,hessian=logistic_objective(coefficients,design,counts,shots)
            if np.max(np.abs(gradient))<=1e-9:
                converged=True;reason='gradient tolerance';break
            try:direction=np.linalg.solve(hessian,gradient)
            except np.linalg.LinAlgError:reason='singular Hessian';break
            improvement=float(gradient@direction)
            if not math.isfinite(improvement) or improvement<=0:reason='invalid Newton direction';break
            step=1.
            while step>=2**-35:
                trial=coefficients-step*direction
                trial_loss=logistic_objective(trial,design,counts,shots)[0]
                if trial_loss<=loss-1e-4*step*improvement+1e-15:break
                step*=.5
            if step<2**-35:reason='line search failed';break
            coefficients=trial
        loss,gradient,_=logistic_objective(coefficients,design,counts,shots)
        if np.max(np.abs(gradient))<=1e-9:converged=True
        probability=expit(design@coefficients);all_prob[indices]=probability
        curves.append(dict(L=L,coefficients=coefficients.tolist(),optimizer_success=converged,
                           iterations=iteration+1,reason=reason,gradient_max=float(np.max(np.abs(gradient))),
                           deviance=binomial_deviance(counts,shots,probability)))
    counts=np.array([r['failures'] for r in points],float);shots=np.array([r['shots'] for r in points],float)
    deviance=binomial_deviance(counts,shots,all_prob)
    return dict(model='independent unconstrained quadratic-logit curve for each size',
        coordinate='t=2*(p-p_min)/(p_max-p_min)-1',p_min=pmin,p_max=pmax,
        parameter_count=3*len(curves),degrees_of_freedom=len(points)-3*len(curves),
        optimizer_success=all(c['optimizer_success'] for c in curves),deviance=deviance,
        fitted_probabilities=all_prob.tolist(),curves=curves)


def fit_common(points,warm_start=None):
    rates=sorted({r['p'] for r in points})
    return H6.ENGINE._candidate_fit(points,'x_only','none',p_min=rates[0],p_max=rates[-1],degree=2,
        nu_bounds=H6.NU_BOUNDS,goodness_alpha=.01,correction=False,warm_start=warm_start)


def common_usable(candidate):
    return bool(candidate and candidate.get('optimizer_success') and
        all(math.isfinite(candidate[k]) for k in ('pc','nu','deviance')) and
        all(math.isfinite(p) and 0<p<1 for p in candidate['fitted_probabilities']))


def checked_common(points,warm_start,global_check=False):
    errors=[];warm=None;global_fit=None
    try:warm=fit_common(points,warm_start)
    except (ValueError,RuntimeError,ArithmeticError) as error:errors.append('warm: '+str(error))
    fallback=not common_usable(warm)
    if global_check or fallback:
        try:global_fit=fit_common(points)
        except (ValueError,RuntimeError,ArithmeticError) as error:errors.append('global: '+str(error))
    valid=[r for r in (warm,global_fit) if common_usable(r)]
    best=min(valid,key=lambda r:r['deviance']) if valid else None
    improvement=(warm['deviance']-global_fit['deviance']) if common_usable(warm) and common_usable(global_fit) else None
    return best,dict(global_checked=bool(global_check or fallback),fallback=fallback,
                     global_improvement=improvement,errors=errors)


def tail_summary(tail,valid,attempted):
    require(0<=tail<=valid<=attempted and attempted>0,'invalid bootstrap counts')
    failed=attempted-valid
    return dict(tail_count=tail,valid_replicates=valid,attempted_replicates=attempted,numerical_failures=failed,
        p_plus_one=(tail+1)/(valid+1) if valid else None,
        mc_binomial_interval=list(clopper_pearson(tail,valid,.05)) if valid else None,
        mc_standard_error=math.sqrt((tail/valid)*(1-tail/valid)/valid) if valid else None,
        all_attempt_failure_bounds=[(tail+1)/(attempted+1),(tail+failed+1)/(attempted+1)],
        interpretation='Plug-in parametric bootstrap calibration, not an exact finite-sample test. Plus-one tail probability uses numerically successful fits; exact 95% binomial interval describes Monte Carlo tail-frequency uncertainty, not physical parameter uncertainty. Failure bounds assign all failed refits to either tail outcome.')


def bootstrap(points,common,alternative):
    rng=np.random.default_rng(SEED);shots=np.array([r['shots'] for r in points],int)
    simulated=rng.binomial(shots,np.array(common['fitted_probabilities']),size=(B,len(points)))
    observed_lr=max(0.,common['deviance']-alternative['deviance'])
    records=[];max_improvement=0.
    for index,counts in enumerate(simulated):
        rows=[dict(r,failures=int(k),none_failures=int(k)) for r,k in zip(points,counts)]
        fitted,validation=checked_common(rows,common['parameters'],global_check=index<5 or index%GLOBAL_CHECK_STRIDE==0)
        if validation['global_improvement'] is not None:max_improvement=max(max_improvement,validation['global_improvement'])
        entry=dict(replicate=index,counts=counts.tolist(),common_fit=fitted,common_validation=validation)
        records.append(entry)
        if (index+1)%100==0:print(f'X common bootstrap {index+1}/{B}',flush=True)
    global_recheck=max_improvement>GLOBAL_DEVIANCE_TOLERANCE
    if global_recheck:
        print('Global-validation discrepancy: globally rechecking every common bootstrap fit.',flush=True)
        for entry in records:
            rows=[dict(r,failures=int(k),none_failures=int(k)) for r,k in zip(points,entry['counts'])]
            fitted,validation=checked_common(rows,common['parameters'],global_check=True)
            entry.update(common_fit=fitted,common_validation=validation)
    output=[]
    for entry in records:
        fitted=entry['common_fit'];index=entry['replicate']
        rows=[dict(r,failures=int(k),none_failures=int(k)) for r,k in zip(points,entry['counts'])]
        record=dict(replicate=index,simulated_failures=entry['counts'],common_valid=common_usable(fitted),alternative_valid=False,
            common_deviance=None,alternative_deviance=None,likelihood_ratio=None,pc=None,nu=None,
            raw_likelihood_ratio=None,
            bound_hits=[],monotone=None,old_crossing_rejected=None,errors=entry['common_validation']['errors'],
            global_checked=entry['common_validation']['global_checked'],fallback=entry['common_validation']['fallback'],
            global_improvement=entry['common_validation']['global_improvement'])
        if record['common_valid']:
            crossings,clusters=H6.crossing_clusters(rows)
            rejected=not (len(clusters)==1 and any(a<=fitted['pc']<=b for a,b in clusters))
            record.update(common_deviance=fitted['deviance'],pc=fitted['pc'],nu=fitted['nu'],
                bound_hits=fitted['bound_hits'],monotone=fitted['monotone_increasing'],old_crossing_rejected=rejected)
        try:
            independent=fit_independent(rows,fitted if record['common_valid'] else None)
            require(independent['optimizer_success'],'alternative optimizer did not converge')
            record.update(alternative_valid=True,alternative_deviance=independent['deviance'])
            if record['common_valid']:
                require(independent['deviance']<=fitted['deviance']+1e-6,'nested alternative has larger deviance')
                record.update(raw_likelihood_ratio=fitted['deviance']-independent['deviance'],
                              likelihood_ratio=max(0.,fitted['deviance']-independent['deviance']))
        except (ValueError,RuntimeError,ArithmeticError) as error:
            record['alternative_valid']=False;record['errors'].append('alternative: '+str(error))
        output.append(record)
    common_valid=[r for r in output if r['common_valid']]
    joint_valid=[r for r in output if r['common_valid'] and r['alternative_valid']]
    gof=tail_summary(sum(r['common_deviance']>=common['deviance']-1e-10 for r in common_valid),len(common_valid),B)
    lr=tail_summary(sum(r['likelihood_ratio']>=observed_lr-1e-10 for r in joint_valid),len(joint_valid),B)
    old_count=sum(r['old_crossing_rejected'] for r in common_valid)
    old_rule=tail_summary(old_count,len(common_valid),B)
    old_rule['rejection_frequency']=old_count/len(common_valid) if common_valid else None
    old_rule['calibration']='Rejection rate of the original raw-sign guard under the fitted common-scaling null, refitting pc in every replicate; not a new circuit experiment.'
    intervals={key:np.quantile([r[key] for r in common_valid],[.025,.975]).tolist() if len(common_valid)>=2 else None for key in ('pc','nu')}
    return dict(replicates=B,seed=SEED,null='fitted common quadratic-logit scaling model',
        goodness_of_fit=gof,likelihood_ratio=lr,old_raw_crossing_rule=old_rule,
        conditional_percentile_intervals=intervals,interval_population='all converged finite common refits, including boundary/nonmonotone flags',
        bound_hit_replicates=sum(bool(r['bound_hits']) for r in common_valid),
        nonmonotone_replicates=sum(not r['monotone'] for r in common_valid),
        optimization_validation=dict(global_check_stride=GLOBAL_CHECK_STRIDE,first_five_checked=True,
            validation_indices=[r['replicate'] for r in output if r['global_checked']],
            maximum_sample_global_improvement=max_improvement,tolerance=GLOBAL_DEVIANCE_TOLERANCE,
            all_replicates_globally_rechecked=global_recheck),records=output)


def write_csv(path,rows):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    start=time.time();snapshot=json.loads((HERE/'input.json').read_text());freeze=json.loads((HERE/'plan-freeze.json').read_text())
    require(freeze['input_sha256']==sha(HERE/'input.json'),'frozen input hash mismatch')
    require(freeze['plan_sha256']==sha(HERE/'analysis-plan.md'),'frozen plan hash mismatch')
    require(freeze['fixed_before_followup_fits_at_unix']<start,'plan not frozen before fit')
    for filename,digest in snapshot['original_selection_validation']['numerical_source_sha256'].items():
        require(sha(ROOT/filename)==digest,'original numerical implementation changed')
    points=[r for r in snapshot['points'] if r['channel']=='x_only']
    require(len(points)==35 and sum(r['shots'] for r in points)==35000,'unexpected X data size')
    crossing=crossing_assessment(points)
    common=fit_common(points);require(common_usable(common),'observed common fit failed')
    original=snapshot['original_fits']['x_only']['candidate']
    reproduction={k:common[k]-original[k] for k in ('pc','nu','deviance')}
    require(all(abs(v)<1e-7 for v in reproduction.values()),'original common fit not reproduced')
    independent=fit_independent(points,common);require(independent['optimizer_success'],'observed alternative fit failed')
    require(independent['deviance']<=common['deviance']+1e-6,'observed models not nested numerically')
    np.testing.assert_allclose(mapped_common_probabilities(points,common),common['fitted_probabilities'],atol=1e-13)
    bootstrap_result=bootstrap(points,common,independent)
    numerical_files=[ENGINE_PATH,H6.ENGINE_PATH]
    report=dict(channel='x_only',exploratory=True,uses_new_circuit_samples=False,number_points=35,number_fault_records=35000,
        crossing=crossing,common_model=common,independent_size_model=independent,
        observed_likelihood_ratio=max(0.,common['deviance']-independent['deviance']),
        bootstrap=bootstrap_result,original_fit_reproduction_differences=reproduction,
        method_deviations=[],secondary_addendum=crossing['secondary_note'],
        provenance=dict(input_sha256=sha(HERE/'input.json'),plan_sha256=sha(HERE/'analysis-plan.md'),
            plan_freeze_sha256=sha(HERE/'plan-freeze.json'),analysis_sha256=sha(__file__),
            secondary_addendum_sha256=sha(HERE/'analysis-addendum.md'),
            numerical_source_sha256={str(p.relative_to(ROOT)):sha(p) for p in numerical_files},
            source_report_sha256=snapshot['source_report_sha256'],source_manifest_sha256=snapshot['source_manifest_sha256'],
            packages={name:importlib.metadata.version(name) for name in ('numpy','scipy')},
            python_executable=sys.executable,python_version=sys.version,
            thread_environment={key:os.environ.get(key) for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
            started_at_unix=start,elapsed_seconds=time.time()-start),
        interpretation='Exploratory reanalysis of previously observed data. Simultaneous crossing compatibility and null-model calibration do not independently establish an asymptotic threshold; original rejection remains unchanged.')
    (HERE/'x-report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    write_csv(HERE/'x-cells.csv',crossing['cells']);write_csv(HERE/'x-pair-bands.csv',crossing['bands'])
    regions=[]
    for pair in crossing['pairs']:
        for index,component in enumerate(pair['components'] or [[None,None]]):
            regions.append(dict(L_small=pair['sizes'][0],L_large=pair['sizes'][1],adjacent=pair['adjacent'],
                component=index,low=component[0],high=component[1],empty=pair['empty'],
                touches_grid_lower=pair['touches_grid_lower'],touches_grid_upper=pair['touches_grid_upper']))
    write_csv(HERE/'x-crossing-regions.csv',regions)
    samples=[dict(r,bound_hits=json.dumps(r['bound_hits']),errors=json.dumps(r['errors'])) for r in bootstrap_result['records']]
    write_csv(HERE/'x-bootstrap.csv',samples)
    print(json.dumps(dict(primary_crossing=crossing['primary_adjacent'],secondary_crossing=crossing['secondary_all_pairs'],
        common_deviance=common['deviance'],alternative_deviance=independent['deviance'],
        likelihood_ratio=report['observed_likelihood_ratio'],gof=bootstrap_result['goodness_of_fit'],
        lr=bootstrap_result['likelihood_ratio'],old_crossing=bootstrap_result['old_raw_crossing_rule'],
        elapsed_seconds=report['provenance']['elapsed_seconds']),indent=2),flush=True)


if __name__=='__main__':main()
