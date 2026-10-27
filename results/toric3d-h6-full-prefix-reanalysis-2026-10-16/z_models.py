#!/usr/bin/env python3
"""Fixed exploratory Z finite-size sensitivity fits and conditional pc profiles."""
import argparse
import csv
import hashlib
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize, brentq
from scipy.special import expit, gammaln, xlogy
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
MODELS = ('quadratic', 'cubic', 'corr_w0p5', 'corr_w1', 'corr_w2', 'corr_linear_w1')
WINDOWS = {'full': (.025, .04), 'narrow': (.0275, .0375)}
OMEGAS = {'corr_w0p5': .5, 'corr_w1': 1., 'corr_w2': 2., 'corr_linear_w1': 1.}
AICC_NOTE = ('Descriptive heuristic with independent binomial cells as sample size; '
             'compare only identical data subsets. Neither AICc nor improved goodness establishes a threshold.')


def select_data(points, min_size, window):
    low, high = WINDOWS[window]
    rows = sorted((row for row in points if row['channel'] == 'z_only' and
                   row['steps'] >= min_size and low <= row['p'] <= high),
                  key=lambda row:(row['steps'], row['p']))
    sizes = np.array([row['steps'] for row in rows], dtype=float)
    reference = float(np.exp(np.mean(np.log(np.unique(sizes)))))
    return {'sizes': sizes, 'p': np.array([row['p'] for row in rows]),
            'shots': np.array([row['shots'] for row in rows], dtype=float),
            'failures': np.array([row['full_prefix_failures'] for row in rows], dtype=float),
            'L_reference': reference, 'relative_log': np.log(sizes/reference),
            'p_min': low, 'p_max': high, 'p_scale': high-low,
            'min_size': min_size, 'window': window}


def parameter_count(model):
    return 7 if model == 'corr_linear_w1' else 5 if model == 'quadratic' else 6


def parameter_names(model):
    names = ['pc', 'nu', 'a0', 'a1', 'a2']
    if model == 'cubic':
        names += ['a3']
    elif model == 'corr_linear_w1':
        names += ['c0', 'c1']
    elif model in OMEGAS:
        names += ['c']
    return names


def parameter_bounds(model):
    bounds = [(0.,1.), (math.log(.3),math.log(5.)), (-15.,15.), (-30.,30.), (-30.,30.)]
    if model == 'cubic':
        bounds += [(-30.,30.)]
    elif model == 'corr_linear_w1':
        bounds += [(-15.,15.), (-15.,15.)]
    elif model in OMEGAS:
        bounds += [(-15.,15.)]
    return bounds


def eta_jacobian(theta, data, model):
    center, lognu = theta[:2]
    nu = math.exp(lognu)
    stretch = np.exp(data['relative_log']/nu)
    x = ((data['p']-data['p_min'])/data['p_scale']-center)*stretch
    degree = 3 if model == 'cubic' else 2
    coefficients = theta[2:degree+3]
    eta = np.polynomial.polynomial.polyval(x, coefficients)
    derivative = np.polynomial.polynomial.polyval(x, np.arange(1,degree+1)*coefficients[1:])
    basis = (data['L_reference']/data['sizes'])**OMEGAS[model] if model in OMEGAS else None
    if model == 'corr_linear_w1':
        eta = eta + (theta[-2]+theta[-1]*x)*basis
        derivative = derivative + theta[-1]*basis
    elif model in OMEGAS:
        eta = eta + theta[-1]*basis
    jacobian = [-stretch*derivative, -x*data['relative_log']/nu*derivative]
    jacobian.extend(x**power for power in range(degree+1))
    if model == 'corr_linear_w1':
        jacobian.extend([basis, x*basis])
    elif model in OMEGAS:
        jacobian.append(basis)
    return eta, np.asarray(jacobian).T


def objective(theta, data, model):
    eta, jacobian = eta_jacobian(theta, data, model)
    n, k = data['shots'], data['failures']
    residual = n*expit(eta)-k
    total = np.sum(n)
    return float(np.sum(n*np.logaddexp(0.,eta)-k*eta)/total), jacobian.T@residual/total


def _linear_start(data, model, center, nu):
    theta = np.zeros(parameter_count(model)); theta[:2] = center, math.log(nu)
    _, jac = eta_jacobian(theta, data, model)
    empirical = (data['failures']+.5)/(data['shots']+1.)
    target = np.log(empirical/(1-empirical))
    weight = np.sqrt(data['shots']*empirical*(1-empirical))
    theta[2:] = np.linalg.lstsq(jac[:,2:]*weight[:,None],target*weight,rcond=None)[0]
    return np.clip(theta,*np.asarray(parameter_bounds(model)).T)


def _starts(data, model, fixed_center=None, warm_starts=()):
    bounds = np.asarray(parameter_bounds(model))
    starts = [np.clip(np.asarray(value,dtype=float),bounds[:,0],bounds[:,1]) for value in warm_starts]
    centers = (.1,.25,.5,.75,.9) if fixed_center is None else (fixed_center,)
    nus = (.35,.6,1.2,2.5,4.5) if fixed_center is None else (.3,.6,1.2,2.5,5.)
    starts += [_linear_start(data,model,center,nu) for center in centers for nu in nus]
    if fixed_center is None:
        fraction = np.clip(np.sum(data['failures'])/np.sum(data['shots']),.001,.999)
        intercept = math.log(fraction/(1-fraction))
        for center in (.25,.5,.75):
            for nu in (.6,1.2,2.5):
                theta = np.zeros(parameter_count(model)); theta[:4] = center,math.log(nu),intercept,2.
                starts.append(theta)
        if model in OMEGAS:
            for sign in (-1.,1.):
                theta = _linear_start(data,model,.5,1.2)
                theta[5] = sign
                starts.append(theta)
    if fixed_center is not None:
        for start in starts:
            start[0] = fixed_center
    return starts


def _optimize(fun, bounds, starts):
    results = [minimize(fun,start,method='L-BFGS-B',jac=True,bounds=bounds,
                       options={'maxiter':1500,'ftol':1e-12,'gtol':1e-8}) for start in starts]
    finite = [result for result in results if np.isfinite(result.fun) and np.all(np.isfinite(result.x))]
    if not finite:
        raise RuntimeError('All likelihood optimizations returned nonfinite values')
    best = min(finite,key=lambda result:result.fun)
    converged = [result for result in finite if result.success]
    if converged:
        converged_best = min(converged,key=lambda result:result.fun)
        if converged_best.fun <= best.fun+1e-12:
            best = converged_best
    polish = minimize(fun,best.x,method='L-BFGS-B',jac=True,bounds=bounds,
                      options={'maxiter':3000,'ftol':1e-14,'gtol':1e-10})
    if np.isfinite(polish.fun) and ((polish.success and polish.fun <= best.fun+1e-14)
                                  or polish.fun < best.fun-1e-12):
        best = polish
    gradient = np.asarray(best.jac).copy()
    for index,(value,(low,high)) in enumerate(zip(best.x,bounds)):
        if value <= low+1e-10 and gradient[index] > 0 or value >= high-1e-10 and gradient[index] < 0:
            gradient[index] = 0.
    return best, {'starts':len(starts),'successful_starts':sum(bool(r.success) for r in results),
        'failed_starts':sum(not bool(r.success) for r in results),
        'nonfinite_starts':len(results)-len(finite),
        'polish_success':bool(polish.success),'polish_message':str(polish.message),
        'projected_gradient_infinity_norm':float(np.max(np.abs(gradient))),
        'finite_loss_range':[float(min(r.fun for r in finite)),float(max(r.fun for r in finite))],
        'start_results':[{'success':bool(r.success),'objective':float(r.fun) if np.isfinite(r.fun) else None,
                          'message':str(r.message)} for r in results]}


def monotonicity(theta, data, model):
    degree = 3 if model == 'cubic' else 2
    derivative = np.arange(1,degree+1)*np.asarray(theta[3:degree+3])
    by_size = []
    nu = math.exp(theta[1]); pc = data['p_min']+theta[0]*data['p_scale']
    for size in np.unique(data['sizes']):
        stretch = (size/data['L_reference'])**(1/nu)
        low,high = [(p-pc)/data['p_scale']*stretch for p in (data['p_min'],data['p_max'])]
        polynomial = derivative.copy()
        if model == 'corr_linear_w1':
            polynomial[0] += theta[-1]*data['L_reference']/size
        candidates = [low,high]
        second = np.arange(1,len(polynomial))*polynomial[1:]
        if len(second)>1:
            candidates += [float(root.real) for root in np.polynomial.polynomial.polyroots(second)
                           if abs(root.imag)<1e-10 and low<=root.real<=high]
        minimum = float(min(np.polynomial.polynomial.polyval(x,polynomial) for x in candidates))
        by_size.append({'size':int(size),'minimum_deta_dx':minimum,
                        'minimum_deta_dp':minimum*stretch/data['p_scale']})
    minimum = min(item['minimum_deta_dx'] for item in by_size)
    return {'monotone_increasing':bool(minimum>1e-7),'derivative_min':minimum,'by_size':by_size}


def _summarize(theta, data, model, success, message, optimizer):
    eta,_ = eta_jacobian(theta,data,model)
    fitted = expit(eta); safe = np.clip(fitted,1e-14,1-1e-14)
    n,k = data['shots'],data['failures']; observed = k/n
    nll = float(np.sum(n*np.logaddexp(0.,eta)-k*eta))
    deviance = float(2*np.sum(xlogy(k,observed/safe)+xlogy(n-k,(1-observed)/(1-safe))))
    loglik = float(np.sum(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)))-nll
    count = len(theta); cells = len(n); df = cells-count
    pearson = (observed-fitted)/np.sqrt(np.maximum(fitted*(1-fitted)/n,1e-30))
    cell_deviance = 2*(xlogy(k,observed/safe)+xlogy(n-k,(1-observed)/(1-safe)))
    aic = -2*loglik+2*count
    bounds = parameter_bounds(model)
    hits = [name for name,value,(low,high) in zip(parameter_names(model),theta,bounds)
            if min(value-low,high-value)<=.002*(high-low)]
    predictions = [{'size':int(size),'p':float(p),'shots':int(shots),'failures':int(failures),
                    'observed':float(obs),'predicted':float(pred),'pearson_residual':float(res),
                    'deviance_residual':float(np.sign(obs-pred)*np.sqrt(max(0.,dev)))}
                   for size,p,shots,failures,obs,pred,res,dev in
                   zip(data['sizes'],data['p'],n,k,observed,fitted,pearson,cell_deviance)]
    return {'model':model,'min_size':data['min_size'],'window':data['window'],
            'pc':float(data['p_min']+theta[0]*data['p_scale']),'nu':float(math.exp(theta[1])),
            'parameters':theta.tolist(),'parameter_names':parameter_names(model),
            'parameter_bounds':bounds,'optimization_parameter_names':['center','log_nu']+parameter_names(model)[2:],
            'physical_parameter_bounds':{'pc':[data['p_min'],data['p_max']],'nu':[.3,5.]},'coefficients':theta[2:6 if model=='cubic' else 5].tolist(),
            'correction_omega':OMEGAS.get(model),'correction_coefficients':theta[5:].tolist() if model in OMEGAS else [],
            'p_min':data['p_min'],'p_max':data['p_max'],'p_scale':data['p_scale'],'L_reference':data['L_reference'],
            'sizes':np.unique(data['sizes']).astype(int).tolist(),'rates':np.unique(data['p']).tolist(),
            'n_cells':cells,'n_parameters':count,'diagnostic_three_sizes_only':len(np.unique(data['sizes']))==3,
            'nll':nll,'log_likelihood':loglik,'deviance':deviance,'degrees_of_freedom':df,
            'deviance_p_value':float(chi2.sf(deviance,df)) if df>0 else None,
            'aic':aic,'aicc':aic+2*count*(count+1)/(cells-count-1) if cells>count+1 else None,
            'aicc_sample_size':cells,'aicc_interpretation':AICC_NOTE,
            'optimizer_success':bool(success),'optimizer_message':message,'optimizer':optimizer,
            'bound_hits':hits,**monotonicity(theta,data,model),
            'fitted_probabilities':fitted.tolist(),'standardized_residuals':pearson.tolist(),
            'largest_absolute_residual':float(np.max(np.abs(pearson))),'predictions':predictions,
            'interpretation':'Exploratory sensitivity candidate; no threshold acceptance claim.'}


def fit_model(data, model, warm_starts=()):
    best,diagnostic = _optimize(lambda theta:objective(theta,data,model),parameter_bounds(model),
                                _starts(data,model,warm_starts=warm_starts))
    return _summarize(best.x,data,model,best.success,str(best.message),diagnostic)


def fit_fixed_pc(data, model, pc, warm_starts=()):
    center = (pc-data['p_min'])/data['p_scale']
    if not -1e-12<=center<=1+1e-12:
        raise ValueError('Profile pc lies outside original parameter bounds')
    center = np.clip(center,0.,1.)
    def fun(nuisance):
        value,gradient = objective(np.r_[center,nuisance],data,model)
        return value,gradient[1:]
    starts = _starts(data,model,fixed_center=center,warm_starts=warm_starts)
    best,diagnostic = _optimize(fun,parameter_bounds(model)[1:],[start[1:] for start in starts])
    theta = np.r_[center,best.x]
    return {'pc':float(pc),'nll':float(best.fun*np.sum(data['shots'])),'parameters':theta.tolist(),
            'optimizer_success':bool(best.success),'optimizer_message':str(best.message),
            'failed_starts':diagnostic['failed_starts'],'nonfinite_starts':diagnostic['nonfinite_starts'],
            'projected_gradient_infinity_norm':diagnostic['projected_gradient_infinity_norm'],
            'polish_success':diagnostic['polish_success'],
            'bound_hits':[name for name,value,(low,high) in zip(parameter_names(model)[1:],theta[1:],parameter_bounds(model)[1:])
                          if min(value-low,high-value)<=.002*(high-low)]}


def connected_regions(grid, values, cutoff, evaluate):
    grid,values = np.asarray(grid),np.asarray(values)
    inside = np.isfinite(values)&(values<=cutoff)
    output = []; index = 0
    while index<len(grid):
        if not inside[index]:
            index += 1; continue
        start = end = index
        while end+1<len(grid) and inside[end+1]: end += 1
        low,high = float(grid[start]),float(grid[end])
        lower_gap = start>0 and not np.isfinite(values[start-1])
        upper_gap = end+1<len(grid) and not np.isfinite(values[end+1])
        refinement_failures = []
        if start>0 and not lower_gap:
            try:
                low = float(brentq(lambda pc:evaluate(pc)-cutoff,grid[start-1],grid[start],xtol=1e-12))
            except ValueError as error:
                lower_gap = True
                refinement_failures.append({'side':'lower','bracket':[float(grid[start-1]),float(grid[start])],'error':str(error)})
        if end+1<len(grid) and not upper_gap:
            try:
                high = float(brentq(lambda pc:evaluate(pc)-cutoff,grid[end],grid[end+1],xtol=1e-12))
            except ValueError as error:
                upper_gap = True
                refinement_failures.append({'side':'upper','bracket':[float(grid[end]),float(grid[end+1])],'error':str(error)})
        output.append({'lower':low,'upper':high,'lower_bound_truncated':start==0,
                       'upper_bound_truncated':end==len(grid)-1,
                       'lower_numerical_gap':bool(lower_gap),'upper_numerical_gap':bool(upper_gap),
                       'refinement_failures':refinement_failures})
        index = end+1
    return output


def profile_pc(data, model, fit, grid_count=301):
    cutoff = float(chi2.ppf(.95,1))
    cache = {}; failures = []
    def evaluate(pc):
        key = float(pc)
        if key not in cache:
            warm = [fit['parameters']]
            if cache:
                finite = [row for row in cache.values() if row['optimizer_success'] and row.get('parameters')]
                if finite:
                    warm.append(min(finite,key=lambda row:abs(row['pc']-key))['parameters'])
            try:
                point = fit_fixed_pc(data,model,key,warm)
                point['lr'] = 2*(point['nll']-fit['nll'])
                if not point['optimizer_success']:
                    failures.append({'pc':key,'message':point['optimizer_message']})
            except (ArithmeticError,RuntimeError,ValueError) as error:
                point = {'pc':key,'nll':None,'lr':None,'parameters':None,'optimizer_success':False,
                         'optimizer_message':str(error),'failed_starts':None,'nonfinite_starts':None,
                         'projected_gradient_infinity_norm':None}
                failures.append({'pc':key,'message':str(error)})
            if point['optimizer_success'] and point['nll'] < fit['nll']-2e-4:
                raise RuntimeError('Profile nuisance optimum improves the free fit; refit required before reporting intervals')
            cache[key] = point
        point = cache[key]
        return point['lr'] if point['optimizer_success'] and point['lr'] is not None else float('nan')
    grid = np.unique(np.r_[np.linspace(data['p_min'],data['p_max'],grid_count),fit['pc']])
    # Begin at the fitted optimum and expand outward to aid continuation.
    for pc in sorted(grid,key=lambda value:abs(value-fit['pc'])): evaluate(float(pc))
    values = [evaluate(float(pc)) for pc in grid]
    intervals = connected_regions(grid,values,cutoff,evaluate)
    finite = [point for point in cache.values() if point['optimizer_success'] and point['nll'] is not None]
    minimum = min(finite,key=lambda point:point['nll'])
    at_free = cache[float(fit['pc'])]
    agreement = abs(at_free['nll']-fit['nll'])<=2e-4 and minimum['nll']>=fit['nll']-2e-4
    return {'fit_id':fit['fit_id'],'model':model,'min_size':data['min_size'],'window':'full',
            'confidence_level':.95,'lr_cutoff':cutoff,'pc_bounds':[data['p_min'],data['p_max']],
            'grid_count':len(grid),'base_grid_spacing':data['p_scale']/(grid_count-1),
            'free_pc':fit['pc'],'free_nll':fit['nll'],'profile_at_free_pc':at_free,
            'minimum_profile_pc':minimum['pc'],'minimum_profile_nll':minimum['nll'],
            'minimum_agreement':bool(agreement),'minimum_nll_difference':minimum['nll']-fit['nll'],
            'regions':intervals,'numerical_failures':failures,
            'evaluations':sorted(cache.values(),key=lambda point:point['pc']),
            'polish_failure_count':sum(point.get('polish_success') is False for point in cache.values()),
            'nuisance_bound_evaluation_count':sum(bool(point.get('bound_hits')) for point in cache.values()),
            'interpretation':'Approximate conditional likelihood-ratio regions under this model and original bounds; no model uncertainty included.',
            'resolution_note':'All detected components at the saved grid resolution; cutoff endpoints refined by nuisance reoptimization and Brent roots. Finer-than-grid disconnected islands cannot be guaranteed absent.'}


def write_csv(path, rows, fields):
    with path.open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore'); writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--profiles-grid',type=int,default=301)
    args = parser.parse_args()
    snapshot = json.loads((HERE/'input.json').read_text()); freeze = json.loads((HERE/'plan-freeze.json').read_text())
    assert hashlib.sha256((HERE/'input.json').read_bytes()).hexdigest()==freeze['input_sha256']
    assert hashlib.sha256((HERE/'analysis-plan.md').read_bytes()).hexdigest()==freeze['plan_sha256']
    report = {'schema_version':1,'channel':'z_only','arm':'full_prefix','exploratory':True,
              'created_at_unix':time.time(),'input_sha256':freeze['input_sha256'],'plan_sha256':freeze['plan_sha256'],
              'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'numerically_equivalent_average_nll_tolerance':1e-12,
              'correction_normalization_note':'Correction coefficients use the geometric-mean L_reference of each subset; raw c values across subsets are not directly comparable.',
              'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
              'aicc_interpretation':AICC_NOTE,'normalization':'x=((p-pc)/(pmax-pmin))*(L/L_reference)^(1/nu); L_reference=geometric mean of included unique sizes.',
              'optimizer_policy':'Deterministic grid and weighted-logit linear starts; original L-BFGS-B tolerances, followed by tighter numerical polish. Retain failed-start counts and finite best candidate.',
              'deviations':['No scientific plan deviation. Linear correction c0 and c1 each retain the original correction coefficient bounds [-15,15]. Numerical polish tightens ftol to 1e-14 and gtol to 1e-10 after original-tolerance multistarts.'],
              'fits':[],'profiles':[],'failures':[],'acceptance':'No threshold accepted by this exploratory analysis; original acceptance policy remains unchanged.'}
    for minimum in (3,5,7):
        for window in WINDOWS:
            data = select_data(snapshot['points'],minimum,window)
            for model in MODELS:
                identity = f'Lmin{minimum}-{window}-{model}'
                try:
                    fit = fit_model(data,model); fit['fit_id']=identity
                    report['fits'].append(fit)
                    print(identity,'pc',fit['pc'],'nu',fit['nu'],'D',fit['deviance'],'p',fit['deviance_p_value'],'bounds',fit['bound_hits'],flush=True)
                except (ArithmeticError,RuntimeError,ValueError) as error:
                    report['failures'].append({'fit_id':identity,'error':str(error)})
                    print(identity,'FAIL',error,flush=True)
    for minimum in (3,5,7):
        for window in WINDOWS:
            subset = [fit for fit in report['fits'] if fit['min_size']==minimum and fit['window']==window]
            floor = min(fit['aicc'] for fit in subset if fit['aicc'] is not None)
            for fit in subset: fit['delta_aicc_same_subset']=fit['aicc']-floor
    (HERE/'z-report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    fields = ['fit_id','model','min_size','window','pc','nu','deviance','degrees_of_freedom','deviance_p_value',
              'log_likelihood','aic','aicc','delta_aicc_same_subset','aicc_sample_size','optimizer_success','bound_hits',
              'monotone_increasing','derivative_min','largest_absolute_residual','diagnostic_three_sizes_only']
    write_csv(HERE/'z-fits.csv',report['fits'],fields)
    predictions = [dict(row,fit_id=fit['fit_id']) for fit in report['fits'] for row in fit['predictions']]
    write_csv(HERE/'z-predictions.csv',predictions,['fit_id','size','p','shots','failures','observed','predicted','pearson_residual','deviance_residual'])
    for minimum in (3,5):
        data = select_data(snapshot['points'],minimum,'full')
        for model in ('quadratic','corr_w1'):
            fit = next(fit for fit in report['fits'] if fit['fit_id']==f'Lmin{minimum}-full-{model}')
            profile = profile_pc(data,model,fit,args.profiles_grid); report['profiles'].append(profile)
            print('PROFILE',fit['fit_id'],'regions',profile['regions'],'minimum agreement',profile['minimum_agreement'],flush=True)
            write_csv(HERE/f'z-profile-{fit["fit_id"]}.csv',profile['evaluations'],
                      ['pc','nll','lr','optimizer_success','optimizer_message','failed_starts','nonfinite_starts','polish_success','bound_hits','projected_gradient_infinity_norm','parameters'])
            (HERE/'z-report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    report['completed_at_unix']=time.time()
    report['fit_count']=len(report['fits']); report['profile_count']=len(report['profiles'])
    (HERE/'z-report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('COMPLETE',len(report['fits']),'fits',len(report['profiles']),'profiles',flush=True)


if __name__=='__main__': main()
