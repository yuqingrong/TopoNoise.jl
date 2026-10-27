#!/usr/bin/env python3
"""Measured paired H+6 X/Z full-prefix results with the existing guarded fitter."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys
import tomllib

import numpy as np
from scipy.special import expit

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'examples'))
OLD_PATH=ROOT/'results/toric3d-h6-z-one-step-2026-10-15/analyze.py'
SPEC=importlib.util.spec_from_file_location('existing_h6_guarded_analysis',OLD_PATH)
OLD=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OLD)
H6=OLD.H6
ARMS=('final_only','full_prefix')
CHANNELS=('x_only','z_only')
COLORS={3:'#7F7F7F',5:'#0072B2',7:'#E69F00',9:'#009E73',11:'#CC79A7'}


def require(condition,message):
    if not condition: raise ValueError(message)


def integer(value):
    require(not isinstance(value,bool),'Boolean is not a count')
    result=int(value)
    require(str(result)==str(value) or result==value,'noninteger count')
    return result


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def expected_protocol(channel):
    """Exact approved runner signature: complete recovery after transfers 1:T."""
    require(channel in CHANNELS,'unknown protocol channel')
    x=channel=='x_only'
    return dict(circuit='legacy H+6-CNOT plaquette transfer',channel=channel,
        geometry='T=L; Ly=Lz=L',arm_order=list(ARMS),
        initial_state='ideal encoded logical |00>' if x else 'noisy seed-H preparation of logical |++>',
        seed_preparation='ideal; zero seed fault locations' if x else '3 L^2 noisy seed-H outputs; no separate r=0 seed correction in either arm',
        final_only='no intermediate physical correction; exactly one final recovery and logical outcome',
        full_prefix='complete physical recovery on all available data slices Q0..Qr after EVERY transfer 1:T, including the last; complete-prefix syndrome zero; exactly one final recovery and logical outcome afterward',
        feedback_timing='after transfers 1:T inclusive; no separate r=0 correction',
        fault_locations='13 L^3 participating transfer outputs' if x else '3 L^2 seed-H outputs plus 13 L^3 participating transfer outputs',
        noise='independent Bernoulli-p faults in the selected channel at its recorded circuit locations',
        sampler='sample_toric3d_faults(p_x=p,p_z=0)' if x else 'sample_toric3d_diagonal_z_faults: seed first, then participating transfer outputs in schedule order',
        shared_final_decoder='family.prefixes[end]; full-plaquette BP+OSD with declared p prior; syndrome-validated recovery' if x else 'family.prefixes[end]; unit-weight data-edge PyMatching MWPM; syndrome-validated recovery',
        decoder_settings=dict(bp_method='minimum_sum',max_iter=50,ms_scaling_factor=1.,omp_thread_count=1,
                              osd_method='OSD_CS',osd_order=2,schedule='parallel') if x else dict(weights=1.,merge_strategy='smallest-weight'),
        z_weights='unit data-edge weights',failure='either residual logical parity after final recovery',
        joint_histogram_encoding='mask = final_only + 2*full_prefix; each term is its state-failure Boolean',
        point_order='sizes outer loop, rates inner loop; both strictly increasing',
        point_seed_generation='MersenneTwister(integer base seed); repeated scalar rand(rng,UInt64); discard duplicate values; 16 lowercase hexadecimal digits',
        point_rng='Random.Xoshiro(point_seed)',
        resume='recreate point Xoshiro and sample/discard completed records without decoding; rebuild raw CSV from authoritative checkpoint')


def validate_runtime(runtime,packages):
    require(runtime.get('julia_threads')==runtime.get('blas_threads')==1,'unexpected thread count')
    environment=runtime.get('environment',{})
    require(all(environment.get(key)=='1' for key in
            ('JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')),
            'unexpected thread environment')
    require(environment.get('JULIA_CONDAPKG_BACKEND')=='Null','unexpected PythonCall backend')
    require(runtime.get('packages')==packages,'runtime package mismatch')
    require(isinstance(runtime.get('python_executable'),str) and
            isinstance(environment.get('JULIA_PYTHONCALL_EXE'),str) and
            Path(runtime['python_executable']).resolve()==Path(environment['JULIA_PYTHONCALL_EXE']).resolve()==Path(sys.executable).resolve(),
            'Python executable mismatch')
    # Embedded libpython can have a different build timestamp than its executable.
    require(isinstance(runtime.get('python_version'),str) and
            runtime['python_version'].split()[:1]==sys.version.split()[:1],'Python release mismatch')
    require(isinstance(runtime.get('julia_version'),str) and runtime['julia_version'].strip(),'missing Julia version')


def paired(joint):
    require(len(joint)==4 and all(type(n) is int and n>=0 for n in joint) and sum(joint)>0,
            'four nonnegative integer joint counts with positive total required')
    neither,a,b,both=joint
    n=sum(joint); delta=(b-a)/n
    def interval(k):
        lo,hi=H6.ENGINE.wilson_interval(k,n,z=2.241402727604947)
        return (0. if k==0 else lo,1. if k==n else hi)
    alo,ahi=interval(a); blo,bhi=interval(b)
    return dict(first='final_only',second='full_prefix',neither=neither,first_only=a,
                second_only=b,both=both,difference=delta,
                standard_error=math.sqrt(max(0,((a+b)/n-delta**2)/(n-1))) if n>1 else None,
                ci_low=max(-1.,blo-ahi),ci_high=min(1.,bhi-alo))


def validate_point(original):
    row=dict(original)
    require(row.get('channel') in CHANNELS,'unknown error channel')
    require(str(row.get('complete')).lower()=='true','incomplete point')
    for key in ('steps','Ly','Lz','shots','target_shots','seed'): row[key]=integer(row[key])
    require(row['steps']==row['Ly']==row['Lz']>=2,'geometry must be T=L')
    require(row['shots']==row['target_shots']>0,'incorrect shot count')
    require(row['seed']>=0,'invalid base seed')
    row['p']=float(row['p']); row['elapsed_seconds']=float(row['elapsed_seconds'])
    require(math.isfinite(row['p']) and 0<=row['p']<.5,'invalid probability')
    require(math.isfinite(row['elapsed_seconds']) and row['elapsed_seconds']>=0,'invalid timing')
    require(re.fullmatch('[0-9a-f]{16}',str(row['point_seed'])) is not None,'invalid exact seed')
    joint=[integer(row[f'joint_{i:02d}']) for i in range(4)]
    require(min(joint)>=0 and sum(joint)==row['shots'],'joint counts do not sum to shots')
    for i,arm in enumerate(ARMS):
        keys=[arm+suffix for suffix in ('_logical_1_failures','_logical_2_failures','_failures')]
        a,b,f=[integer(row[key]) for key in keys]; row.update(zip(keys,(a,b,f)))
        require(0<=max(a,b)<=f<=min(row['shots'],a+b) and min(a,b)>=0,'invalid logical union')
        require(f==sum(n for mask,n in enumerate(joint) if mask&(1<<i)),'joint marginal mismatch')
    row['joint']=joint
    return row


def load_points(directories):
    import importlib.metadata
    sources=[ROOT/'Project.toml',HERE/'run_study.jl',*sorted((ROOT/'src').rglob('*.jl'))]
    if (ROOT/'Manifest.toml').exists(): sources.append(ROOT/'Manifest.toml')
    expected_sources={str(p.relative_to(ROOT)):sha(p) for p in sources}
    packages={name:importlib.metadata.version(name) for name in ('numpy','scipy','pymatching','ldpc')}
    points=[]; inputs=[]; common_runtime=None
    for directory in map(Path,directories):
        require(not (directory/'failed-record.toml').exists(),'failed record in input')
        with (directory/'raw.csv').open(newline='') as stream:
            rows=[validate_point(r) for r in csv.DictReader(stream)]
        checkpoint=tomllib.loads((directory/'checkpoint.toml').read_text())
        metadata=tomllib.loads((directory/'metadata.toml').read_text())
        for key in ('config','runtime','source_hashes','protocol','point_seeds_hex','generation','started_at_unix','updated_at_unix'):
            require(key in checkpoint and checkpoint[key]==metadata.get(key),f'metadata/checkpoint mismatch: {key}')
        require(checkpoint.get('format_version')==metadata.get('format_version')==1,'unsupported format')
        config=checkpoint['config']; channel=config['channel']; runtime=checkpoint['runtime']
        require(checkpoint['protocol']==expected_protocol(channel),'unexpected full-prefix protocol signature')
        require(checkpoint['source_hashes']==expected_sources,'simulation source hashes changed')
        validate_runtime(runtime,packages)
        require(common_runtime is None or common_runtime==runtime,'runtime differs across shards')
        common_runtime=runtime
        sizes,rates=config['sizes'],config['rates']
        require(sizes==sorted(set(sizes)) and all(type(L) is int and L>=2 for L in sizes),'invalid size grid')
        require(rates==sorted(set(rates)) and all(math.isfinite(p) and 0<=p<.5 for p in rates),'invalid rate grid')
        require(config['dimensions']==[[L,L,L] for L in sizes],'configured geometry mismatch')
        grid=[(channel,L,p) for L in sizes for p in rates]
        require([(r['channel'],r['steps'],r['p']) for r in rows]==grid,'raw grid identity/order mismatch')
        require(len(rows)==len(checkpoint['states'])==len(checkpoint['point_seeds_hex'])>0,'grid/state count mismatch')
        require(metadata.get('complete') is True and metadata['expected_points']==metadata['completed_points']==len(rows)
                and metadata['completed_shots']==len(rows)*config['shots'],'incomplete metadata')
        require(metadata['raw_sha256']==sha(directory/'raw.csv') and metadata['checkpoint_sha256']==sha(directory/'checkpoint.toml'),
                'raw/checkpoint hash mismatch')
        blob=(directory/'checkpoint.toml').read_bytes()
        payload,n=re.subn(rb'^payload_sha256 = "[0-9a-f]{64}"\n',b'',blob,flags=re.M)
        require(n==1 and hashlib.sha256(payload).hexdigest()==checkpoint['payload_sha256'],'payload checksum mismatch')
        name=f"checkpoints/checkpoint-{checkpoint['generation']:06d}-{checkpoint['payload_sha256'][:12]}.toml"
        require(metadata['current_snapshot']==name and (directory/name).read_bytes()==blob,'immutable snapshot mismatch')
        require(0<checkpoint['started_at_unix']<=checkpoint['updated_at_unix'],'invalid timestamps')
        for row,state,seed in zip(rows,checkpoint['states'],checkpoint['point_seeds_hex']):
            require(row['seed']==config['seed'] and row['shots']==state['completed_shots']==config['shots'],'shot/base seed mismatch')
            require(row['point_seed']==seed and row['elapsed_seconds']==state['elapsed_seconds'],'point seed/time mismatch')
            require(row['joint']==state['joint_counts'],'checkpoint histogram mismatch')
            for i,arm in enumerate(ARMS):
                require(row[arm+'_failures']==state['failures'][i],'checkpoint marginal mismatch')
                for logical in (1,2):
                    require(row[f'{arm}_logical_{logical}_failures']==state[f'logical_{logical}_counts'][i],'checkpoint logical mismatch')
        points.extend(rows)
        inputs.append(dict(directory=str(directory.resolve()),metadata=metadata,
                           hashes={n:sha(directory/n) for n in ('raw.csv','metadata.toml','checkpoint.toml',name)}))
    require(points,'no inputs')
    keys=[(r['channel'],r['steps'],r['p']) for r in points]
    require(len(keys)==len(set(keys)) and len({r['point_seed'] for r in points})==len(points),'duplicate points/seeds')
    for channel in {r['channel'] for r in points}:
        subset=[r for r in points if r['channel']==channel]
        require(len(subset)==len({r['steps'] for r in subset})*len({r['p'] for r in subset}),'nonrectangular channel grid')
    return sorted(points,key=lambda r:(r['channel'],r['steps'],r['p'])),inputs


def validate_selection(points,inputs,selection,path):
    study=Path(path).resolve().parent
    require(selection.get('stage')=='production' and selection.get('pilot_excluded') is True
            and selection.get('frozen_before_sampling') is True,'frozen production selection required')
    frozen=selection.get('frozen_at_unix')
    require(type(frozen) in (int,float) and math.isfinite(frozen) and frozen>0,'invalid selection freeze time')
    bootstrap=selection.get('bootstrap_replicates')
    require(type(bootstrap) is int and bootstrap>=H6.GUARD_POLICY['minimum_bootstrap_replicates'],'invalid bootstrap count')
    require(set(selection.get('channels',{}))==set(CHANNELS)=={r['channel'] for r in points},'production must declare both measured channels')
    require(selection.get('protocols')=={channel:expected_protocol(channel) for channel in CHANNELS},'frozen protocol signature mismatch')
    shards=selection['shards']
    require(isinstance(shards,dict) and shards and all(isinstance(name,str) and Path(name).name==name
            and name not in ('.','..','pilot-x','pilot-z','pilot-z-refine') for name in shards),'invalid production shard names')
    actual={Path(x['directory']).resolve():x['metadata'] for x in inputs}
    require(len(actual)==len(inputs) and set(actual)=={(study/name).resolve() for name in shards},'input shard set differs from selection')
    for name,declared in shards.items():
        meta=actual[(study/name).resolve()]
        require(all(meta['config'][key]==declared[key] for key in ('channel','sizes','rates','shots','seed')),'frozen shard config mismatch')
        require(frozen<meta['started_at_unix'],'production started before selection freeze')
        require(meta['source_hashes']==selection['source_hashes'],'sources differ from selection')
        require(meta['protocol']==selection['protocols'][declared['channel']],'protocol differs from selection')
    numerical_files=('results/toric3d-h6-z-one-step-2026-10-15/analyze.py',
                     'examples/fit_toric3d_h6_threshold.py','examples/fit_toric3d_diagonal_feedback.py')
    numerical_hashes={name:sha(ROOT/name) for name in numerical_files}
    require(selection.get('numerical_source_sha256')==numerical_hashes,'numerical manifest incomplete or changed')
    pilot_names={'pilot-x':'x_only','pilot-z':'z_only'}
    pilot_hashes=selection.get('pilot_input_sha256',{})
    require(set(pilot_hashes)==set(pilot_names),'pilot manifest must include both channels')
    additional_names={'pilot-z-refine':'z_only'}
    additional_hashes=selection.get('additional_pilot_input_sha256',{})
    require(isinstance(additional_hashes,dict) and set(additional_hashes)==set(additional_names),
            'additional pilot manifest must include the frozen Z refinement')
    all_pilot_names={**pilot_names,**additional_names}
    all_pilot_hashes={**pilot_hashes,**additional_hashes}
    pilot_points=[]; pilot_inputs=[]
    for name in all_pilot_names:
        hashes={filename:sha(study/name/filename) for filename in ('raw.csv','metadata.toml','checkpoint.toml')}
        require(all_pilot_hashes[name]==hashes,'pilot manifest incomplete or evidence changed')
        # Refinement may revisit pilot grid points, but its independent records
        # must each pass the complete raw/checkpoint/protocol/source validator.
        loaded_points,loaded_inputs=load_points([study/name])
        pilot_points.extend(loaded_points); pilot_inputs.extend(loaded_inputs)
    for item in pilot_inputs:
        pilot_meta=item['metadata']; name=Path(item['directory']).name
        require(pilot_meta['config']['channel']==all_pilot_names[name] and pilot_meta['updated_at_unix']<=frozen,
                'pilot channel/completion differs from selection')
        require(pilot_meta['source_hashes']==selection['source_hashes'],'pilot/production source mismatch')
        require(pilot_meta['runtime']==inputs[0]['metadata']['runtime'],'pilot/production runtime mismatch')
    pilot_seeds={r['point_seed'] for r in pilot_points}
    require(len(pilot_seeds)==len(pilot_points),'overlapping pilot point seeds')
    require(not pilot_seeds&{r['point_seed'] for r in points},'pilot/production seed overlap')
    for channel,settings in selection['channels'].items():
        subset=[r for r in points if r['channel']==channel]
        require(sorted({r['steps'] for r in subset})==settings['sizes'] and sorted({r['p'] for r in subset})==settings['rates'],
                'frozen channel grid mismatch')
        require(all(r['shots']==settings['shots'] for r in subset),'frozen sample count mismatch')
        fit_rates=settings.get('fit_rates',[])
        require(fit_rates and all(type(p) in (int,float) and math.isfinite(p) and p>0 for p in fit_rates)
                and fit_rates==sorted(set(fit_rates)) and set(fit_rates)<=set(settings['rates']),
                'fit rates must be an ordered positive subset of sampled rates')
        require(type(settings.get('bootstrap_seed')) is int and settings['bootstrap_seed']>=0,'invalid bootstrap seed')
    require(selection['total_fault_records']==sum(r['shots'] for r in points),'frozen total mismatch')
    return dict(validated=True,selection_sha256=sha(path),numerical_source_sha256=numerical_hashes,
                pilot_input_sha256=pilot_hashes,additional_pilot_input_sha256=additional_hashes,
                bootstrap_replicates=bootstrap,
                bootstrap_seeds={channel:selection['channels'][channel]['bootstrap_seed'] for channel in CHANNELS})


def scaling(points,rates,bootstrap,seed):
    # Only relabel measured counts; the existing numerical model/guards are reused unchanged.
    adapted=[dict(r,one_step_failures=r['full_prefix_failures']) for r in points]
    return OLD.scaling_diagnostics(adapted,rates,bootstrap,seed)


def write_csv(path,rows):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def plot_figures(report,output):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/h6-full-prefix-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from plotting_theme import paper_theme,save_figure,style_axes
    channels=sorted({r['channel'] for r in report['points']})
    def curves(ax,channel,arm):
        subset=[r for r in report['marginals'] if r['channel']==channel and r['arm']==arm]
        for L in sorted({r['L'] for r in subset}):
            rows=[r for r in subset if r['L']==L]; y=np.array([r['rate'] for r in rows])
            ax.errorbar([r['p'] for r in rows],y,
                yerr=np.maximum(0,[y-np.array([r['ci_low'] for r in rows]),np.array([r['ci_high'] for r in rows])-y]),
                color=COLORS.get(L,'#D55E00'),marker='o',label=fr'$L={L}$',elinewidth=.85)
        style_axes(ax); ax.set(xlabel=r'Gate-output fault probability $p$',ylabel=f'Logical {channel[0].upper()} state failure rate',ylim=(0,1))
        ax.margins(x=.04)
    with paper_theme():
        fig,axes=plt.subplots(len(channels),2,figsize=(12.4,4.5*len(channels)),squeeze=False)
        for i,channel in enumerate(channels):
            for j,arm in enumerate(ARMS):
                ax=axes[i,j]; curves(ax,channel,arm)
                ax.set_title(f"{channel[0].upper()} errors · "+('Final-only recovery' if arm=='final_only' else 'Full-prefix complete recovery'))
                ax.legend(loc='upper left',ncol=2,fontsize=9)
        fig.subplots_adjust(left=.08,right=.98,bottom=.1,top=.94,wspace=.23,hspace=.32)
        fig.text(.08,.022,'H+6 · T = L · paired fault records within each channel · 95% Wilson intervals',fontsize=9)
        save_figure(fig,output/'full-prefix-comparison'); plt.close(fig)
        for channel in channels:
            fig,ax=plt.subplots(figsize=(8.2,5.5)); curves(ax,channel,'full_prefix')
            ax.set_title(f'H+6 · Full-prefix complete {channel[0].upper()} recovery')
            ax.legend(loc='upper right',ncol=2,fontsize=9)
            fit=report['fits'].get(channel); candidate=fit['candidate'] if fit else None
            if candidate and not candidate['bound_hits']:
                inset=ax.inset_axes([.12,.69,.38,.23]); xs=[]
                for L in fit['sizes']:
                    rows=[r for r in report['marginals'] if r['channel']==channel and r['arm']=='full_prefix' and r['L']==L and r['p'] in fit['rates']]
                    x=np.array([(r['p']-candidate['pc'])*L**(1/candidate['nu']) for r in rows]); xs.extend(x)
                    y=np.array([r['rate'] for r in rows])
                    inset.errorbar(x,y,yerr=np.maximum(0,[y-np.array([r['ci_low'] for r in rows]),np.array([r['ci_high'] for r in rows])-y]),
                                   color=COLORS.get(L,'#D55E00'),marker='o',ms=2.3,lw=.8,elinewidth=.5,capsize=1)
                grid=np.linspace(min(xs),max(xs),300)
                xnorm=grid/(candidate['p_scale']*candidate['L_reference']**(1/candidate['nu']))
                inset.plot(grid,expit(np.polynomial.polynomial.polyval(xnorm,candidate['coefficients'])),color='0.2',lw=1,
                           ls='-' if fit['status']=='accepted' else '--')
                inset.set(xlabel=r'$(p-p_c)L^{1/\nu}$',ylabel='Logical error')
                inset.xaxis.set_major_locator(MaxNLocator(4)); style_axes(inset,inset=True)
                inset.set_title('Accepted fit' if fit['status']=='accepted' else 'Rejected fit: diagnostic only',fontsize=8)
            annotation='Threshold not established' if not fit or fit['status']!='accepted' else fr'$p_c={fit["pc"]:.5f}$; $\nu={fit["nu"]:.2f}$'
            ax.text(.97,.13,annotation,transform=ax.transAxes,ha='right',fontsize=9)
            fig.subplots_adjust(left=.12,right=.97,bottom=.16,top=.90)
            fig.text(.12,.025,'Full-prefix correction on all available data slices after every transfer, including the last · 95% Wilson bars',fontsize=8.5)
            save_figure(fig,output/f'{channel[0]}-full-prefix-local-fit'); plt.close(fig)
        fig,axes=plt.subplots(1,len(channels),figsize=(7*len(channels),4.7),squeeze=False)
        for ax,channel in zip(axes[0],channels):
            subset=[r for r in report['paired_differences'] if r['channel']==channel]
            for L in sorted({r['L'] for r in subset}):
                rows=[r for r in subset if r['L']==L]; y=np.array([r['difference'] for r in rows])
                ax.errorbar([r['p'] for r in rows],y,yerr=np.maximum(0,[y-np.array([r['ci_low'] for r in rows]),np.array([r['ci_high'] for r in rows])-y]),
                            color=COLORS.get(L,'#D55E00'),marker='o',label=fr'$L={L}$',elinewidth=.85)
            ax.axhline(0,color='0.4',ls='--',lw=.8); style_axes(ax)
            ax.set(title=f'{channel[0].upper()} errors · Full-prefix minus final-only',xlabel=r'Gate-output fault probability $p$',ylabel='Difference in failure probability')
            ax.legend(loc='best',ncol=2,fontsize=8)
        fig.subplots_adjust(left=.07,right=.98,bottom=.19,top=.89,wspace=.25)
        fig.text(.07,.035,'Negative favors full-prefix recovery · approximate pointwise 95% paired intervals (Bonferroni–Wilson)',fontsize=9)
        save_figure(fig,output/'full-prefix-paired-differences'); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',action='append',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--selection',type=Path)
    args=parser.parse_args()
    points,inputs=load_points(args.input)
    selection=json.loads(args.selection.read_text()) if args.selection else None
    evidence=validate_selection(points,inputs,selection,args.selection) if selection else None
    marginals=[]; differences=[]; fits={}
    for row in points:
        for arm in ARMS:
            k,n=row[arm+'_failures'],row['shots']; lo,hi=H6.wilson_interval(k,n)
            marginals.append(dict(channel=row['channel'],L=row['steps'],p=row['p'],arm=arm,shots=n,failures=k,rate=k/n,ci_low=lo,ci_high=hi))
        differences.append(dict(channel=row['channel'],L=row['steps'],p=row['p'],shots=row['shots'],**paired(row['joint'])))
    if selection:
        for channel,settings in selection['channels'].items():
            fits[channel]=scaling([r for r in points if r['channel']==channel],settings['fit_rates'],selection['bootstrap_replicates'],settings['bootstrap_seed'])
    output=args.output; output.mkdir(parents=True,exist_ok=True)
    write_csv(output/'marginals.csv',marginals); write_csv(output/'paired-differences.csv',differences)
    report=dict(arms=list(ARMS),points=points,inputs=inputs,selection=selection,selection_validation=evidence,
                marginals=marginals,paired_differences=differences,fits=fits,
                total_fault_records=sum(r['shots'] for r in points),total_protocol_outcomes=2*sum(r['shots'] for r in points),
                analysis_source_sha256=sha(__file__),
                interpretation='Separate X-only/Z-only protocols; different logical-state seed circuits and final decoders. Every full-prefix correction includes the last transfer. Thresholds are withheld if any frozen guard fails.')
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    plot_figures(report,output)
    print(json.dumps(dict(points=len(points),records=report['total_fault_records'],fits={k:v['status'] for k,v in fits.items()})))


if __name__=='__main__': main()
