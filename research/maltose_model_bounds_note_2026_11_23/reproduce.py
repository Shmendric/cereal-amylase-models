"""Reproduce retrospective maltose-model bound crossovers (2026).

No new laboratory data or statistical confidence intervals. The two time series
are digitised values already published and analysed in the 2026 monograph.
Run from this directory in the repository with: python reproduce.py
"""
from __future__ import annotations
import csv, hashlib, json, platform, sys
from functools import lru_cache
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize, brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'zertifizierte_fehlergrenzen/data/original/Hackenberg.csv'
OUT = HERE / 'generated'
OUT.mkdir(exist_ok=True)
EXPECTED_SHA256 = 'd40867e07ec1333931a46df1cbe88a0eef42833782d1714f9dab0f7e5e'
if hashlib.sha256(SOURCE.read_bytes()).hexdigest() != EXPECTED_SHA256:
    raise ValueError('Source data SHA256 changed; this analysis uses the archived technical readings')
with SOURCE.open(encoding='utf-8-sig', newline='') as f:
    source_rows=list(csv.DictReader(f))
DATA={}
for name in ('H_STD','H_HIGH'):
    rr=[r for r in source_rows if r['curve_id']==name]
    t=np.array([float(r['time_min']) for r in rr]); y=np.array([float(r['maltose_mg_per_g_dry_flour']) for r in rr])
    if len(t)!=7 or not np.all(np.diff(t)>0): raise ValueError('Unexpected time series')
    DATA[name]=(t,y)

def prediction(t,p0,delta,k,n):
    return p0+delta*(-np.expm1(-np.power(k*t,n)))

@lru_cache(maxsize=4096)
def fit(name,scenario,model,delta_max,n_max):
    t,y=DATA[name]; yy=y.copy()
    if scenario=='rounded': yy=np.round(yy,2)
    if scenario=='late_plus': yy[5:]+=.236
    if scenario=='late_minus': yy[5:]-=.236
    p0=float(yy[0]); tt=t[1:5]; obs=yy[1:5]; test=yy[5:]
    def objective(z):
        k=float(np.exp(z[0])); n=float(z[1]) if model=='W' else 1.
        h=-np.expm1(-np.power(k*tt,n)); denom=float(h@h)
        delta=float(np.clip(((obs-p0)@h)/denom if denom else 0.,0.,delta_max))
        residual=p0+delta*h-obs
        return float(residual@residual),delta,k,n
    bounds=[(np.log(1e-6),np.log(10.))]+([(0.1,n_max)] if model=='W' else [])
    starts=[(np.log(k),) for k in np.logspace(-5,-1,7)] if model=='E' else [
        (np.log(k),n) for k in np.logspace(-4,-1,6) for n in np.linspace(.1,n_max,5)]
    candidates=[]
    for start in starts:
        opt=minimize(lambda x:objective(x)[0],np.asarray(start),bounds=bounds,method='L-BFGS-B',options={'maxiter':500,'ftol':1e-14,'gtol':1e-12})
        if np.isfinite(opt.fun): candidates.append(opt)
    if not candidates: raise RuntimeError('No finite optimisation result')
    best=min(candidates,key=lambda z:z.fun)
    loss,delta,k,n=objective(best.x)
    later=float(np.sqrt(np.mean((prediction(t[5:],p0,delta,k,n)-test)**2)))
    return dict(delta=delta,k=k,n=n,early_RMSE=float(np.sqrt(loss/4)),late_RMSE=later,
                delta_max=delta_max,n_max=n_max,at_delta_bound=bool(abs(delta-delta_max)<1e-4),
                at_shape_bound=bool(model=='W' and abs(n-n_max)<1e-4))

def diff_std(scenario,n_max):
    return fit('H_STD',scenario,'W',100.,float(n_max))['late_RMSE']-fit('H_STD',scenario,'E',100.,4.)['late_RMSE']

def diff_high(scenario,delta_max):
    d=float(delta_max)
    return fit('H_HIGH',scenario,'W',d,4.)['late_RMSE']-fit('H_HIGH',scenario,'E',d,4.)['late_RMSE']

def scan_roots(function,grid,xtol):
    vals=[function(float(x)) for x in grid]
    return [float(brentq(function,float(a),float(b),xtol=xtol)) for a,b,u,v in
            zip(grid[:-1],grid[1:],vals[:-1],vals[1:]) if u*v<0]

SCENARIOS=() if '--figures-only' in sys.argv else ('nominal','rounded','late_plus','late_minus')
results=[]
for scenario in SCENARIOS:
    roots_std=scan_roots(lambda x:diff_std(scenario,x),[1.85,1.9,1.94,1.98,2.,2.02,2.06,2.10,2.16],5e-5)
    roots_high=scan_roots(lambda x:diff_high(scenario,x),[35.,40.,45.,50.,55.,60.,70.,80.,100.],2e-3)
    results.append(dict(scenario=scenario,nmax_crossovers_MSM_4_78=roots_std,
                        delta_max_crossovers_MSM_8_37_mg_g=roots_high))
    print(scenario,results[-1],flush=True)
if not SCENARIOS: print('Only regenerating plots; scenario outputs remain unchanged',flush=True)
if SCENARIOS:
    with (OUT/'crossover_scenarios.json').open('w',encoding='utf-8') as f:
        json.dump(dict(source_sha256=EXPECTED_SHA256, results=results,
            interpretation='Deterministic perturbations of technical digitisation; no confidence intervals',
            software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)),f,ensure_ascii=False,indent=2)
    with (OUT/'crossover_scenarios.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f);writer.writerow(['scenario','nmax_lower','nmax_upper','delta_max_mg_per_g'])
        for r in results:
            ns=r['nmax_crossovers_MSM_4_78'];ds=r['delta_max_crossovers_MSM_8_37_mg_g']
            writer.writerow([r['scenario'],ns[0] if len(ns)>0 else '',ns[1] if len(ns)>1 else '',ds[0] if ds else ''])

# Plots corresponding to three figures in the English manuscript.
plt.rcParams.update({'font.size':9})
fig,axes=plt.subplots(1,2,figsize=(9.2,3.3),constrained_layout=True)
for ax,name,xs,fun in [(axes[0],'MSM 4.78',np.linspace(1.76,2.5,22),lambda x:diff_std('nominal',x)),
                      (axes[1],'MSM 8.37',np.linspace(35,100,22),lambda x:diff_high('nominal',x))]:
    ax.plot(xs,[fun(float(x)) for x in xs]);ax.axhline(0,color='0.4',lw=.8,ls='--')
    ax.set(title=name, ylabel='Late RMSE: W − E (mg/g)',xlabel=('Upper shape bound nmax' if name.endswith('4.78') else 'Upper amplitude bound (mg/g)'))
    ax.grid(alpha=.15)
fig.savefig(OUT/'figure_1_crossover_curves.svg');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(9.2,3.5),constrained_layout=True)
for ax,name,variants in [(axes[0],'H_STD',[('E','E',100.,4.),('W nmax 1.97','W',100.,1.97),('W nmax 4','W',100.,4.)]),
                         (axes[1],'H_HIGH',[('E','E',100.,4.),('W cap 50','W',50.,4.),('W cap 100','W',100.,4.)])]:
    t,y=DATA[name]; ax.scatter(t[1:5],y[1:5],marker='o',s=18,label='Fit data');ax.scatter(t[5:],y[5:],marker='s',s=30,label='Later points');ax.scatter(t[0],y[0],marker='x',s=30,label='Fixed initial value')
    grid=np.linspace(0,180,350)
    for label,model,cap,nmax in variants:
        v=fit(name,'nominal',model,cap,nmax);ax.plot(grid,prediction(grid,y[0],v['delta'],v['k'],v['n']),label=label)
    ax.axvline(60,color='0.4',ls=':',lw=.8)
    ax.set(title=name,xlabel='Time (min)',ylabel='Maltose (mg/g dry flour)');ax.legend(fontsize=7,frameon=False)
fig.savefig(OUT/'figure_2_model_continuations.svg');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(9.2,3.3),constrained_layout=True)
for ax,name,xx in [(axes[0],'H_STD',[1.915,1.97,2.,2.04,2.2,4.]),(axes[1],'H_HIGH',[35.,50.,53.45,60.,75.,100.])]:
    fits=[fit(name,'nominal','W',100. if name=='H_STD' else x,x if name=='H_STD' else 4.) for x in xx]
    ax.plot(xx,[z['early_RMSE'] for z in fits],'o-',label='Early RMSE');ax.plot(xx,[z['late_RMSE'] for z in fits],'s-',label='Late RMSE')
    ax.set(title=name,ylabel='RMSE (mg/g)',xlabel=('Upper shape bound nmax' if name=='H_STD' else 'Upper amplitude bound (mg/g)'));ax.legend(frameon=False)
fig.savefig(OUT/'figure_3_early_late_rmse.svg');plt.close(fig)
print('Finished. Files written to',OUT,flush=True)
