"""Recheck the archived version 14 data for the corrected manuscript version 16.

No original project functions are imported. Data and archived results are read
from the supplied package. These are technical checks and counterexamples, not
new laboratory observations or a proof of scientific priority.
"""
from pathlib import Path
from decimal import Decimal
import argparse, csv, hashlib, json, math, platform
import numpy as np
import scipy
from scipy.optimize import differential_evolution, lsq_linear, minimize_scalar
from scipy.special import wrightomega
from scipy.integrate import solve_ivp

ap=argparse.ArgumentParser()
ap.add_argument('--package',type=Path,required=True)
ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parent)
a=ap.parse_args();R=a.package.resolve();O=a.output.resolve();O.mkdir(parents=True,exist_ok=True)
checks={}
def readcsv(p):
    with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def savecsv(name,rows):
    with (O/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def rmse(x):return float(np.sqrt(np.mean(np.asarray(x)**2)))
def response(t,p0,b,d):return p0+b*(-np.expm1(-d*np.asarray(t)))/d

# Linear parameters P0 and delta are profiled under the original limits.
# The original fits optimised all model parameters jointly.
def fit_profile(t,y,kind):
    p0cap=float(y.min())
    def one(z,details=False):
        k=float(np.exp(z[0]));n=1. if kind=='E' else float(z[1])
        f=-np.expm1(-(k*t)**n)
        X=np.column_stack((np.ones(len(t)),f))
        lin=lsq_linear(X,y,bounds=([0.,1e-5],[p0cap,.8]),method='bvls',tol=1e-13)
        err=X@lin.x-y
        if details:return {'P0':float(lin.x[0]),'delta':float(lin.x[1]),'k_per_h':k,'n':n,'SSE':float(err@err)}
        return float(err@err)
    bounds=[(math.log(.001),math.log(30.))]
    if kind=='E':
        grid=np.linspace(*bounds[0],151);vals=np.array([one([g]) for g in grid]);i=int(vals.argmin())
        z=minimize_scalar(lambda z:one([z]),bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),method='bounded',options={'xatol':1e-13})
        best=[float(z.x)]
    else:
        sol=differential_evolution(one,bounds+[(.1,3.)],seed=2026110515,popsize=16,maxiter=450,tol=1e-11,polish=True)
        best=sol.x
    return one(best,True)
def pred(p,t):return p['P0']+p['delta']*(-np.expm1(-(p['k_per_h']*t)**p['n']))
def isotonic(y):
    # Pool adjacent violators, equal weights. Gives a lower SSE bound for every
    # increasing curve on these discrete time points, including E and W.
    blocks=[]
    for v in y:
        blocks.append([float(v),1])
        while len(blocks)>1 and blocks[-2][0]>blocks[-1][0]:
            u,n=blocks.pop();v,m=blocks.pop();blocks.append([(u*n+v*m)/(n+m),n+m])
    return np.concatenate([np.repeat(v,n) for v,n in blocks])

forecasts=[];fitout=[];isotonic_out={}
for sample,tf in [('A','time_digitized_h'),('F','time_h')]:
    rows=readcsv(R/f'fallstudie_shintani/digitized_{sample}.csv')
    t=np.array([float(r[tf]) for r in rows]);y=np.array([float(r['maltose_g_per_g_input_flour']) for r in rows])
    archived=json.loads((R/('fallstudie_shintani/case_results.json' if sample=='A' else 'fallstudie_shintani/F_results.json')).read_text())
    naive=np.repeat(y[5],2)
    forecasts.append({'source':'Shintani','material':sample,'model':'letzter_frueher_Wert','unit':'g/g input flour','late_RMSE':rmse(naive-y[6:])})
    iso=isotonic(y);isotonic_out[sample]={'all_point_monotone_RMSE_floor':rmse(iso-y),'nondecreasing_projection':iso.tolist()}
    for kind in ['E','W']:
        q=fit_profile(t[:6],y[:6],kind);predicted=pred(q,t[6:]);late=rmse(predicted-y[6:])
        am=archived['models'][kind]['early']
        archived_late=am['late_RMSE'] if sample=='A' else archived['models'][kind]['late_RMSE']
        checks[f'{sample}_{kind}_training_objective_reproduced']=abs(q['SSE']-am['SSE'])<1e-9
        checks[f'{sample}_{kind}_late_error_reproduced']=abs(late-archived_late)<1e-5
        forecasts.append({'source':'Shintani','material':sample,'model':kind,'unit':'g/g input flour','late_RMSE':late})
        fitout.append({'sample':sample,'model':kind,**q,'late_RMSE':late,'late_predictions':predicted.tolist()})

H=R/'datenpruefung_2026_11_05/hackenberg_pilot'
hr=readcsv(H/'digitized_maltose.csv');hf=json.loads((H/'fit_results.json').read_text())['curves']
for cid in ['H_STD','H_HIGH']:
    rr=[r for r in hr if r['curve_id']==cid];v={float(r['time_min']):float(r['maltose_mg_per_g_dry_flour']) for r in rr}
    ty=np.array([v[120.],v[180.]])
    forecasts.append({'source':'Hackenberg','material':cid,'model':'letzter_frueher_Wert','unit':'mg/g dry flour','late_RMSE':rmse(np.repeat(v[60.],2)-ty)})
    for kind in ['E','W']:
        p=hf[cid]['models'][kind]['parameters'];p0=v[0.];shape=1. if kind=='E' else p[2]
        yhat=p0+p[0]*(-np.expm1(-(p[1]*np.array([120.,180.]))**shape))
        forecasts.append({'source':'Hackenberg','material':cid,'model':kind,'unit':'mg/g dry flour','late_RMSE':rmse(yhat-ty)})
savecsv('vorhersage_vergleich.csv',forecasts)
scenario_benchmarks={}
for sample,name in [('A','pixel_scenarios.npz'),('F','F_pixel_scenarios.npz')]:
    with np.load(R/'fallstudie_shintani'/name,allow_pickle=False) as ar:
        yy=ar['input_y'];naive=np.sqrt(np.mean((yy[:,5,None]-yy[:,6:])**2,axis=1))
        scenario_benchmarks[sample]={'n':len(naive),'persistence_better_than_E_count':int(np.sum(naive<ar['late_E'])),'persistence_better_than_W_count':int(np.sum(naive<ar['late_W'])),'scope':'Existing bounded graphical scenarios; archived E/W fits reused, not new experimental replicates or probabilities.'}

# Assay reanalysis through independent linear-design-matrix fits in log space.
ar=readcsv(R/'externe_pruefung/aktivitaeten.csv');products=['amylase_M','amylase_S','pancreatin','saliva']
v20=np.array([[float(next(r for r in ar if r['lab']==lab and r['product']==p)['activity_20']) for p in products] for lab in 'ABCDE'])
v37=np.array([[float(next(r for r in ar if r['lab']==lab and r['product']==p)['activity_37']) for p in products] for lab in 'ABCDE'])
logratio=np.log(v37/v20);assay_rows=[];errors={m:[] for m in ['identity','common','product']}
for hold,lab in enumerate('ABCDE'):
    train=np.arange(5)!=hold;target=logratio[train].ravel()
    z1=np.linalg.lstsq(np.ones((16,1)),target,rcond=None)[0]
    z4=np.linalg.lstsq(np.tile(np.eye(4),(4,1)),target,rcond=None)[0]
    for model,z in [('identity',np.zeros(4)),('common',np.repeat(z1[0],4)),('product',z4)]:
        e=z-logratio[hold];errors[model].append(e)
        assay_rows.append({'held_lab':lab,'model':model,'log_RMSE':rmse(e),'absolute_relative_error':float(np.abs(np.expm1(e)).mean())})
assay={m:{'log_RMSE':rmse(e),'MAPE':float(np.abs(np.expm1(e)).mean())} for m,e in errors.items()}
archived=json.loads((R/'externe_pruefung/ergebnisse.json').read_text())
checks['assay_results_reproduced']=all(abs(assay[m]['log_RMSE']-archived['prediction_errors'][m]['log_RMSE'])<1e-12 and abs(assay[m]['MAPE']-archived['prediction_errors'][m]['mean_absolute_relative_error'])<1e-12 for m in assay)
savecsv('assay_laborfehler.csv',assay_rows)

# Structural reconstruction, variance and finite-domain restrictions.
P0,b,d,t1,h=.015,.004,.04,5.,10.
t=t1+h*np.arange(3);y=response(t,P0,b,d);A,B=np.diff(y)
dr=float(np.log(A/B)/h);br=float(A*dr*np.exp(dr*t1)/(-np.expm1(-dr*h)));p0r=float(y[0]-br*(-np.expm1(-dr*t1))/dr)
checks['three_point_reconstruction']=np.allclose([p0r,br,dr],[P0,b,d],rtol=1e-12,atol=1e-12)
g=np.array([-1/A,1/A+1/B,-1/B])/h
def d_from_y(Y):return np.log((Y[1]-Y[0])/(Y[2]-Y[1]))/h
gd=[]
for j in range(3):
    yy=y.astype(complex);yy[j]+=1e-25j;gd.append(float(d_from_y(yy).imag/1e-25))
checks['three_point_gradient']=np.allclose(g,gd,rtol=1e-11)
precision=[]
def rel_sd(step):
    a1,b1=np.diff(response(t1+step*np.arange(3),P0,b,d))
    grad=np.array([-1/a1,1/a1+1/b1,-1/b1])/step
    return float(.001*np.linalg.norm(grad)/d)
opt=minimize_scalar(rel_sd,bounds=(1.,100.),method='bounded',options={'xatol':1e-10})
for step in [5.,10.,20.,opt.x,60.,90.]:precision.append({'h_min':float(step),'relative_SD_d':rel_sd(step)})
checks['local_design_optimum_reproduced']=abs(d*opt.x-1.590407)<1e-5
crit=float(min(B/2,(A-B)/4));eps=.45*crit
corners=np.array([d_from_y(y+eps*np.array([a1,a2,a3])) for a1 in [-1,1] for a2 in [-1,1] for a3 in [-1,1]])
lo=math.log((A-2*eps)/(B+2*eps))/h;hi=math.log((A+2*eps)/(B-2*eps))/h
checks['sharp_box_bounds_reproduced']=abs(corners.min()-lo)<1e-12 and abs(corners.max()-hi)<1e-12
L=np.array([[0,1,-1],[1,-1,1]])
cov=L@np.diag(np.square([.05,.03,.10]))@L.T
checks['calibration_product_variance_preserved']=abs(np.ones(2)@cov@np.ones(2)-.05**2)<1e-14

# The substrate-pool model and the minimal model have equal initial slopes here,
# but they are not interchangeable over an arbitrary duration.
S0,c,rho,dec=.280,1.1,1.5,.12
b_initial=c*S0/(1+rho)
tt=np.array([0.,.5,1.,2.,10.]);G=-np.expm1(-dec*tt)/dec
u=wrightomega(np.log(rho)+rho-c*G)/rho
Ppool=P0+S0*(1-u);Pminimal=response(tt,P0,b_initial,dec)
pool_ode=solve_ivp(lambda t,y:[-c*np.exp(-dec*t)*y[0]/(1+rho*y[0])],[0,10],[1.],t_eval=tt,method='DOP853',rtol=1e-12,atol=1e-13)
checks['pool_closed_and_ODE_agree']=np.max(abs(pool_ode.y[0]-u))<1e-9
checks['pool_bound_reproduced']=bool(np.all((Ppool>=P0)&(Ppool<=P0+S0)))
missing_condition={'P0_counterexample':-.01,'S0':S0,'K':S0/rho,'c':c,'d':dec,'P_at_zero':-.01,'message':'The counterexample with P0<0 illustrates why the corrected manuscript explicitly requires P0>=0 for a nonnegative product. It is not a counterexample within that corrected domain.'}
mass_limit={'b_from_pool_initial_slope':b_initial,'d':dec,'minimal_asymptotic_increment':b_initial/dec,'available_pool':S0,'finite_pool_condition_for_all_time':'b/d <= available pool for d>0; for d=0,b>0 a closed finite-pool system has only finite-time validity.','time_h':tt.tolist(),'minimal_P':Pminimal.tolist(),'pool_P':Ppool.tolist(),'message':'Counterexample to treating the two model classes as interchangeable, not to their own conditional reconstruction formulas.'}
checks['minimal_finite_pool_counterexample']=b_initial/dec>S0
savecsv('modell_gueltigkeit.csv',[{'time_h':float(ti),'minimal_P':float(pm),'pool_P':float(pp)} for ti,pm,pp in zip(tt,Pminimal,Ppool)])

# Local information in an observation and a deterministic transform of it.
tx=np.array([5.,10.,20.,60.]);plateau=100.;k=.02;sig=.5
mu=plateau*(-np.expm1(-k*tx));J=(plateau*tx*np.exp(-k*tx))[:,None]
D=1/(tx*(plateau-mu));C=sig**2*np.eye(4);Bmat=np.vstack([np.eye(4),np.diag(D)])
joint=Bmat@C@Bmat.T;JJ=Bmat@J
I0=float((J.T@np.linalg.solve(C,J)).item());I1=float((JJ.T@np.linalg.pinv(joint)@JJ).item())
naive=np.diag(np.diag(joint));Iwrong=float((JJ.T@np.linalg.solve(naive,JJ)).item())
checks['joint_information_preserved']=abs(I1/I0-1)<1e-10
checks['false_independence_doubles_local_information']=abs(Iwrong/I0-2)<1e-10

# Exact scalar ACS boundary; scope remains conditional on transcription.
gap=(Decimal(1)-Decimal('-.322011'))*Decimal('161.4')-Decimal('213.49')
checks['ACS_exact_negative_source_reproduced']=str(gap)=='-0.1174246'

checks['units_retained_for_each_source']=all(r['unit']==('g/g input flour' if r['source']=='Shintani' else 'mg/g dry flour') for r in forecasts)
out={'source_data_version':14,'target_manuscript_version':17,'calculation_status':'fresh_implementations_without_importing_project_functions','runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'checks':{k:bool(v) for k,v in checks.items()},'all_checks_pass':bool(all(checks.values())),'nominal_forecasts':forecasts,'profiled_Shintani_fits':fitout,'monotone_lower_bounds':isotonic_out,'assay':assay,'reconstruction':{'P0':p0r,'b':br,'d':dr,'gradient':g.tolist(),'critical_absolute_error':crit,'optimal_h':float(opt.x),'precision':precision},'calibration_log_covariance':cov.tolist(),'minimal_finite_pool':mass_limit,'pool_initial_condition':missing_condition,'dependent_information':{'covariance_rank':int(np.linalg.matrix_rank(joint)),'correct_information_ratio':I1/I0,'naive_information_ratio':Iwrong/I0},'ACS_source_term':str(gap),'limitations':['All forecast comparisons are retrospective and use published processed values.','No missing experimental variances are inferred.','The counterexamples delimit model domains; they are not measured cereal observations.','These checks are internal technical checks, not an independent expert opinion or a proof of novelty.']}
out['existing_pixel_scenarios_with_persistence']=scenario_benchmarks
(O/'rechenpruefung.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
assert all(checks.values()),checks
print(json.dumps({'checks_passed':len(checks),'forecasts':forecasts,'assay':assay,'minimal_pool_limit':mass_limit['minimal_asymptotic_increment'],'finite_pool':S0},ensure_ascii=False,indent=2))
