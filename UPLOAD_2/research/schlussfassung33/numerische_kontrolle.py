"""Recompute reported errors, outward rounding and the new reading sensitivity.

No import of the original fit implementation. Reuses the published working
data and archived fit parameters for nominal arithmetic; refits only the new
technical contour/axis variant. This is not an independent experiment.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING
import csv,json,platform
import numpy as np
import scipy
from scipy.optimize import least_squares

S=Path(__file__).resolve().parent;R=S.parent
M=json.loads((R/'monografie/manuskript.json').read_text())
TABLES={b['table_number']:b['table'] for c in M for b in c['blocks'] if 'table' in b}
checks=[]
def check(name,condition,detail):
    checks.append(dict(name=name,passed=bool(condition),detail=detail))
    if not condition:raise AssertionError(name+': '+str(detail))
def csvrows(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def err(pred,y):return float(np.sqrt(np.mean((np.asarray(pred)-np.asarray(y))**2)))
def pred(t,p,kind,p0=None):
    if p0 is None:p0,amplitude,k,*shape=p
    else:amplitude,k,*shape=p
    n=1 if kind=='E' else shape[0]
    return p0+amplitude*(-np.expm1(-(k*np.asarray(t))**n))
nominal=[];curves={}
for label,fname in [('A','case_results.json'),('F','F_results.json')]:
    old=json.loads((R/'fallstudie_shintani'/fname).read_text())
    prot=json.loads((R/'fallstudie_shintani'/('digitization_protocol.json' if label=='A' else 'digitization_F.json')).read_text())
    pts=prot['points'];t=np.array([6*(p['pixel_x']-223.5)/722.5 for p in pts])
    y=np.array([p['maltose_g_per_g_input_flour'] for p in pts])
    curves[label]=dict(t=t,y=y,train=np.arange(6),test=np.arange(6,8),unit='g/g',p0=None,old=old['models'])
h=csvrows(R/'datenpruefung_2026_11_05/hackenberg_pilot/digitized_maltose.csv')
old=json.loads((R/'datenpruefung_2026_11_05/hackenberg_pilot/fit_results.json').read_text())
for label in ['H_STD','H_HIGH']:
    pts=[p for p in h if p['curve_id']==label]
    curves[label]=dict(t=np.array([float(p['time_min']) for p in pts]),
        y=np.array([float(p['maltose_mg_per_g_dry_flour']) for p in pts]),
        train=np.arange(1,5),test=np.arange(5,7),unit='mg/g',
        p0=float(pts[0]['maltose_mg_per_g_dry_flour']),old=old['curves'][label]['models'])
for index,(label,c) in enumerate(curves.items()):
    t,y,train,test=c['t'],c['y'],c['train'],c['test']
    table=TABLES[22 if index<2 else 27]
    group=index%2
    persistence=err(np.full(len(test),y[train[-1]]),y[test])
    printed=float(table['rows'][group*3][3].replace(',','.'))
    check(label+' persistence vs manuscript',abs(persistence-printed)<=.0000051,
          dict(calculated=persistence,printed=printed))
    for j,kind in enumerate(['E','W']):
        fit=c['old'][kind]['early'] if index<2 else c['old'][kind]
        pars=fit['parameters'];yp=pred(t,pars,kind,c['p0'])
        early,late=err(yp[train],y[train]),err(yp[test],y[test])
        row=table['rows'][group*3+j+1]
        check(label+' '+kind+' early/late vs manuscript',
            all(abs(a-float(b.replace(',','.')))<=.0000051 for a,b in zip([early,late],row[2:])),
            dict(early=early,late=late,printed=row[2:]))
        nominal.append(dict(curve=label,model=kind,unit=c['unit'],early_RMSE=early,
                            late_RMSE=late,persistence_late_RMSE=persistence))

# Compute exact rational endpoints and verify the displayed decimal intervals.
N,Z,L,dN,dZ,dL=map(F,['.004','.012','.020','.0002','.0002','.001'])
D=Z-N;dD=dN+dZ
intervals={
    'a_hardbounds':((D-dD)/(L+dL),(D+dD)/(L-dL),'.36190','.44211'),
    'v_hardbounds':((L-dL)*(N-dN)/(D+dD),(L+dL)*(N+dN)/(D-dD),'.00859','.01161'),
    'v_kappa':(F('.0078')/F('.82'),F('.0082')/F('.78'),'.009512','.010513'),
    'a_kappa':(F('.0038')*F('.78')/F('.0082'),F('.0042')*F('.82')/F('.0078'),'.36146','.44154')}
rounding=[]
for name,(lo,hi,shown_lo,shown_hi) in intervals.items():
    check(name+' outward rounding',F(shown_lo)<=lo<=hi<=F(shown_hi),
          dict(exact_lower=str(lo),exact_upper=str(hi),display_lower=shown_lo,display_upper=shown_hi))
    rounding.append(dict(name=name,exact_lower=str(lo),exact_upper=str(hi),
        decimal_lower=float(lo),decimal_upper=float(hi),display_lower=shown_lo,display_upper=shown_hi))

# New contour positions plus newly measured common axis anchors. No nominal data
# files are overwritten and the historical global certificates stay nominal.
positions=csvrows(S/'Konturkontrolle_28_Symbole.csv')
axes=csvrows(S/'Achsenkontrolle.csv')
def anchor(study,axis,value):
    return float(next(p['new_pixel'] for p in axes if p['study']==study and p['axis']==axis and float(p['label'])==value))
fresh=[];summaries=[]
for label,c in curves.items():
    study='Shintani' if label in ['A','F'] else 'Hackenberg'
    pts=[p for p in positions if p['point_id'].startswith(label) and p['study']==study]
    if study=='Shintani':
        x0,x6=anchor(study,'x',0),anchor(study,'x',6)
        y0,yt=anchor(study,'y',0),anchor(study,'y',.3)
        t=np.array([6*(float(p['new_x'])-x0)/(x6-x0) for p in pts])
        y=np.array([.3*(y0-float(p['new_y']))/(y0-yt) for p in pts])
    else:
        y0,yt=anchor(study,'y',0),anchor(study,'y',50)
        t=c['t'];y=np.array([c['p0']]+[50*(y0-float(p['new_y']))/(y0-yt) for p in pts])
    tr,te=c['train'],c['test'];p0=c['p0']
    persist=err(np.full(len(te),y[tr[-1]]),y[te]);late={}
    for kind in ['E','W']:
        old=c['old'][kind]['early'] if study=='Shintani' else c['old'][kind]
        initial=np.array(old['parameters']);lb=np.array(old['bounds_lower']);ub=np.array(old['bounds_upper'])
        if study=='Shintani':ub[0]=min(y[tr])
        starts=[initial]
        kidx=2 if p0 is None else 1
        for factor in [.5,2]:
            z=initial.copy();z[kidx]*=factor;starts.append(z)
        fits=[]
        for start in starts:
            start=np.maximum(lb+1e-12,np.minimum(ub-1e-12,start))
            fit=least_squares(lambda p:pred(t[tr],p,kind,p0)-y[tr],start,bounds=(lb,ub),
                xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=20000)
            fits.append(fit)
        fit=min(fits,key=lambda z:float(np.dot(z.fun,z.fun)))
        check(label+' '+kind+' new contour refit converged',fit.success,fit.message)
        late[kind]=err(pred(t[te],fit.x,kind,p0),y[te])
        fresh.append(dict(curve=label,model=kind,unit=c['unit'],parameters=fit.x.tolist(),
            early_RMSE=err(pred(t[tr],fit.x,kind,p0),y[tr]),late_RMSE=late[kind],
            persistence_late_RMSE=persist,starts=3,
            status='technical sensitivity variant; no global search certificate'))
    nominal_pair=[v for v in nominal if v['curve']==label]
    old_sign=np.sign(nominal_pair[1]['late_RMSE']-nominal_pair[0]['late_RMSE'])
    check(label+' contour variant retains E/W order',np.sign(late['W']-late['E'])==old_sign,late)
    summary=dict(curve=label,unit=c['unit'],max_abs_value_difference=float(max(abs(y-c['y']))),
        max_abs_symbol_dy=max(abs(float(p['dy'])) for p in pts),
        persistence_RMSE=persist,E_late_RMSE=late['E'],W_late_RMSE=late['W'])
    summaries.append(summary)
    if study=='Shintani':
        check(label+' persistence still better in contour variant',persist<min(late.values()),summary)
        check(label+' late decrease remains in contour variant',y[5]>y[6],dict(fall=float(y[5]-y[6])))
    else:check(label+' E still better than persistence',late['E']<persist,summary)
check('28 symbols retained inside original local scenarios',all(
    abs(float(p['dy']))<=float(p['scenario_halfwidth_px']) and
    abs(float(p['dx']))<=float(p['scenario_halfwidth_px']) for p in positions),len(positions))
result=dict(status='PASS',version=33,checks=checks,nominal_errors=nominal,
    rounding=rounding,contour_refits=fresh,contour_summary=summaries,
    note='Same AI-assisted project; no human independent rereading or new experiment.',
    versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
(S/'Numerische_Kontrolle.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status='PASS',checks=len(checks),contour_summary=summaries),ensure_ascii=False,indent=2))
