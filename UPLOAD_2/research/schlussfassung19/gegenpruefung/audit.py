from pathlib import Path
import csv,json,hashlib
import mpmath as mp
import numpy as np
from scipy.optimize import linprog
import sympy as s
R=Path(__file__).resolve().parent
mp.mp.dps=100
rows=list(csv.DictReader((R/'nominal_data.csv').open()))
candidates=json.loads((R/'candidates.json').read_text())['candidates']
cert={x['curve']:x for x in csv.DictReader((R/'certificates.csv').open())}
results=[]
for c in candidates:
 rr=[x for x in rows if x['curve']==c['curve']]
 t=[mp.mpf(x['time']) for x in rr];y=[mp.mpf(x['value']) for x in rr]
 p,a,d=[mp.mpf(c[k]) for k in ['p0','amplitude','rate']]
 residual=[yi-p-a*(-mp.expm1(-d*ti)) for ti,yi in zip(t,y)]
 eps=max(abs(x) for x in residual)
 assert mp.mpf(cert[c['curve']]['lower'])<=eps<=mp.mpf(cert[c['curve']]['upper'])
 ids=c['suggested_active_indices_zero_based'];sign=np.array([1 if residual[i]>0 else -1 for i in ids])
 assert all(sign[j]!=sign[j+1] for j in range(3))
 J=np.array([[1,float(-mp.expm1(-d*t[i])),float(a*t[i]*mp.exp(-d*t[i]))] for i in ids])
 mat=np.vstack([(sign[:,None]*J).T,np.ones(4)])
 w=np.linalg.solve(mat,np.array([0,0,0,1.]));assert min(w)>0
 tf=np.array([float(x) for x in t]);yf=np.array([float(x) for x in y]);lp=[];failed=0
 for rate in float(d)*np.logspace(-4,4,81):
  u=-np.expm1(-rate*tf)
  M=np.vstack([np.column_stack([-np.ones(len(t)),-u,-np.ones(len(t))]),np.column_stack([np.ones(len(t)),u,-np.ones(len(t))])])
  fit=linprog([0,0,1],A_ub=M,b_ub=np.r_[-yf,yf],bounds=[(0,None)]*3,method='highs')
  if fit.success:
   lp.append(fit.fun)
   assert fit.fun>=float(cert[c['curve']]['lower'])-1e-7*max(1,float(eps))
  else:failed+=1
 results.append({'curve':c['curve'],'max_residual_100_digits':mp.nstr(eps,70),'dual_weights':w.tolist(),'stationarity_residual':float(np.max(np.abs(mat@w-[0,0,0,1]))),'fixed_rate_LP_success':len(lp),'fixed_rate_LP_failures':failed,'best_grid_error':min(lp),'scope':'High precision and stationary weights are crosschecks; only the original interval/alternation proof gives the global certificate.'})
# Exact algebraic reconstruction, gradient and origin criterion.
A,B,Y,h,t=s.symbols('A B Y h t',positive=True)
d=s.log(A/B)/h
amp=A**2/(A-B)*(A/B)**(t/h)
p0=Y-amp*(1-s.exp(-d*t))
expected=Y-A**2/(A-B)*((A/B)**(t/h)-1)
assert s.simplify(s.expand_power_base(p0-expected,force=True))==0
x1,x2,x3=s.symbols('x1 x2 x3');expr=s.log((x2-x1)/(x3-x2))/h
g=[s.diff(expr,x) for x in [x1,x2,x3]]
assert s.simplify(sum(g))==0
subs={x1:0,x2:A,x3:A+B}
var=s.factor(sum(z.subs(subs)**2 for z in g))
assert s.simplify(var-2/h**2*(1/A**2+1/B**2+1/(A*B)))==0
# Recompute every leave-one-labor prediction from its own training-only ratios.
assay=list(csv.DictReader((R/'aktivitaeten.csv').open(encoding='utf-8-sig')))
errors={k:[] for k in ['identity','common','product']}
for lab in sorted({x['lab'] for x in assay}):
 train=[x for x in assay if x['lab']!=lab];test=[x for x in assay if x['lab']==lab]
 q=mp.exp(mp.fsum(mp.log(mp.mpf(x['activity_37'])/mp.mpf(x['activity_20'])) for x in train)/len(train))
 for x in test:
  tt=[z for z in train if z['product']==x['product']]
  qp=mp.exp(mp.fsum(mp.log(mp.mpf(z['activity_37'])/mp.mpf(z['activity_20'])) for z in tt)/len(tt))
  for name,factor in [('identity',1),('common',q),('product',qp)]:
   rat=mp.mpf(x['activity_20'])*factor/mp.mpf(x['activity_37']);errors[name].append((mp.log(rat)**2,abs(rat-1)))
metrics={k:{'log_RMSE':float(mp.sqrt(mp.fsum(x[0] for x in v)/len(v))),'MAPE':float(mp.fsum(x[1] for x in v)/len(v))} for k,v in errors.items()}
assert abs(metrics['common']['log_RMSE']-.1646696659111359)<1e-12
assert abs(metrics['product']['log_RMSE']-.1644684553487647)<1e-12
# Compare all 80 numeric activity/SD fields with archived original JATS table.
import xml.etree.ElementTree as ET
from decimal import Decimal
original={}
for tr in ET.parse(R/'quelle_tabelle1.xml').findall('.//tbody/tr'):
 cells=[''.join(c.itertext()).strip() for c in tr]
 if len(cells)==33 and cells[0].startswith('Lab '):
  lab=cells[0][-1]
  if lab not in original: original[lab]=cells
assert set(original)==set('ABCDE')
products=['amylase_M','amylase_S','pancreatin','saliva']
checked=0
for row in assay:
 j=products.index(row['product']); cells=original[row['lab']]
 for field,offset in [('activity_20',1),('published_SD_20',3),('activity_37',5),('published_SD_37',7)]:
  assert Decimal(row[field])==Decimal(cells[8*j+offset])
  checked+=1
assert checked==80
out={'original_table_numeric_fields_matched':checked,'status' :'PASS','precision_decimal_digits':100,'curves':results,'symbolic_checks':['nonnegative_origin_identity','common_offset_gradient_zero','three_point_variance_with_shared_middle_observation'],'assay':metrics,'human_independent_review':False,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in R.iterdir() if p.is_file() and p.name not in ['Ergebnisse.json']}}
(R/'Ergebnisse.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
