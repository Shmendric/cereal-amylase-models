"""Fresh symbolic and numerical checks of the manuscript's central relations.

No project model functions are imported. Supplied data and source transcriptions
remain inputs and are not made experimentally independent by this calculation.
"""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
import csv,json,math
import numpy as np
import sympy as s
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.special import wrightomega

S=Path(__file__).resolve().parent;R=S.parent;checks={};details={}
def check(name,value):
 checks[name]=bool(value)
 if not value:raise AssertionError(name)
read=lambda p:json.loads(p.read_text())
rows=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig')))
m=read(R/'monografie/manuskript.json');tables={b['table_number']:b['table'] for c in m for b in c['blocks'] if 'table'in b}

t,b,d,p0,h,t1=s.symbols('t b d p0 h t1',positive=True)
P=p0+b*(1-s.exp(-d*t))/d
check('Eq22 first and second initial derivatives',s.simplify(s.diff(P,t).subs(t,0)-b)==0 and s.simplify(-s.diff(P,t,2).subs(t,0)/b-d)==0)
check('Eq23 two-time quotient',s.simplify((P.subs(t,2*t1)-p0)/(P.subs(t,t1)-p0)-(1+s.exp(-d*t1)))==0)
check('Minimal model d=0 continuous limit',s.simplify(s.limit(P,d,0)-(p0+b*t))==0)
AA=s.simplify(P.subs(t,t1+h)-P.subs(t,t1));BB=s.simplify(P.subs(t,t1+2*h)-P.subs(t,t1+h))
check('Eq27 successive increment ratio',s.simplify(BB/AA-s.exp(-d*h))==0)
check('Eq29 Aitken plateau',s.simplify(P.subs(t,t1)+AA**2/(AA-BB)-(p0+b/d))==0)
check('Eq30 nonnegative initial condition identity',s.simplify(P.subs(t,t1)-AA**2/(AA-BB)*(s.exp(d*t1)-1)-p0)==0)
y1,y2,y3=s.symbols('Y1 Y2 Y3',real=True);A=y2-y1;B=y3-y2
dr=(s.log(A)-s.log(B))/h
grad=s.Matrix([s.diff(dr,z) for z in [y1,y2,y3]])
expected=s.Matrix([-1/A,1/A+1/B,-1/B])/h
check('Eq31 derivative gradient',all(s.simplify(x)==0 for x in grad-expected))
check('Constant additive offset cancels',s.simplify(sum(grad))==0)
check('Eq31 equal independent error variance',s.simplify((grad.T*grad)[0]-2/h**2*(1/A**2+1/B**2+1/(A*B)))==0)

def G(tt,dd):return np.asarray(tt) if dd==0 else -np.expm1(-dd*np.asarray(tt))/dd
recon_error=0.;corner_cases=0
for p,b0,dd,tt,hh in product([0.,.015],[.001,.004],[.01,.04,.2],[0.,5.],[5.,20.]):
 yy=p+b0*G(np.array([tt,tt+hh,tt+2*hh]),dd);a1=yy[1]-yy[0];a2=yy[2]-yy[1]
 dd1=math.log(a1/a2)/hh;bb1=dd1*a1*math.exp(dd1*tt)/(-math.expm1(-dd1*hh));p1=yy[0]-bb1*G(tt,dd1)
 recon_error=max(recon_error,abs(dd1-dd),abs(bb1-b0),abs(p1-p))
 eps=.2*min(a2/2,(a1-a2)/4)
 lo=math.log((a1-2*eps)/(a2+2*eps))/hh;hi=math.log((a1+2*eps)/(a2-2*eps))/hh
 vals=[]
 for signs in product([-1,1],repeat=3):
  zz=yy+eps*np.array(signs);v=math.log((zz[1]-zz[0])/(zz[2]-zz[1]))/hh;vals.append(v);corner_cases+=1
 check('Eq32 sharp corners case '+str(corner_cases),abs(min(vals)-lo)<2e-10 and abs(max(vals)-hi)<2e-10)
check('Three-point reconstruction deterministic grid',recon_error<1e-9)
xx=s.symbols('x',positive=True);FF=2*(1+s.exp(xx)+s.exp(2*xx))/(xx**2*(1-s.exp(-xx))**2)
df=s.lambdify(xx,s.diff(FF,xx),'numpy');opt=brentq(df,1.,2.5)
check('Reported local spacing optimum',abs(opt-1.59041)<.000005 and abs(opt/.04-39.76)<.005)
fi=t*t*s.exp(-2*d*t)
check('Eq34 information derivative',s.simplify(s.diff(fi,t)-2*t*s.exp(-2*d*t)*(1-d*t))==0)
details['reconstruction']=dict(max_absolute_error=recon_error,hard_bound_corners=corner_cases,spacing_x=opt,spacing_minutes=opt/.04)

u,c,rho,gg=s.symbols('u c rho G',positive=True)
check('Pool implicit derivative gives ODE',s.simplify((-c*s.exp(-d*t))/(1/u+rho)+c*s.exp(-d*t)*u/(1+rho*u))==0)
max_pool_error=0.;pool_cases=0
grid=np.linspace(0,25,101)
for cc,dd,rr in product([.01,.4,1.1],[0.,.12,.7],[0.,.25,1.5,12.]):
 g=G(grid,dd);exact=np.exp(-cc*g) if rr==0 else wrightomega(np.log(rr)+rr-cc*g)/rr
 sol=solve_ivp(lambda tt,z:-cc*np.exp(-dd*tt)*z/(1+rr*z),[0,25],[1.],t_eval=grid,method='DOP853',rtol=2e-12,atol=2e-13)
 assert sol.success
 max_pool_error=max(max_pool_error,float(np.max(np.abs(exact-sol.y[0]))));pool_cases+=1
 check('Pool range and monotonicity '+str(pool_cases),np.all(exact>0) and np.all(exact<=1+1e-12) and np.all(np.diff(exact)<=1e-12))
check('Pool closed form versus new ODE grid',max_pool_error<1e-9)
pp=.015+.280*(1-wrightomega(np.log(1.5)+1.5-1.1*G(np.array([.5,1,2]),.12))/1.5)
check('Chapter4 printed pool example',np.max(np.abs(pp-np.array([.07209,.12037,.19251])))<.0000051)
v_inf=float(wrightomega(np.log(1.5)+1.5-1.1/.12)/1.5)
check('Chapter4 positive residual pool limit',abs(v_inf-.000467847)<5.1e-10)
deadline=-math.log(1-.12*.280/.1232)/.12
check('Chapter5 finite pool horizon example',abs(deadline-2.654)<.00051)
details['pool']=dict(cases=pool_cases,max_absolute_ode_error=max_pool_error,remaining_pool_limit=v_inf,finite_horizon_hours=deadline,rho_zero='Continuous mathematical limit only, not finite S0/K with S0>0 and K finite.')

v,a,eta,k,e,lam,mu,ar,an=s.symbols('v a eta k e lambda mu ar an',positive=True)
obs=s.Matrix([a*v,a*(v+eta*e),k*v]);T={v:lam*v,a:a/lam,eta:lam*eta,k:k/lam}
check('Common channels full curve invariant',all(s.simplify(x)==0 for x in obs.subs(T,simultaneous=True)-obs))
J=obs.jacobian([v,a,eta,k]);check('Common channels positive rank-three minor',s.simplify(J[:,[0,2,3]].det()-a*a*e*v)==0)
N,D,E=s.symbols('N D E',positive=True)
fiber={a:N/v,eta:D*v/(e*N),k:E/v}
check('Common channels inverse family exhaustive by substitution',all(s.simplify(x)==0 for x in obs.subs(fiber,simultaneous=True)-s.Matrix([N,N+D,E])))
obs2=s.Matrix([an*v,ar*eta*e,k*v]);T2={v:lam*v,an:an/lam,ar:ar/mu,eta:mu*eta,k:k/lam}
check('Separate channels two invariances',all(s.simplify(x)==0 for x in obs2.subs(T2,simultaneous=True)-obs2))
J2=obs2.jacobian([v,an,ar,eta,k]);check('Separate channels positive rank-three minor',s.simplify(J2[:,[0,2,4]].det()-an*eta*e*v)==0)
fiber2={an:N/v,ar:D/(eta*e),k:E/v}
check('Separate inverse family substitution',all(s.simplify(x)==0 for x in obs2.subs(fiber2,simultaneous=True)-s.Matrix([N,D,E])))
rn=s.symbols('rn',positive=True)
cal={v:rn*eta*e*N/D,an:D/(rn*eta*e),ar:D/(eta*e),k:E*D/(rn*eta*e*N)}
check('Separate calibrated reference inversion',all(s.simplify(x)==0 for x in obs2.subs(cal,simultaneous=True)-s.Matrix([N,D,E])))
check('Separate calibrated native assay inversion',s.simplify((an*v).subs({v:E/k,an:N*k/E},simultaneous=True)-N)==0)
check('Effective factor reparametrization',s.simplify(an*(v+ar/an*eta*e)-(an*v+ar*eta*e))==0)
chi=s.symbols('chi',nonnegative=True)
check('Native change can preserve linear dose response',s.expand((1+chi*e)*N+ar*eta*e-(N+(ar*eta+N*chi)*e))==0)
for ev in [0,.001,.02,.1]:
 check('Native constancy counterexample dose '+str(ev),abs((.004+.4*ev)-((1+50*ev)*.004+.2*ev))<1e-15)

zz,ll=s.symbols('Z L',positive=True);dd=zz-N
av=dd/ll;vv=ll*N/dd
g1=s.Matrix([s.diff(av,x) for x in [N,zz,ll]]);g2=s.Matrix([s.diff(vv,x) for x in [N,zz,ll]])
check('Eq62 reference covariance gradients',all(s.simplify(x)==0 for x in g1-s.Matrix([-1/ll,1/ll,-dd/ll**2])) and all(s.simplify(x)==0 for x in g2-s.Matrix([ll*zz/dd**2,-ll*N/dd**2,N/dd])))
lo,hi=s.symbols('lo hi',positive=True);hmean=2*lo*hi/(lo+hi)
check('Eq78 minimax endpoint equality',s.simplify(hmean/lo-1-(1-hmean/hi))==0)
check('Eq78 minimax maximum formula',s.simplify(hmean/lo-1-(hi-lo)/(hi+lo))==0)

labdata=rows(R/'externe_pruefung/aktivitaeten.csv');assay={};labmetrics={}
for kind in ['identity','common','product']:
 residual=[];relative=[];perlab={}
 for lab in sorted({x['lab'] for x in labdata}):
  train=[x for x in labdata if x['lab']!=lab];test=[x for x in labdata if x['lab']==lab]
  assert len(train)==16 and len(test)==4
  rr=[]
  for x in test:
   tr=train if kind=='common' else [z for z in train if z['product']==x['product']]
   logfactor=0 if kind=='identity' else np.mean([math.log(float(z['activity_37'])/float(z['activity_20'])) for z in tr])
   pred=float(x['activity_20'])*math.exp(logfactor);err=math.log(pred/float(x['activity_37']))
   rr.append(err);residual.append(err);relative.append(abs(pred/float(x['activity_37'])-1))
  perlab[lab]=float(np.sqrt(np.mean(np.square(rr))))
 assay[kind]=dict(Elog=float(np.sqrt(np.mean(np.square(residual)))),MARE=100*float(np.mean(relative)))
 labmetrics[kind]=perlab
tab=tables[13]['rows']
for kind,row in zip(['identity','common','product'],tab):
 check('Table13 '+kind,abs(assay[kind]['Elog']-float(row[2].replace(',','.')))<.0000051 and abs(assay[kind]['MARE']-float(row[3].replace('%','').replace(',','.')))<.0051)
for row in tables[14]['rows']:
 check('Table14 lab '+row[0],all(abs(labmetrics[k][row[0]]-float(row[j].replace(',','.')))<.0000051 for k,j in [('common',1),('product',2)]))
details['assay']=dict(metrics=assay,lab_Elog=labmetrics)

bench=read(R/'schlussfassung28/modellpruefung.json')['results'];curve={'Shintani_A':'Shintani A','Shintani_F':'Shintani F','H_STD':'MSM 4,78','H_HIGH':'MSM 8,37'}
for item,row in zip(bench,tables[35]['rows']):
 check('Table35 identity '+item['curve']+item['model'],row[0].split('\n')[0]==curve[item['curve']] and row[1]==item['model'])
 check('Table35 arithmetic '+item['curve']+item['model'],all(abs(float(row[j].replace(',','.'))-item[key])<.0000051 for j,key in [(2,'early_RMSE'),(3,'later_RMSE'),(4,'full_numerical_minimax')]))

certs=rows(R/'zertifizierte_fehlergrenzen/results/certificates.csv');data=rows(R/'zertifizierte_fehlergrenzen/data/nominal_data.csv')
for cert,row in zip(certs,tables[33]['rows']):
 yy=[F(x['value']) for x in data if x['curve']==cert['curve']]
 floor=max([F(0),*[yy[i]-yy[j] for i in range(len(yy)) for j in range(i+1,len(yy))]])/2
 check('Monotone floor '+cert['curve'],F(cert['monotone_floor_lower'])<=floor<=F(cert['monotone_floor_upper']))
 upper=float(cert['upper']);printed=float(row[1].replace('≈','').replace(',','.'))
 check('Table33 certified value '+cert['curve'],abs(printed-upper)<5.1e-11)
 gain=100*(1-float(floor)/upper);shown=float(row[3].replace('%','').replace(',','.'))
 check('Table33 outward improvement '+cert['curve'],gain<=shown+1e-12 and shown-gain<.010001)

L=read(R/'quellen/teilmodell_parameter.json')['liquefaction']
term=(1-F(str(L['h'][2])))*F(str(L['S0']))-F(str(L['x0'][3]))
check('Eq82 exact negative source term',term==F('-.1174246'))
boundlo=(1-F('-.3220105'))*F('161.35')-F('213.495')
boundhi=(1-F('-.3220115'))*F('161.45')-F('213.485')
check('Chapter8 source rounding sign',boundlo==F('-.188605825') and boundhi==F('-.046243325') and boundhi<0)
details['source_term']=dict(exact=str(term),rounding=[str(boundlo),str(boundhi)],scope='Conditional transcribed source term. No reconstruction of original full fit.')

# Independent unit changes in E, W and the minimal model.
tt=np.linspace(0,4,15)
check('Minimal model hours to minutes',np.max(np.abs(.004*G(tt,.04)-(.004/60)*G(60*tt,.04/60)))<1e-14)
for shape in [1.,.7,2.]:
 orig=.01+.2*(-np.expm1(-(.8*tt)**shape));converted=10+200*(-np.expm1(-((.8/60)*(60*tt))**shape))
 check('E W consistent time and product units '+str(shape),np.max(np.abs(converted/1000-orig))<1e-14)

out=dict(version=36,public_package_version='1.0',status='PASS',check_count=len(checks),checks=checks,details=details,scope='Neu formulierte interne Symbolik, ODE- und Tabellenkontrollen. Keine neue Messung und keine personell unabhängige Begutachtung. Die Alternationsbegründung und Vollständigkeit der positiven Familien sind zusätzlich im Abnahmebericht argumentativ geprüft.')
(S/'Formel_und_Tabellenabnahme.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status='PASS',checks=len(checks),pool_max_error=max_pool_error,assay=assay),ensure_ascii=False))
