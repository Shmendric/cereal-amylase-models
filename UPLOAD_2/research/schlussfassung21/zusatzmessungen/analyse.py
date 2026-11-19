"""Reproduzierbare mathematische und synthetische Prüfung der Zusatzmessmodelle."""
from pathlib import Path
import json,csv,sys,platform,hashlib
import numpy as np
import sympy as s
import scipy
from scipy.special import ndtr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
out=ROOT/'ergebnisse';out.mkdir(exist_ok=True)
checks=[]
def check(name,condition,detail=''):
 if not bool(condition):raise RuntimeError(name)
 checks.append({'name':name,'status':'PASS','detail':detail})
v,a,eta,kappa,e,lam=s.symbols('v a eta kappa e lam',positive=True)
N,Z,E,L=s.symbols('N Z E L',positive=True)
D=Z-N
f=s.Matrix([a*v,a*(v+eta*e),kappa*v]);theta=s.Matrix([v,a,eta,kappa]);J=f.jacobian(theta)
trans={v:lam*v,a:a/lam,eta:lam*eta,kappa:kappa/lam}
check('Skalierungssymmetrie',all(s.simplify(x)==0 for x in f.subs(trans,simultaneous=True)-f))
null=s.Matrix([v,-a,eta,-kappa]);check('Jacobi Nullrichtung',s.simplify(J*null)==s.zeros(3,1))
check('Generischer Rang drei',J.rank()==3)
Ji=s.Matrix([a*v,a*(v+L)]).jacobian([v,a]);check('Idealer Rang zwei',s.simplify(Ji.det()-a*L)==0)
inv={a:D/L,v:L*N/D};check('Globale ideale Inverse',all(s.simplify(x)==0 for x in s.Matrix([a*v,a*(v+L)]).subs(inv,simultaneous=True)-s.Matrix([N,Z])))
cal_k={v:E/kappa,a:N*kappa/E,eta:D*E/(e*N*kappa)}
check('Inverse bei bekanntem kappa',all(s.simplify(x)==0 for x in f.subs(cal_k,simultaneous=True)-s.Matrix([N,Z,E])))
cal_eta={v:eta*e*N/D,a:D/(eta*e),kappa:E*D/(eta*e*N)}
check('Inverse bei bekanntem eta',all(s.simplify(x)==0 for x in f.subs(cal_eta,simultaneous=True)-s.Matrix([N,Z,E])))
# Vollständige Kurven: bekannte oder separat parametrisierte Formfunktionen.
gn,gr,ge=s.symbols('G_N G_R G_E');curves=s.Matrix([a*v*gn,a*v*gn+a*eta*e*gr,kappa*v*ge])
check('Zeitkurven Symmetrie',all(s.simplify(x)==0 for x in curves.subs(trans,simultaneous=True)-curves))
m0,m1=s.symbols('m0 m1');check('Dosisgerade',s.simplify(a*(v+eta*e)-(a*v+a*eta*e))==0)
# Unabhängige symbolische Ableitung und handnotierte Gradienten.
ga=s.Matrix([-1/L,1/L,-D/L**2]);gv=s.Matrix([L*Z/D**2,-L*N/D**2,N/D])
check('Gradient a',s.simplify(s.Matrix([s.diff(D/L,x) for x in [N,Z,L]])-ga)==s.zeros(3,1))
check('Gradient v',s.simplify(s.Matrix([s.diff(L*N/D,x) for x in [N,Z,L]])-gv)==s.zeros(3,1))
check('Gemeinsamer additiver Fehler',s.simplify(ga[0]+ga[1])==0 and s.simplify(gv[0]+gv[1]-L/D)==0)
check('Gemeinsamer multiplikativer Fehler',s.simplify(gv[0]*N+gv[1]*Z)==0 and s.simplify(ga[0]*N+ga[1]*Z-D/L)==0)
# Numerische Rundreise auf deterministischem Raster.
rng=np.random.default_rng(20261111)
maxerr=0.; ngrid=300
for _ in range(ngrid):
 vv,aa,ee,hh=10**rng.uniform(-2.5,-.5),rng.uniform(.1,.9),10**rng.uniform(-2.5,-.5),rng.uniform(.3,1.5)
 ll=ee*hh;nn=aa*vv;zz=aa*(vv+ll)
 recovered=np.array([ll*nn/(zz-nn),(zz-nn)/ll]); truth=np.array([vv,aa])
 maxerr=max(maxerr,float(np.max(np.abs(recovered/truth-1))))
check('Dreihundert numerische Inversionen',maxerr<1e-12,str(maxerr))
# Zwei physikalisch verschiedene synthetische Zustände mit identischen Signalen.
example=[]
for vv,aa,hh,rr in [(.01,.4,.5,.8),(.02,.2,1.,.4)]:
 example.append(dict(v=vv,a=aa,eta=hh,r=rr,q=1.,e=.02,N=aa*vv,Z=aa*(vv+hh*.02),E=rr*vv))
check('Gegenbeispiel',all(abs(example[0][k]-example[1][k])<1e-15 for k in ['N','Z','E']))
# Partielle Identifizierung aus harten Kalibrierintervallen, fehlerfreie Signale.
nn,zz,ee,obs=.004,.008,.02,.008
eta_bounds=(.5,1.);k_bounds=(.4,.8);C=nn*ee/(zz-nn)
vl=max(C*eta_bounds[0],obs/k_bounds[1]);vu=min(C*eta_bounds[1],obs/k_bounds[0]);al,au=nn/vu,nn/vl
for vv in np.linspace(vl,vu,101):
 aa=nn/vv;hh=(zz-nn)*vv/(nn*ee);kk=obs/vv
 check_cond=eta_bounds[0]-1e-14<=hh<=eta_bounds[1]+1e-14 and k_bounds[0]-1e-14<=kk<=k_bounds[1]+1e-14
 if not check_cond:raise RuntimeError('Intervallkonstruktion')
check('Intervall ist exakt konstruierbar',abs(vl-.01)<1e-15 and abs(vu-.02)<1e-15)
# Monte Carlo ist ausschließlich ein konstruiertes Fehlerexperiment.
sigma=.0002; rel_L=.05; M=200000; scenarios=[]
rng=np.random.default_rng(20261112)
for ll in [.0005,.002,.01,.02,.04]:
 aa=.4; vv=.01; nn=aa*vv; zz=aa*(vv+ll); dd=zz-nn
 G=np.array([[-1/ll,1/ll,-dd/ll**2],[ll*zz/dd**2,-ll*nn/dd**2,nn/dd]])
 cov=np.diag([sigma**2,sigma**2,(ll*rel_L)**2]);covout=G@cov@G.T
 x=rng.normal([nn,zz,ll],[sigma,sigma,ll*rel_L],size=(M,3));d=x[:,1]-x[:,0]
 valid=(x[:,0]>0)&(d>0)&(x[:,2]>0)
 ar=d[valid]/x[valid,2];vr=x[valid,2]*x[valid,0]/d[valid]
 row=dict(L=ll,N=nn,Z=zz,D=dd,signal_to_noise_D=dd/(2**.5*sigma),invalid=int((~valid).sum()),draws=M,invalid_percent=100*float((~valid).mean()),delta_rel_sd_a=100*float(np.sqrt(covout[0,0])/aa),delta_rel_sd_v=100*float(np.sqrt(covout[1,1])/vv),conditional_a_quantiles=np.quantile(ar,[.025,.5,.975]).tolist(),conditional_v_quantiles=np.quantile(vr,[.025,.5,.975]).tolist())
 scenarios.append(row)
check('Synthetischer schwacher Kontrast',scenarios[0]['invalid_percent']>20 and scenarios[3]['invalid']==0)
# Ein endlicher additiver Versatz: Differenz bleibt, die v-Rekonstruktion ändert sich.
nn,zz,ll=.004,.012,.02;shift=.001
check('Endlicher additiver Gegenfall',abs(((zz+shift)-(nn+shift))/ll-.4)<1e-14 and abs(ll*(nn+shift)/(zz-nn)-.0125)<1e-14)
# Fehlergrenzen: garantierter endlicher Bereich, wenn D - delta_D positiv bleibt.
dn,dz,dl=.0002,.0002,.001;dd=zz-nn;deltaD=dn+dz
alo=(dd-deltaD)/(ll+dl);ahi=(dd+deltaD)/(ll-dl)
vlo=(ll-dl)*(nn-dn)/(dd+deltaD);vhi=(ll+dl)*(nn+dn)/(dd-deltaD)
check('Konservative Fehlergrenzen',alo<.4<ahi and vlo<.01<vhi)
# Lokaler Ranggewinn durch kalibrierte Zusatzgröße.
Jcal=f.col_join(s.Matrix([kappa])).jacobian(theta)
check('Kalibrierkanal hebt Rang auf vier',Jcal.rank()==4)
# Positive Gewichtung kann eine Nullrichtung nicht aufheben.
w=s.symbols('w1:4',positive=True);I=J.T*s.diag(*w)*J
check('Informationsmatrix bleibt singulär',s.simplify(I*null)==s.zeros(4,1))
# Zusätzliche Kontrollen der im Text verwendeten Rang- und Skalenargumente.
check('Positiver Dreierminor',s.simplify(J[:,[0,2,3]].det()-a**2*e*v)==0)
fiber={a:N/v,eta:D*v/(e*N),kappa:E/v}
check('Vollständige Faser erfüllt alle Signale',all(s.simplify(x)==0 for x in f.subs(fiber,simultaneous=True)-s.Matrix([N,Z,E])))
check('Invariantes Kalibrierprodukt',s.simplify((eta*kappa-D*E/(e*N)).subs({N:a*v,Z:a*(v+eta*e),E:kappa*v},simultaneous=True))==0)
pp=s.symbols('p1:5',real=True);H=v**pp[0]*a**pp[1]*eta**pp[2]*kappa**pp[3]
DH=v*s.diff(H,v)-a*s.diff(H,a)+eta*s.diff(H,eta)-kappa*s.diff(H,kappa)
check('Potenzkriterium',s.simplify(DH/H-sum([pp[0],-pp[1],pp[2],-pp[3]]))==0)
for scale in [.001,1000.]:
 n2,z2,l2=.004*scale,.012*scale,.02*scale
 check('Einheitenwechsel '+str(scale),abs(((z2-n2)/l2)-.4)<1e-14 and abs((l2*n2/(z2-n2))/scale-.01)<1e-14)
# Abgesicherte äußere Grenzen gegen die acht Ecken der Eingangsbücher.
import itertools
for n3,z3,l3 in itertools.product([nn-dn,nn+dn],[zz-dz,zz+dz],[ll-dl,ll+dl]):
 av=(z3-n3)/l3;vvv=l3*n3/(z3-n3)
 if not (alo-1e-14<=av<=ahi+1e-14 and vlo-1e-14<=vvv<=vhi+1e-14):raise RuntimeError('Grenzenprüfung')
check('Fehlergrenzen an acht Eingangsecken',True)
result=dict(status='PASS',scope='Definierte Beobachtungsmodelle und synthetische Fehlerfälle; keine experimentellen Getreidedaten',checks=checks,check_count=len(checks),example=example,partial_identification=dict(v=[vl,vu],a=[al,au],eta=list(eta_bounds),kappa=list(k_bounds)),error_bounds=dict(a=[alo,ahi],v=[vlo,vhi]),scenarios=scenarios,seed_inversion=20261111,seed_monte_carlo=20261112,units='N,Z,L,v in definierter Rateneinheit R; E in Assayeinheit U; a,eta,r dimensionslos; kappa in U/R',versions=dict(python=sys.version.split()[0],numpy=np.__version__,sympy=s.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,platform=platform.platform()))
(out/'pruefergebnisse.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
with (out/'szenarien.csv').open('w',newline='',encoding='utf-8-sig') as f:
 writer=csv.DictWriter(f,fieldnames=[k for k in scenarios[0] if 'quantiles' not in k]);writer.writeheader();writer.writerows([{k:v for k,v in row.items() if 'quantiles' not in k} for row in scenarios])
# Figuren: Trennbarkeit und Präzisionsgrenze, keine empirischen Messkurven.
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(9,3.4));xx=np.linspace(.006,.025,300)
ax[0].plot(xx,.004/xx,label='N = 0,004');ax[0].scatter([.01,.02],[.4,.2],c=['#18617c','#b15b32'],zorder=3);ax[0].set(xlabel='Enzymrate v₀ / R',ylabel='Zugänglichkeit a',title='Identisches natives Signal')
vv=np.linspace(.01,.02,160);ax[1].plot(vv,.004/vv,label='a');ax[1].plot(vv,(.008-.004)*vv/(.004*.02),label='eta');ax[1].plot(vv,.008/vv,label='kappa / (U/R)');ax[1].set(xlabel='Enzymrate v₀ / R',title='Auch drei Kanäle lassen Spielraum');ax[1].legend(frameon=False)
fig.tight_layout();fig.savefig(out/'symmetrie.png',dpi=190);plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(9,3.5));ss=np.array([x['signal_to_noise_D'] for x in scenarios]);ax[0].loglog(ss,[x['delta_rel_sd_a'] for x in scenarios],'o-',label='a');ax[0].loglog(ss,[x['delta_rel_sd_v'] for x in scenarios],'s-',label='v₀');ax[0].set(xlabel='D / Standardunsicherheit von D',ylabel='Lokale relative Standardunsicherheit / %',title='Kleine Differenz verstärkt Fehler');ax[0].legend(frameon=False)
ax[1].plot(ss,[x['invalid_percent'] for x in scenarios],'o-',color='#b15b32');ax[1].set(xlabel='D / Standardunsicherheit von D',ylabel='Nicht auswertbare Ziehungen / %',title='Synthetische Gaußfehler');ax[1].set_xscale('log');fig.tight_layout();fig.savefig(out/'praezision.png',dpi=190);plt.close(fig)
print(json.dumps({'status':result['status'],'checks':len(checks),'intervals':result['partial_identification'],'scenarios':[{k:x[k] for k in ['L','delta_rel_sd_a','delta_rel_sd_v','invalid_percent']} for x in scenarios]},ensure_ascii=False,indent=2))
