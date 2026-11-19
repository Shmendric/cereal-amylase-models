"""Härtetest getrennter Zugänglichkeiten: Algebra, Gegenfälle, exakte Bereiche.

Alle Zahlen sind gesetzte Modellfälle, keine Getreidemessungen.
"""
from pathlib import Path
import json,csv,sys,hashlib
import numpy as np
import sympy as s
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parent
O=R/'ergebnisse';O.mkdir(exist_ok=True)
checks=[]
def check(name,condition,detail=''):
    if not bool(condition):raise AssertionError(name)
    checks.append(dict(name=name,status='PASS',detail=detail))

v,an,ar,eta,k,e,lam,mu=s.symbols('v a_N a_R eta kappa e lambda mu',positive=True)
n,d,ex,rho,L=s.symbols('N D E rho L',positive=True)
theta=s.Matrix([v,an,ar,eta,k])
f=s.Matrix([an*v,ar*eta*e,k*v])
J=f.jacobian(theta)
trans={v:lam*v,an:an/lam,ar:ar/mu,eta:mu*eta,k:k/lam}
check('Zwei unabhängige positive Skalierungen',all(s.simplify(x)==0 for x in f.subs(trans,simultaneous=True)-f))
u=s.Matrix([v,-an,0,0,-k]);w=s.Matrix([0,0,-ar,eta,0])
check('Erste Jacobi-Nullrichtung',s.simplify(J*u)==s.zeros(3,1))
check('Zweite Jacobi-Nullrichtung',s.simplify(J*w)==s.zeros(3,1))
check('Unabhängige Nullrichtungen',s.Matrix.hstack(u,w).rank()==2)
check('Rang drei',J.rank()==3)
check('Positiver Dreierminor',s.simplify(J[:,[0,2,4]].det()-an*eta*e*v)==0)
fiber={an:n/v,ar:d/(eta*e),k:ex/v}
check('Vollständige positive Lösungsfamilie erfüllt Beobachtungen',all(s.simplify(x)==0 for x in f.subs(fiber,simultaneous=True)-s.Matrix([n,d,ex])))
check('Separate Referenzwirkung fixiert Referenzzugang',s.simplify((ar-d/(eta*e)).subs(fiber,simultaneous=True))==0)
check('Bekanntes eta lässt native Skalierung frei',all(s.simplify(x)==0 for x in f.subs({v:lam*v,an:an/lam,k:k/lam},simultaneous=True)-f))
check('Kalibrierter Assay bestimmt native Faktoren',s.simplify((k*v).subs(v,ex/k)-ex)==0 and s.simplify((an*v).subs({v:ex/k,an:n*k/ex},simultaneous=True)-n)==0)
check('Bekanntes kappa lässt Referenzaufteilung frei',all(s.simplify(x)==0 for x in f.subs({ar:ar/mu,eta:mu*eta},simultaneous=True)-f))
inv={v:rho*eta*e*n/d,an:d/(rho*eta*e),ar:d/(eta*e),k:ex*d/(rho*eta*e*n)}
check('Inverse bei bekannter Referenzwirkung und Verhältnis',all(s.simplify(x)==0 for x in f.subs(inv,simultaneous=True)-s.Matrix([n,d,ex])))
check('Verhältnis der inversen Faktoren',s.simplify((ar/an).subs(inv,simultaneous=True)-rho)==0)
eff=ar*eta/an
check('Altes Modell enthält effektives Produkt rho eta',s.simplify(an*(v+eff*e)-(an*v+ar*eta*e))==0)
va=eta*e*(an*v)/(ar*eta*e);aa=(ar*eta*e)/(eta*e)
check('Naiver Ratenfehler ist inverses Zugangsverhältnis',s.simplify(va/v-an/ar)==0)
check('Naiver Zugangsfehler ist Zugangsverhältnis',s.simplify(aa/an-ar/an)==0)
check('Naives Produkt bleibt richtig',s.simplify(va*aa-an*v)==0)
gn,gr,ge=s.symbols('G_N G_R G_E')
curve=s.Matrix([an*v*gn,an*v*gn+ar*eta*e*gr,k*v*ge])
check('Bekannte Formfunktionen beseitigen beide Skalierungen nicht',all(s.simplify(x)==0 for x in curve.subs(trans,simultaneous=True)-curve))

rng=np.random.default_rng(2026111122)
max_inverse=0.;max_orbit=0.;cases=400
for _ in range(cases):
    vv,ee=10**rng.uniform(-3,-.3,2);a,b,hh,kk=rng.uniform(.05,.9,4)
    nn=a*vv;dd=b*hh*ee;eeobs=kk*vv;rr=b/a
    reco=np.array([rr*hh*ee*nn/dd,dd/(rr*hh*ee),dd/(hh*ee),eeobs*dd/(rr*hh*ee*nn)])
    max_inverse=max(max_inverse,float(np.max(np.abs(reco/np.array([vv,a,b,kk])-1))))
    lv,mv=10**rng.uniform(-1,1,2)
    obs2=np.array([(a/lv)*(lv*vv),(b/mv)*(mv*hh)*ee,(kk/lv)*(lv*vv)])
    max_orbit=max(max_orbit,float(np.max(np.abs(obs2/np.array([nn,dd,eeobs])-1))))
check('400 numerische Inversionen',max_inverse<1e-12,str(max_inverse))
check('400 numerische Skalierungsfälle',max_orbit<1e-12,str(max_orbit))

example=[]
for name,vv,a,b,kk in [('A',.01,.4,.2,.8),('B',.02,.2,.2,.4)]:
    example.append(dict(state=name,v0=vv,a_N=a,a_R=b,eta=1.,e=.02,kappa=kk,rho=b/a,N=a*vv,D=b*.02,Z=a*vv+b*.02,E=kk*vv))
check('Bekanntes eta und identische drei Signale trotz verschiedener nativer Parameter',all(abs(example[0][key]-example[1][key])<1e-15 for key in ['N','Z','E','eta','e']))
check('Gegenfall auch bei Zugängen und Wiederfindung höchstens eins',all(0<x['a_N']<=1 and 0<x['a_R']<=1 and 0<x['kappa']<=1 for x in example),'q = 1 U/R, daher kappa numerisch gleich r0')

nn,dd,eo,ee=.004,.004,.008,.02
rmin,rmax=.8,1.25;ell=.02
vref=ell*nn/dd
known_ratio=dict(rho=[rmin,rmax],ell=ell,v0=[rmin*vref,rmax*vref],a_N=[dd/(rmax*ell),dd/(rmin*ell)],naive_v0=vref,
                 relative_naive_error=[1/rmax-1,1/rmin-1])
check('Exakte Verhältnisintervalle konstruktiv',all(abs((dd/(rr*ell))*(rr*ell*nn/dd)-nn)<1e-15 for rr in np.linspace(rmin,rmax,251)))
tol=.1;rlo=1/(1+tol);rhi=1/(1-tol)
check('Zehnprozentbedingung an beiden Grenzen exakt',abs((1/rlo-1)-tol)<1e-14 and abs((1/rhi-1)+tol)<1e-14)
check('Außerhalb der Zehnprozentgrenzen kann Ziel verletzt sein',1/(rlo*.99)-1>tol and 1-1/(rhi*1.01)>tol)

# Minimax hier ausschließlich für relativen Fehler von v0 bei bekanntem ell und fehlerfreien N,D.
rr0=2*rmin*rmax/(rmin+rmax);bound=(rmax-rmin)/(rmax+rmin)
worst=lambda x:max(abs(x/rmin-1),abs(x/rmax-1))
numeric=minimize_scalar(worst,bounds=(rmin,rmax),method='bounded',options={'xatol':1e-14})
check('Minimax-Kandidat gleicht Randfehler aus',abs((rr0/rmin-1)-(1-rr0/rmax))<1e-14)
check('Minimax-Wert stimmt mit beiden Rändern überein',abs(worst(rr0)-bound)<1e-14)
check('Numerische Gegenoptimierung des relativen Minimaxwertes',numeric.success and abs(numeric.fun-bound)<1e-7 and abs(numeric.x-rr0)<1e-7)
minimax=dict(rho0=rr0,worst_relative_error=bound,v0_point=rr0*vref,criterion='Maximaler relativer Fehler von v0; keine optimale Enzymdosierung',numeric_optimizer=float(numeric.x))

eta_min,eta_max=.8,1.2;kmin,kmax=.4,.8;C=nn*ee/dd
vl=max(C*rmin*eta_min,eo/kmax);vu=min(C*rmax*eta_max,eo/kmin)
construction=[]
for vv in np.linspace(vl,vu,401):
    hhprod=vv/C
    lo=max(eta_min,hhprod/rmax);hi=min(eta_max,hhprod/rmin)
    hh=np.sqrt(lo*hi);rr=hhprod/hh;a=nn/vv;b=dd/(hh*ee);kk=eo/vv
    construction.append(dict(v0=float(vv),a_N=float(a),a_R=float(b),eta=float(hh),rho=float(rr),kappa=float(kk)))
    if not (lo<=hi+1e-13 and eta_min-1e-13<=hh<=eta_max+1e-13 and rmin-1e-13<=rr<=rmax+1e-13 and kmin-1e-13<=kk<=kmax+1e-13):raise AssertionError('Konstruktion Intervall')
    if max(abs(a*vv-nn),abs(b*hh*ee-dd),abs(kk*vv-eo))>1e-14:raise AssertionError('Konstruktion Signale')
check('401 Konstruktionen im exakten Schnittbereich',abs(vl-.0128)<1e-14 and abs(vu-.02)<1e-14)
check('Engere unvereinbare Kalibrierung ergibt leeren Schnitt',max(C*.8*.8,eo/.2)>min(C*1.25*1.2,eo/.1))
partial=dict(N=nn,D=dd,E=eo,e=ee,rho=[rmin,rmax],eta=[eta_min,eta_max],kappa=[kmin,kmax],v0=[vl,vu],a_N=[nn/vu,nn/vl],construction_points=len(construction))

# Sichere äußere Grenzen mit gemeinsamem Eingang N in D = Z-N.
zz=nn+dd;dn=dz=.0002;dl=.001;delta=dn+dz
avbounds=[(dd-delta)/(rmax*(ell+dl)),(dd+delta)/(rmin*(ell-dl))]
vvbounds=[rmin*(ell-dl)*(nn-dn)/(dd+delta),rmax*(ell+dl)*(nn+dn)/(dd-delta)]
draws=rng.uniform([-dn,-dz,-dl,rmin],[dn,dz,dl,rmax],size=(20000,4))
nt=nn+draws[:,0];zt=zz+draws[:,1];lt=ell+draws[:,2];rt=draws[:,3];dt=zt-nt
vv=rt*lt*nt/dt;aa=dt/(rt*lt)
check('20000 gemeinsame Fehlereingänge innerhalb sicherer Grenzen',np.all((aa>=avbounds[0])&(aa<=avbounds[1])&(vv>=vvbounds[0])&(vv<=vvbounds[1])))
check('Positive Nenner der harten Grenzen',nn>dn and dd>delta and ell>dl and rmin>0)
error_bounds=dict(delta_N=dn,delta_Z=dz,delta_ell=dl,rho=[rmin,rmax],a_N=avbounds,v0=vvbounds,numerical_cases=20000,scope='Harte gesetzte Grenzen; keine Standardunsicherheiten oder Konfidenzintervalle')

chi=s.symbols('chi',nonnegative=True)
zpert=(1+chi*e)*an*v+ar*eta*e
check('Native Veränderung kann perfekte Dosisgerade ergeben',s.expand(zpert-(an*v+(ar*eta+an*v*chi)*e))==0)
doses=np.array([0,.0025,.005,.01,.02]);counter=[]
for state,vv,a,ch,kk in [('C',.01,.4,0.,.8),('D',.02,.2,50.,.4)]:
    zzs=(1+ch*doses)*a*vv+a*doses
    counter.append(dict(state=state,v0=vv,a_N=a,a_R=a,eta=1.,rho=1.,chi=ch,kappa=kk,N=a*vv,E=kk*vv,doses=doses.tolist(),Z=zzs.tolist()))
check('Bekanntes rho und eta genügen bei ungesicherter nativer Konstanz nicht',np.allclose(counter[0]['Z'],counter[1]['Z'],rtol=0,atol=1e-15) and counter[0]['v0']!=counter[1]['v0'])
check('Gegenfall zeigt verschiedene native Faktoren bei gleichem Extraktsignal',abs(counter[0]['E']-counter[1]['E'])<1e-15)
check('Native Modifikation bleibt im Gegenfall positiv',np.all(1+50*doses>0))

result=dict(status='PASS',check_count=len(checks),checks=checks,numerical_inverse_cases=cases,seed=2026111122,example=example,known_ratio=known_ratio,
            ten_percent_ratio=[rlo,rhi],minimax=minimax,partial_identification=partial,error_bounds=error_bounds,dose_counterexample=counter,
            scope='Getrennte operationelle Zugänglichkeiten und ausdrücklich definierte Beobachtungsmodelle; synthetische Zahlen, keine Laborvalidierung.',
            versions=dict(python=sys.version.split()[0],numpy=np.__version__,sympy=s.__version__,matplotlib=matplotlib.__version__))
(O/'pruefergebnisse.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
with (O/'intervall_konstruktionen.csv').open('w',newline='') as f0:
    cw=csv.DictWriter(f0,fieldnames=list(construction[0]));cw.writeheader();cw.writerows(construction)
with (O/'gegenfaelle.csv').open('w',newline='') as f0:
    cw=csv.DictWriter(f0,fieldnames=list(example[0]));cw.writeheader();cw.writerows(example)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(9,3.8))
rr=np.linspace(.5,1.5,301)
ax[0].plot(rr,1/rr,label=r'$\widehat v_0/v_0$',color='#18617c')
ax[0].plot(rr,rr,label=r'$\widehat a/a_N$',color='#b15b32')
ax[0].axhline(1,color='.6',lw=.8);ax[0].axvspan(.8,1.25,color='.85',alpha=.4)
ax[0].set(xlabel=r'Zugangsverhältnis $\rho=a_R/a_N$',ylabel='Verhältnis',title='Fehler bei angenommener Gleichheit');ax[0].legend(frameon=False,fontsize=10)
rr0s=np.linspace(.65,1.4,301)
ax[1].plot(rr0s,[100*worst(x) for x in rr0s],color='#18617c');ax[1].scatter([rr0],[100*bound],color='#b15b32',zorder=4)
ax[1].axvline(1,color='.6',lw=.8,ls='--');ax[1].set(xlabel='Für die Auswertung gesetztes ρ₀',ylabel='Größter relativer Ratenfehler / %',title='Nur für 0,8 ≤ ρ ≤ 1,25')
fig.tight_layout();fig.savefig(O/'zugangsverhaeltnis.png',dpi=190);plt.close(fig)
print(json.dumps(dict(status='PASS',checks=len(checks),known_ratio=known_ratio,partial=partial,minimax=minimax),ensure_ascii=False,indent=2))
