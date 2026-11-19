"""Zusätzliche interne Gegenrechnungen, ohne Originaloptimierer zu importieren.

Die LP-Suche liefert numerische zulässige Kandidaten, kein globales Zertifikat.
Die Fehlerbeispiele sind exakt konstruierte Zufallsvariablen, keine Messdaten.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import copy
import csv
import hashlib
import json
import math
import numpy as np
import sympy as sp
from scipy.optimize import linprog, minimize_scalar

S = Path(__file__).resolve().parent
R = S.parent
checks = {}
details = {}

def check(name, value):
    if name in checks:
        raise ValueError('Doppelt vergebener Prüfname: ' + name)
    checks[name] = bool(value)
    if not value:
        raise AssertionError(name)

def expect_rejected(name, operation):
    try:
        operation()
    except (ValueError, AssertionError):
        check(name, True)
    else:
        check(name, False)

# Algebraische Kontrolle der Beobachtungsäquivalenz. Keine historische Neuheit.
t, p0, b, d, A, f, k1, k2 = sp.symbols('t p0 b d A f k1 k2', real=True)
curve = p0 + b / d * (1 - sp.exp(-d*t))
check('Inaktivierung ergibt dieselbe beobachtete Exponentialkurve',
      sp.simplify(sp.diff(curve,t) - b*sp.exp(-d*t)) == 0)
check('Stetiger Grenzfall bei verschwindender Inaktivierung',
      sp.simplify(sp.limit(curve,d,0) - p0 - b*t) == 0)
p2 = p0 + A*(f*(1-sp.exp(-k1*t)) + (1-f)*(1-sp.exp(-k2*t)))
check('Diskrete Vertauschung der Poolbezeichnungen',
      sp.simplify(p2-p2.xreplace({f:1-f,k1:k2,k2:k1})) == 0)
check('Gleiche Raten machen den Anteil unbestimmbar',
      sp.simplify(sp.diff(p2.subs(k2,k1),f)) == 0)
check('Leerer erster Pool lässt seine Rate frei', sp.simplify(sp.diff(p2.subs(f,0),k1)) == 0)
check('Leerer zweiter Pool lässt seine Rate frei', sp.simplify(sp.diff(p2.subs(f,1),k2)) == 0)
check('Verschwindende Amplitude entfernt sämtliche Raten', p2.subs(A,0) == p0)
# Eine volle lokale Sensitivitätsrangzahl an einem inneren Punkt widerlegt die
# falsche Folgerung, allein die diskrete Symmetrie erzwinge überall Rangverlust.
sens = sp.Matrix([[sp.diff(p2,x).subs({p0:1,A:1,f:sp.Rational(1,3),k1:sp.log(2),k2:sp.log(3),t:i})
                   for x in [p0,A,f,k1,k2]] for i in range(5)])
det = sp.simplify(sens.det())
check('Poollabel Symmetrie erzwingt keinen lokalen Rangverlust', det != 0)
details['local_sensitivity_example'] = dict(times=list(range(5)), parameters=['P0','A','f','k1','k2'],
    values=['1','1','1/3','ln(2)','ln(3)'],determinant=str(det), scope='Ein exakter innerer Beispielpunkt; keine generelle Identifizierbarkeit endlicher Daten bewiesen.')

# Exakte synthetische Kalibrierung mit harten Intervallen, keine Konfidenzintervalle.
N,D,E,e = map(F,['.004','.004','.008','.020'])
C=N*e/D
wide=(C*F('.4'),C)
ref=(C*F('.48')*F('.98'),C*F('.52')*F('1.02'))
assay=(E/F('.82'),E/F('.78'))
check('Breite Referenzgrenze exakt', wide == (F('.008'),F('.020')))
check('Enge Referenzgrenze exakt', ref == (F('.009408'),F('.010608')))
check('Selektiver Assay Kehrwertreihenfolge korrekt', assay == (F(2,205),F(2,195)))
check('Referenz verengt um Faktor zehn', (wide[1]-wide[0])/(ref[1]-ref[0]) == 10)
check('Assay verengt um Faktor 23.985', (wide[1]-wide[0])/(assay[1]-assay[0]) == F('23.985'))
check('Schnitt ist assaybegrenzt', ref[0] <= assay[0] <= assay[1] <= ref[1])
vsignal=((E-F('.0002'))/F('.82'),(E+F('.0002'))/F('.78'))
asignal=((N-F('.0002'))/vsignal[1],(N+F('.0002'))/vsignal[0])
saved=json.loads((R/'schlussfassung28/kalibrierpruefung.json').read_text())
check('Harte Signalfehler native Leistung nachgerechnet',
      np.allclose(list(map(float,vsignal)),saved['signal_error_case']['v0'],rtol=0,atol=1e-15))
check('Harte Signalfehler nativen Zugang nachgerechnet',
      np.allclose(list(map(float,asignal)),saved['signal_error_case']['aN'],rtol=0,atol=1e-15))
check('Signalfehler verbreitern beide Intervallenden', vsignal[0]<assay[0] and vsignal[1]>assay[1])
def intersect(a,b):
    if max(a[0],b[0]) > min(a[1],b[1]):
        raise ValueError('Leere Schnittmenge')
    return max(a[0],b[0]), min(a[1],b[1])
expect_rejected('Unvereinbare Kalibriergrenzen werden abgelehnt',lambda:intersect((F(0),F(1)),(F(2),F(3))))

# Vorhersage gegen einen künftigen beobachteten Wert versus den latenten Wert.
u,v,w,bb,dd,tt2,ssm2,ssj2,cm0,cr0=sp.symbols('u v w B delta tau2 sm2 sj2 cm cr',real=True)
expansion=sp.Poly(sp.expand((bb+u-v)**2-(w-dd-v)**2),u,v,w)
moments={(0,0,0):1,(1,0,0):0,(0,1,0):0,(0,0,1):0,
         (2,0,0):tt2,(0,2,0):ssj2,(0,0,2):ssm2,(1,1,0):cm0,(0,1,1):cr0}
expected=sum(coef*moments[power] for power,coef in expansion.terms())
check('Allgemeine Kovarianzformel aus Momenten hergeleitet',
      sp.simplify(expected-(bb**2+tt2-dd**2-ssm2-2*cm0+2*cr0))==0)
B,tau,sm,sj = map(F,['.003','.005','.002','.004'])
risk_cases=[]
for delta in [F('.020'),F('.001')]:
    vals=[(B+a*tau,q*sm-delta,c*sj) for a,q,c in product([-1,1],repeat=3)]
    mod=sum((m-z)**2 for m,r,z in vals)/len(vals)
    ref_r=sum((r-z)**2 for m,r,z in vals)/len(vals)
    check('Unabhängige Risikodifferenz delta='+str(delta),mod-ref_r == B*B+tau*tau-delta*delta-sm*sm)
    risk_cases.append(dict(delta=str(delta),correlated=False,model_observed=str(mod),constant_observed=str(ref_r)))
delta=F('.001')
vals=[(B+a*tau,q*sm-delta,a*sj) for a,q in product([-1,1],repeat=2)]
rm=sum((m-z)**2 for m,r,z in vals)/4
rr=sum((r-z)**2 for m,r,z in vals)/4
lm=sum(m*m for m,r,z in vals)/4
lr=sum(r*r for m,r,z in vals)/4
cm=sum((m-B)*z for m,r,z in vals)/4
cr=sum((r+delta)*z for m,r,z in vals)/4
check('Korrelierte Risikodifferenz exakt',rm-rr == B*B+tau*tau-delta*delta-sm*sm-2*cm+2*cr)
check('Korrelierter beobachteter Modellfehler exakt', rm==F('.000010'))
check('Korrelierter beobachteter Konstantenfehler exakt', rr==F('.000021'))
check('Latenter Modellfehler exakt', lm==F('.000034'))
check('Latenter Konstantenfehler exakt', lr==F('.000005'))
check('Beobachteter und latenter Vorteil können entgegengesetzt sein',rm<rr and lm>lr)
risk_cases.append(dict(delta=str(delta),correlated=True,model_observed=str(rm),constant_observed=str(rr),
    model_latent=str(lm),constant_latent=str(lr),cov_model=str(cm),cov_constant=str(cr)))
details['exact_synthetic_risk_cases']=risk_cases

# Vertrag für Quelle, Einheiten, Modellparameter und zeitlichen Prüfplan.
data_path=R/'zertifizierte_fehlergrenzen/data/nominal_data.csv'
raw=data_path.read_bytes()
data=list(csv.DictReader(data_path.open()))
bench=json.loads((R/'schlussfassung28/modellpruefung.json').read_text())
def verify_source(content):
    if hashlib.sha256(content).hexdigest()!=bench['source_data_sha256']:
        raise ValueError('Quelldaten verändert')
verify_source(raw)
check('Nominale Quelldaten unverändert',True)
expect_rejected('Negativkontrolle veränderte Quelldaten',lambda:verify_source(raw+b'\n'))

def validate_contract(item, rr):
    cid,kind=item['curve'],item['model']
    tt=np.array([float(x['time']) for x in rr])
    if not np.isfinite(tt).all() or not np.all(np.diff(tt)>0):
        raise ValueError('Zeitachse')
    sh=cid.startswith('Shintani')
    expected_time='h' if sh else 'min'
    if {x['time_unit'] for x in rr}!={expected_time}:
        raise ValueError('Zeiteinheit')
    use=list(range(len(rr))) if sh else list(range(1,len(rr)))
    if item['training_indices']!=use[:-2] or item['test_indices']!=use[-2:]:
        raise ValueError('Zeitaufteilung')
    for task in ['early','all_points']:
        n= {'E':3,'W':4,'P2':5}[kind] - (1 if task=='early' and not sh else 0)
        pp=item[task]['parameters'];bb=item[task]['bounds']
        if len(pp)!=n or len(bb)!=n:
            raise ValueError('Parameteranzahl')
        if not all(math.isfinite(x) and lo-1e-9<=x<=hi+1e-9 for x,(lo,hi) in zip(pp,bb)):
            raise ValueError('Parametergrenze')
    return True
for item in bench['details']:
    rr0=[x for x in data if x['curve']==item['curve']]
    check(item['curve']+' '+item['model']+' unabhängiger Eingabevertrag',validate_contract(item,rr0))

first=bench['details'][0]
rr0=[x for x in data if x['curve']==first['curve']]
bad=copy.deepcopy(first);bad['training_indices']=bad['training_indices']+bad['test_indices'][:1]
expect_rejected('Negativkontrolle zeitliches Datenleck',lambda:validate_contract(bad,rr0))
bad=copy.deepcopy(first);bad['early']['parameters'].append(0)
expect_rejected('Negativkontrolle zusätzliche Parameter',lambda:validate_contract(bad,rr0))
bad=copy.deepcopy(first);bad['early']['parameters'][0]=-1
expect_rejected('Negativkontrolle verletzte Parametergrenze',lambda:validate_contract(bad,rr0))
wrong=copy.deepcopy(rr0)
for x in wrong:x['time_unit']='min'
expect_rejected('Negativkontrolle verwechselte Zeiteinheit',lambda:validate_contract(first,wrong))

# Andere Optimierungsstruktur: Für festes k ist E linear in P0 und Amplitude.
# Eine LP minimiert den maximalen absoluten Fehler. Die eindimensionale Suche
# prüft die gespeicherten E-Kandidaten ohne differential_evolution oder SLSQP.
lp_results=[]
for cid in ['Shintani_A','Shintani_F','H_STD','H_HIGH']:
    rr0=[x for x in data if x['curve']==cid]
    tt=np.array([float(x['time']) for x in rr0]);yy=np.array([float(x['value']) for x in rr0])
    item=next(x for x in bench['details'] if x['curve']==cid and x['model']=='E')
    bb=item['all_points']['bounds'];stored=item['all_points']['parameters']
    def lp(logk, want_parameters=False):
        u=-np.expm1(-math.exp(logk)*tt)
        mat=np.vstack([np.column_stack([np.ones(len(tt)),u,-np.ones(len(tt))]),
                       np.column_stack([-np.ones(len(tt)),-u,-np.ones(len(tt))])])
        res=linprog([0,0,1],A_ub=mat,b_ub=np.r_[yy,-yy],bounds=[bb[0],bb[1],(0,None)],method='highs')
        if not res.success:raise ValueError('LP nicht gelöst: '+res.message)
        params=[res.x[0],res.x[1],logk]
        residual=max(abs(params[0]+params[1]*u-yy))
        return (float(residual),list(map(float,params))) if want_parameters else float(residual)
    grid=sorted(set(np.linspace(*bb[2],181).tolist()+[stored[2]]))
    scores=[lp(z) for z in grid]
    candidates=[(scores[i],grid[i]) for i in range(len(grid))]
    for i in range(1,len(grid)-1):
        if scores[i]<=scores[i-1] and scores[i]<=scores[i+1]:
            opt=minimize_scalar(lp,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-12})
            candidates.append((lp(opt.x),float(opt.x)))
    score,lk=min(candidates)
    value,params=lp(lk,True)
    original=next(x for x in bench['results'] if x['curve']==cid and x['model']=='E')['full_numerical_minimax']
    tol=max(2e-7,abs(original)*2e-5)
    check(cid+' E LP Kandidat innerhalb vorab gesetzter Toleranz',abs(value-original)<=tol)
    check(cid+' E LP Kandidat zulässig',all(lo-1e-9<=x<=hi+1e-9 for x,(lo,hi) in zip(params,bb)))
    lp_results.append(dict(curve=cid,stored_error=original,lp_error=value,difference=value-original,tolerance=tol,parameters=params,
        rate_grid_points=len(grid),scope='Zulässige numerische Kandidaten; kein globaler Optimalitätsnachweis.'))
details['independent_lp_search']=lp_results

report=dict(version=29,status='PASS',checks=checks,check_count=len(checks),details=details,
    scope='Interne Gegenrechnung mit anderer LP-Struktur und negativen Kontrollen; keine personell unabhängige Prüfung oder neue Messung.')
(S/'wissenschaft_pruefung.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status='PASS',checks=len(checks),lp=lp_results),ensure_ascii=False))
