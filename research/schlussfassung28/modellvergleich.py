"""Bounded nominal E/W/P2 comparison. Numerical candidates, no new global certificates."""
from pathlib import Path
import csv
import hashlib
import json
import time
import numpy as np
from scipy.optimize import differential_evolution, least_squares, minimize

S = Path(__file__).resolve().parent
R = S.parent
DATA = R / 'zertifizierte_fehlergrenzen/data/nominal_data.csv'
SEED = 20261112
rows = list(csv.DictReader(DATA.open()))
checks = {}

def check(k, v):
    checks[k] = bool(v)
    if not v:
        raise AssertionError(k)

def pred(p, t, kind, fixed):
    if fixed is None:
        p0, delta, *q = p
    else:
        p0, (delta, *q) = fixed, p
    if kind == 'P2':
        f, log1, log2 = q
        return p0 + delta * (f * (-np.expm1(-np.exp(log1) * t)) +
                             (1 - f) * (-np.expm1(-np.exp(log2) * t)))
    logk = q[0]
    n = q[1] if kind == 'W' else 1.
    return p0 + delta * (-np.expm1(-(np.exp(logk) * t) ** n))

def bounds(kind, yy, cfg, fixed):
    cap, kmin, kmax, nmax = cfg
    b = ([] if fixed is not None else [(0., float(min(yy)))]) + [(1e-5 if fixed is None else 1e-6, cap)]
    if kind == 'P2':
        return b + [(0., 1.), (np.log(kmin), np.log(kmax)), (np.log(kmin), np.log(kmax))]
    return b + [(np.log(kmin), np.log(kmax))] + ([(.1, nmax)] if kind == 'W' else [])

def p2_from_e(p, fixed):
    return [*p[:-1], .5, p[-1], p[-1]]

def fit(kind, t, y, cfg, fixed, objective, embeddings):
    b = bounds(kind, y, cfg, fixed)
    # Normalize residual scale only; original units are retained in reporting.
    scale = max(float(np.ptp(y)), float(np.max(np.abs(y))) * .1, 1e-8)
    residual = lambda p: (pred(p, t, kind, fixed) - y) / scale
    fun = lambda p: float(residual(p) @ residual(p)) if objective == 'least_squares' else float(np.max(abs(residual(p))))
    candidates = []
    runs = []
    lo, hi = np.array(b).T
    for seed in [SEED, SEED + 1]:
        z = differential_evolution(fun, b, seed=seed, popsize=10, maxiter=500,
                                   tol=1e-9, atol=1e-12, polish=False)
        candidates.append(np.asarray(z.x));runs.append(dict(seed=seed, objective=float(z.fun), success=bool(z.success), iterations=int(z.nit)))
    candidates += [np.clip(np.asarray(p), lo, hi) for p in embeddings]
    polished = []
    for p in candidates:
        if objective == 'least_squares':
            z = least_squares(residual, np.clip(p, lo + 1e-12, hi - 1e-12),
                              bounds=(lo, hi), x_scale='jac', max_nfev=5000,
                              ftol=1e-12, xtol=1e-12, gtol=1e-12)
            polished.append(np.asarray(z.x))
        else:
            # Epigraph polishing for the maximum absolute residual.
            x0 = np.r_[p, max(abs(residual(p))) + 1e-10]
            con = lambda x: np.r_[x[-1] - residual(x[:-1]), x[-1] + residual(x[:-1])]
            z = minimize(lambda x: x[-1], x0, method='SLSQP',
                         bounds=[*b, (0, None)], constraints={'type': 'ineq', 'fun': con},
                         options={'maxiter': 2000, 'ftol': 1e-12})
            if np.all(np.isfinite(z.x)) and np.min(con(z.x)) >= -1e-7:
                polished.append(np.clip(z.x[:-1], lo, hi))
    # Embeddings remain admissible even if numerical polishing makes them worse.
    all_candidates = candidates + polished
    best = min(all_candidates, key=fun)
    return best, dict(parameters=best.tolist(), bounds=b, objective=objective,
                      objective_normalized=fun(best), global_search_runs=runs,
                      finite_numerically_usable=bool(np.isfinite(fun(best))),
                      active_bounds=[i for i, x in enumerate(best) if abs(x-lo[i]) < 1e-6 or abs(x-hi[i]) < 1e-6])

results = []
records = []
start = time.monotonic()
for cid in ['Shintani_A', 'Shintani_F', 'H_STD', 'H_HIGH']:
    rr = [r for r in rows if r['curve'] == cid]
    t = np.array([float(r['time']) for r in rr]); y = np.array([float(r['value']) for r in rr])
    is_shintani = cid.startswith('Shintani')
    fixed = None if is_shintani else float(y[0])
    use = np.arange(len(t)) if is_shintani else np.arange(1, len(t))
    train = use[:-2]; test = use[-2:]
    cfg = (.8, .001, 30., 3.) if is_shintani else (100., 1e-6, 10., 4.)
    early_solutions = {}; full_solutions = {}; full_reports = {}
    for kind in ['E', 'W', 'P2']:
        embedding_early = [] if kind == 'E' else ([*early_solutions['E'], 1.] if kind == 'W' else p2_from_e(early_solutions['E'], fixed))
        pe, er = fit(kind, t[train], y[train], cfg, fixed, 'least_squares', [] if kind=='E' else [embedding_early])
        early_solutions[kind] = pe
        # For minimax only, all points and a freely fitted nonnegative P0 are used for every curve.
        # This is a different task from fixed-initial-value continuation.
        embedding_all = [] if kind=='E' else ([*full_solutions['E'], 1.] if kind=='W' else p2_from_e(full_solutions['E'], None))
        pa, ar = fit(kind, t, y, cfg, None, 'minimax', [] if kind=='E' else [embedding_all])
        full_solutions[kind] = pa; full_reports[kind] = ar
        ye = pred(pe, t, kind, fixed); ya = pred(pa, t, kind, None)
        result = dict(curve=cid, model=kind, time_unit=rr[0]['time_unit'], value_unit=rr[0]['value_unit'],
                      training_points=len(train), later_points=2,
                      early_RMSE=float(np.sqrt(np.mean((ye[train]-y[train])**2))),
                      later_RMSE=float(np.sqrt(np.mean((ye[test]-y[test])**2))),
                      full_numerical_minimax=float(np.max(abs(ya-y))),
                      amplitude_cap=cfg[0], initial_value_policy='free_nonnegative' if fixed is None else 'reported_initial_fixed_for_continuation',
                      training_parameters=len(pe), full_parameters=len(pa),
                      training_residual_degrees_of_freedom=len(train)-len(pe))
        results.append(result)
        records.append(dict(curve=cid,model=kind,early=er,all_points=ar,training_indices=train.tolist(),test_indices=test.tolist(),
                            early_predictions=ye.tolist(),all_points_predictions=ya.tolist()))
        check(cid+'_'+kind+'_finite', all(np.isfinite([result['early_RMSE'], result['later_RMSE'], result['full_numerical_minimax']])))
        check(cid+'_'+kind+'_nonnegative', min(ya) >= 0 and min(ye) >= 0)
        check(cid+'_'+kind+'_monotone', np.all(np.diff(ya) >= -1e-12))
    # A parallel-pool family contains E exactly; numerical results cannot be allowed to lose its candidate.
    ee, pp = [next(x for x in results if x['curve']==cid and x['model']==k) for k in ['E','P2']]
    check(cid+'_P2_contains_E_training', pp['early_RMSE'] <= ee['early_RMSE'] + 1e-8)
    check(cid+'_P2_contains_E_minimax', pp['full_numerical_minimax'] <= ee['full_numerical_minimax'] + 1e-8)
    print(json.dumps(dict(curve=cid,completed=True)),flush=True)

# Mechanistic nonuniqueness of a single exponential and the forecasting risk identity.
import sympy as sp
t,b,d,P0 = sp.symbols('t b d P0', positive=True)
PI=P0+b*(1-sp.exp(-d*t))/d
check('inactivation_integrated_rate', sp.simplify(sp.diff(PI,t)-b*sp.exp(-d*t))==0)
delta, bias, tau2, sigma2 = sp.symbols('delta bias tau2 sigma2', real=True)
check('forecast_risk_difference', sp.expand((bias*bias+tau2+sigma2)-(delta*delta+sp.Symbol('sigma_m2')+sigma2)) == bias*bias+tau2-delta*delta-sp.Symbol('sigma_m2'))

with (S/'modellvergleich.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
report=dict(status='PASS',version=28,seed=SEED,seconds=round(time.monotonic()-start,2),
            source_data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),
            scope='Nominal bounded comparison of E, W and parallel nonnegative pools. Numerical minimax candidates are upper bounds, not globally certified optima. No new material group, independent validation or parameter identification is claimed.',
            new_measurements=0,results=results,details=records,checks=checks,check_count=len(checks))
(S/'modellpruefung.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status='PASS',checks=len(checks),results=results),ensure_ascii=False),flush=True)
