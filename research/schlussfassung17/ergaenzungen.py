"""Check retained sensitivity results and derive a physical three-point guard.

No new measured observations, parameter fits or experimental error estimates.
"""
from pathlib import Path
import csv, hashlib, itertools, json, math, platform, sys
import importlib.metadata as metadata
import numpy as np
from scipy.stats import qmc

R = Path(__file__).resolve().parent.parent
O = Path(__file__).resolve().parent
H = R / 'datenpruefung_2026_11_05/hackenberg_pilot'
S = R / 'robustheitspruefung_2026_11_05/robustness_runs.csv'

def read_csv(path):
    with path.open(encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n')

raw = read_csv(H / 'digitized_maltose.csv')
runs = read_csv(S)
sobol = qmc.Sobol(d=16, scramble=True, seed=20261105).random_base2(6)
times = np.array([5., 10., 20., 60., 120., 180.])
configurations = ['D100_N4', 'D50_N4', 'D100_N2']
source_errors = []
groups = {}
for row in runs:
    cid, mode = row['curve_id'], row['mode']
    j = int(row['scenario']); ic = ['H_STD', 'H_HIGH'].index(cid)
    rr = [r for r in raw if r['curve_id'] == cid]
    if mode == 'baseline_only':
        s = np.full(16, [0., .5, 1.][j])
        y = np.array([float(r['maltose_mg_per_g_dry_flour']) for r in rr[1:]])
    else:
        s = sobol[j]
        y0, y50 = 590 + 2*s[0]-1, 166 + 2*s[1]-1
        yp = np.array([float(r['pixel_y']) for r in rr[1:]]) + (2*s[2+6*ic:8+6*ic]-1)*2.
        y = 50*(y0-yp)/(y0-y50)
    p0 = float(rr[0]['maltose_mg_per_g_dry_flour'])
    if mode != 'pixels_only':
        p0 += (2*s[14+ic]-1)*float(rr[0]['reported_SEM_mg_per_g'])
    delta, k, shape = (float(row[n]) for n in ['delta', 'k_per_min', 'shape'])
    pred = p0 + delta*(-np.expm1(-(k*times)**shape))
    early = float(np.sqrt(np.mean((pred[:4]-y[:4])**2)))
    late = float(np.sqrt(np.mean((pred[4:]-y[4:])**2)))
    source_errors += [abs(p0-float(row['initial_product'])),
                      abs(early-float(row['early_RMSE'])),
                      abs(late-float(row['late_RMSE']))]
    key = mode, cid, j, row['configuration']
    groups.setdefault(key, {})[row['model']] = late

paired = []; summary = []; consistency = []
for mode, cid in itertools.product(['baseline_only', 'pixels_only', 'combined'], ['H_STD', 'H_HIGH']):
    count = 3 if mode == 'baseline_only' else 64
    tally = {'E_all': 0, 'W_all': 0, 'changes_with_bounds': 0}
    for j in range(count):
        diffs = []
        for conf in configurations:
            pair = groups[mode, cid, j, conf]
            diffs.append(pair['E'] - pair['W'])
        status = 'E_all' if all(d < 0 for d in diffs) else 'W_all' if all(d > 0 for d in diffs) else 'changes_with_bounds'
        tally[status] += 1
        paired.append({'mode': mode, 'curve_id': cid, 'scenario': j,
                       **dict(zip(configurations, diffs)), 'decision': status})
    summary.append({'mode': mode, 'curve_id': cid, 'scenario_count': count, **tally})

for row in runs:
    key = row['mode'], row['curve_id'], int(row['scenario']), row['configuration']
    pair = groups[key]
    consistency.append(abs(pair['E']-pair['W']-float(row['E_late_minus_W_late'])))

def inverse(y1, a, bdiff, t1, h):
    d = math.log1p((a-bdiff)/bdiff)/h
    scale = a/(-math.expm1(-d*h))
    p0 = y1-scale*math.expm1(d*t1)
    initial_rate = d*scale*math.exp(d*t1)
    plateau = y1+scale
    return p0, initial_rate, d, plateau, scale*math.expm1(d*t1)

cases = []; errors = []
for p0, rate, d, t1, h in itertools.product([0., .015], [.001, .004], [.001, .01, .04, .2], [0., 1., 5., 10.], [2., 5., 20., 40.]):
    y = [p0+rate*(-math.expm1(-d*t))/d for t in [t1, t1+h, t1+2*h]]
    aa, bb = y[1]-y[0], y[2]-y[1]
    back = inverse(y[0], aa, bb, t1, h)
    aitken = y[0]-(y[1]-y[0])**2/(y[2]-2*y[1]+y[0])
    expected_plateau = p0+rate/d
    errs = [abs(back[0]-p0), abs(back[1]-rate), abs(back[2]-d), abs(aitken-expected_plateau)]
    errors.extend(errs)
    cases.append({'P0': p0, 'b': rate, 'd': d, 't1': t1, 'h': h,
                  'Y1': y[0], 'Y2': y[1], 'Y3': y[2],
                  'reconstructed_P0': back[0], 'aitken_plateau': aitken,
                  'admissibility_threshold': back[4], 'max_absolute_error': max(errs)})

negative = inverse(.001, .010, .005, 5., 20.)
example = {'kind': 'synthetic_not_measured', 'unit': 'g maltose per g flour',
           'times_min': [5., 25., 45.], 'product_values': [.001, .011, .016],
           'A': .010, 'B': .005, 'reconstructed_P0': negative[0],
           'reconstructed_b': negative[1], 'reconstructed_d_per_min': negative[2],
           'plateau': negative[3], 'nonnegative_P0_requires_Y1_at_least': negative[4]}
checks = {
    'all_1572_stored_predictions_checked_without_refitting': len(runs) == 1572 and max(source_errors) < 1e-10,
    'all_786_model_pairs_complete': len(groups) == 786 and all(set(g) == {'E', 'W'} for g in groups.values()),
    'all_recorded_RMSE_differences_reconstructed': max(consistency) < 1e-10,
    'all_262_paired_scenarios_have_three_boundary_choices': len(paired) == 262 and len(summary) == 6,
    'all_counts_conserve_scenarios': all(r['E_all']+r['W_all']+r['changes_with_bounds'] == r['scenario_count'] for r in summary),
    '256_exact_model_cases_reconstruct_parameters_and_plateau': len(cases) == 256 and max(errors) < 2e-10,
    'increasing_positive_observations_can_imply_negative_P0': negative[0] < 0 and .010 > .005 > 0,
    'zero_time_reduces_nonnegativity_condition_to_Y1': inverse(.015, .010, .005, 0., 20.)[4] == 0.,
}
checks = {name: bool(value) for name, value in checks.items()}
write_csv(O/'rangvergleich_gepaart.csv', paired)
write_csv(O/'rangvergleich_zusammenfassung.csv', summary)
write_csv(O/'rekonstruktion_kontrollfaelle.csv', cases)
dump(O/'gegenbeispiel_anfangsmenge.json', example)
dump(O/'ergaenzungen_ergebnisse.json', {
    'target_manuscript_version': 17, 'source_data_version': 14,
    'scope': 'Deterministic reuse of historical sensitivity scenarios and synthetic algebra checks; no new measurements or novelty finding.',
    'historical_runs_sha256': hashlib.sha256(S.read_bytes()).hexdigest(),
    'historical_parameters_not_refitted': True,
    'max_reconstructed_RMSE_difference': max(source_errors),
    'max_algebra_control_absolute_error': max(errors),
    'paired_boundary_comparison': summary, 'synthetic_example': example,
    'checks': checks, 'all_pass': all(checks.values())})
packages = ['numpy', 'scipy', 'matplotlib', 'python-docx', 'Pillow', 'PyMuPDF', 'lxml']
dump(O/'laufzeitumgebung.json', {'python': platform.python_version(), 'executable': sys.executable,
    'platform': platform.platform(), 'packages': {p: metadata.version(p) for p in packages},
    'new_analysis_random_seed': None,
    'inherited_sobol_seed_for_reconstruction': 20261105,
    'historical_seeds': {'Shintani': 20261105, 'Hackenberg_Sobol': 20261105,
                         'synthetic_precision': 2026110513, 'independent_nominal_fits': 2026110515},
    'provenance_note': 'Environment observed during version 17 verification; not a timestamped execution log of every historical script.'})
print(json.dumps({'checks': checks, 'summary': summary, 'example': example,
                  'max_source_reconstruction_error': max(source_errors)}, ensure_ascii=False, indent=2))
assert all(checks.values()), [k for k, v in checks.items() if not v]
