"""Exact synthetic calibration intervals. No experimental calibration data."""
from fractions import Fraction as F
from pathlib import Path
import csv
import json

S = Path(__file__).resolve().parent
N, D, E, e = map(F, ['0.004', '0.004', '0.008', '0.020'])
C = N * e / D
checks = {}

def check(name, value):
    checks[name] = bool(value)
    if not value:
        raise AssertionError(name)

def from_reference(rho, eta):
    return C * rho[0] * eta[0], C * rho[1] * eta[1]

def from_assay(kappa):
    return E / kappa[1], E / kappa[0]

def intersection(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    if hi < lo:
        raise ValueError('Incompatible calibration bounds')
    return lo, hi

def access(v):
    return N / v[1], N / v[0]

wide = from_reference(tuple(map(F, ['0.4', '1.0'])), (F(1), F(1)))
narrow = from_reference(tuple(map(F, ['0.48', '0.52'])), tuple(map(F, ['0.98', '1.02'])))
assay = from_assay(tuple(map(F, ['0.78', '0.82'])))
joint = intersection(narrow, assay)
point = from_assay((F('0.8'), F('0.8')))
check('reference_wide_exact', wide == (F('.008'), F('.020')))
check('reference_narrow_exact', narrow == (F('.009408'), F('.010608')))
check('assay_exact', assay == (F(2, 205), F(2, 195)))
check('joint_intersection_assay_limited', joint == assay)
check('known_native_assay_single_point', point == (F('.010'), F('.010')))
check('reference_width_reduction_tenfold', (wide[1] - wide[0]) / (narrow[1] - narrow[0]) == 10)
check('assay_width_reduction_exact', (wide[1] - wide[0]) / (assay[1] - assay[0]) == F(4797, 200))

rows = []
for name, v in [('eta_known_rho_wide', wide), ('reference_tight', narrow),
                ('native_assay_tight', assay), ('joint_tight', joint), ('native_assay_exact', point)]:
    a = access(v)
    rows.append(dict(case=name, v0_lower=float(v[0]), v0_upper=float(v[1]),
                     aN_lower=float(a[0]), aN_upper=float(a[1]),
                     width=float(v[1] - v[0]),
                     width_reduction_vs_wide=None if v[0] == v[1] else float((wide[1] - wide[0]) / (v[1] - v[0]))))
    # Construct every intermediate v with a common eta=1 where possible.
    if name in ['reference_tight', 'joint_tight']:
        for j in range(401):
            v0 = v[0] + (v[1] - v[0]) * F(j, 400)
            h = v0 / C
            eta_lo, eta_hi = max(F('.98'), h / F('.52')), min(F('1.02'), h / F('.48'))
            check(f'{name}_construction_{j}', eta_lo <= eta_hi)
            eta = (eta_lo + eta_hi) / 2
            rho = h / eta
            aN = N / v0
            aR = rho * aN
            kappa = E / v0
            check(f'{name}_signals_{j}', aR * eta * e == D and aN * v0 == N and kappa * v0 == E)

# Hard signal errors are kept distinct from calibration intervals.
dN, dE = F('.0002'), F('.0002')
measurement_v = ((E - dE) / F('.82'), (E + dE) / F('.78'))
measurement_a = ((N - dN) * F('.78') / (E + dE), (N + dN) * F('.82') / (E - dE))
check('signal_error_widens_intervals', measurement_v[0] < assay[0] and measurement_v[1] > assay[1])
check('synthetic_truth_reconstructed', N / point[0] == F('.4') and D / (F(1) * e) == F('.2'))

with (S / 'kalibrierintervalle.csv').open('w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
report = dict(status='PASS', kind='synthetic_exact_interval_calculation', new_measurements=0,
              inputs=dict(N=float(N), D=float(D), E=float(E), e=float(e)),
              cases=rows, signal_error_case=dict(v0=list(map(float, measurement_v)), aN=list(map(float, measurement_a))),
              check_count=len(checks), checks=checks,
              scope='Fixed positive signals for interval contraction; hard signal errors examined separately. Intervals are assumptions, not confidence intervals.')
(S / 'kalibrierpruefung.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
print(json.dumps(dict(status=report['status'], check_count=len(checks), cases=rows), ensure_ascii=False))
