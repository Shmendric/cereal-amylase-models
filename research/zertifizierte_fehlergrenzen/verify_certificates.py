"""Independent certificate verification using only the Python standard library.

No optimization, SciPy, NumPy or floating point operations enter the bounds.
Decimal CSV entries and candidate parameters are interpreted as exact rationals.
The mathematical justification is in PROOF.md. A successful verification proves
global minimax brackets for these nominal data, not historical scientific novelty
and not experimental adequacy of the model.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import argparse, csv, hashlib, json

HERE = Path(__file__).resolve().parent
if not __debug__:
    raise RuntimeError('Run without -O or -OO so that every certificate assertion is active.')

def verify_source_mapping(root_path, rows):
    expected = {}
    for name, time_field in [('Shintani_A', 'time_digitized_h'), ('Shintani_F', 'time_h')]:
        with (root_path / f'data/original/{name}.csv').open(encoding='utf-8-sig', newline='') as f:
            original = list(csv.DictReader(f))
        for row in original:
            expected[(name, row['point_id'])] = (row[time_field], row['maltose_g_per_g_input_flour'])
    with (root_path / 'data/original/Hackenberg.csv').open(encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            name = row['curve_id']
            expected[(name, name+'_'+row['time_min'])] = (row['time_min'], row['maltose_mg_per_g_dry_flour'])
    actual = {(row['curve'], row['point_id']): (row['time'], row['value']) for row in rows}
    assert len(actual) == len(rows), 'Duplicate nominal point'
    assert expected == actual, 'Mapping to original technical CSVs differs'

def exp_negative_interval(z, tolerance=Q(1, 10**42)):
    """For rational z>=0, enclose exp(-z) by reciprocal positive Taylor bounds."""
    if z < 0:
        raise ValueError('z must be nonnegative')
    if z == 0:
        return Q(1), Q(1), 0
    total, term, n = Q(1), Q(1), 0
    while True:
        next_term = term * z / (n + 1)
        ratio = z / (n + 2)
        if ratio < 1:
            upper_exp = total + next_term / (1 - ratio)
            lower, upper = 1 / upper_exp, 1 / total
            if upper - lower <= tolerance:
                return lower, upper, n
        n += 1
        total += next_term
        term = next_term
        if n > 10000:
            raise RuntimeError('Taylor enclosure did not converge')

def decimal_out(q, upwards=False, places=25):
    scaled = q * 10**places
    integer = -((-scaled.numerator) // scaled.denominator) if upwards else scaled.numerator // scaled.denominator
    sign = '-' if integer < 0 else ''
    magnitude = abs(integer)
    return f'{sign}{magnitude // 10**places}.{magnitude % 10**places:0{places}d}'

def verify(root_path=HERE):
    data_path = root_path / 'data/nominal_data.csv'
    data_bytes = data_path.read_bytes()
    data_hash = hashlib.sha256(data_bytes).hexdigest()
    provenance = json.loads((root_path / 'data/provenance.json').read_text(encoding='utf-8'))
    witness = json.loads((root_path / 'results/candidates.json').read_text(encoding='utf-8'))
    assert data_hash == provenance['nominal_data_sha256'] == witness['data_sha256'], 'Nominal data hash mismatch'
    for source in provenance['source_files']:
        assert hashlib.sha256((root_path / source['local_path']).read_bytes()).hexdigest() == source['sha256'], 'Source hash mismatch'
    with data_path.open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    verify_source_mapping(root_path, rows)
    expected = {'Shintani_A': 8, 'Shintani_F': 8, 'H_STD': 7, 'H_HIGH': 7}
    assert set(row['curve'] for row in rows) == set(expected)
    assert [c['curve'] for c in witness['candidates']] == list(expected)
    summaries = []
    point_records = []
    for c in witness['candidates']:
        rr = [row for row in rows if row['curve'] == c['curve']]
        assert len(rr) == c['n_points'] == expected[c['curve']], 'Point count mismatch'
        t = [Q(row['time']) for row in rr]
        y = [Q(row['value']) for row in rr]
        assert t[0] >= 0 and all(t[i] < t[i+1] for i in range(len(t)-1))
        assert all(v >= 0 for v in y)
        p0, amp, d = Q(c['p0']), Q(c['amplitude']), Q(c['rate'])
        assert p0 >= 0 and amp > 0 and d > 0, 'Inadmissible witness'
        residuals = []
        for i, (ti, yi) in enumerate(zip(t, y)):
            e_lo, e_hi, n = exp_negative_interval(d * ti)
            # amp>0; f in [p0+amp*(1-e_hi), p0+amp*(1-e_lo)]
            r_lo, r_hi = yi - p0 - amp * (1-e_lo), yi - p0 - amp * (1-e_hi)
            assert r_lo <= r_hi
            residuals.append((r_lo, r_hi))
            point_records.append(dict(curve=c['curve'], point_id=rr[i]['point_id'], time=rr[i]['time'], value=rr[i]['value'], residual_lower=decimal_out(r_lo), residual_upper=decimal_out(r_hi, True), taylor_degree=n))
        M = max(max(abs(a), abs(b)) for a, b in residuals)
        alternatives = []
        for indices in combinations(range(len(t)), 4):
            bounds, signs = [], []
            for i in indices:
                a, b = residuals[i]
                if a > 0:
                    bounds.append(a); signs.append(1)
                elif b < 0:
                    bounds.append(-b); signs.append(-1)
                else:
                    break
            if len(signs) == 4 and all(signs[j] != signs[j+1] for j in range(3)):
                alternatives.append((min(bounds), indices, signs))
        assert alternatives, 'No certified four-point alternation'
        m, indices, signs = max(alternatives, key=lambda item: item[0])
        assert m > 0 and M >= m and M-m < Q(1, 10**10), 'Certificate gap exceeds 1e-10'
        max_drop, drop_indices = max((max(Q(0), y[i]-y[j]), (i, j)) for i in range(len(t)) for j in range(i+1, len(t)))
        mono = max_drop / 2
        # Construct the exact monotone minimax sampled sequence as a second witness.
        monotone_values = [max(Q(0), max(y[:i+1])-mono) for i in range(len(y))]
        assert all(monotone_values[i] <= monotone_values[i+1] for i in range(len(y)-1))
        assert max(abs(a-b) for a, b in zip(y, monotone_values)) <= mono
        assert mono <= M
        summary = dict(curve=c['curve'], n_points=len(t), lower=decimal_out(m), upper=decimal_out(M, True), certificate_gap_upper=decimal_out(M-m, True), monotone_floor_lower=decimal_out(mono), monotone_floor_upper=decimal_out(mono, True), monotone_to_exponential_ratio_lower=decimal_out(mono/M), monotone_to_exponential_ratio_upper=decimal_out(mono/m, True), active_indices_zero_based=list(indices), active_point_ids=[rr[i]['point_id'] for i in indices], active_times=[rr[i]['time'] for i in indices], signs=signs, monotonicity_drop_pair=[rr[i]['point_id'] for i in drop_indices], p0=c['p0'], amplitude=c['amplitude'], rate=c['rate'], time_unit=c['time_unit'], value_unit=c['value_unit'])
        summary['maximum_absolute_monotone_gain_upper'] = decimal_out(M-mono, True)
        summary['maximum_relative_monotone_gain_upper'] = decimal_out(1-mono/M, True)
        summaries.append(summary)
    return dict(status='VERIFIED', arithmetic='fractions.Fraction; outward rounding to 25 decimal places', certificate_gap_tolerance='0.0000000001', input_sha256=data_hash, results=summaries), point_records

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=HERE)
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args()
    result, points = verify(args.root)
    if not args.check_only:
        out = args.root / 'results'
        (out / 'certificates.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
        for filename, records in [('certificates.csv', result['results']), ('point_enclosures.csv', points)]:
            with (out / filename).open('w', encoding='utf-8', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
    print(json.dumps({'status': result['status'], 'curves': len(result['results']), 'max_gap': max(r['certificate_gap_upper'] for r in result['results'])}))

if __name__ == '__main__':
    main()
