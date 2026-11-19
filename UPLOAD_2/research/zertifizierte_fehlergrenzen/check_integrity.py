"""Meaningful negative checks: bad data/witnesses must not pass certification."""
from pathlib import Path
import json, shutil, tempfile
from verify_certificates import verify, exp_negative_interval, decimal_out, Q, HERE

def rejected_change(mutator):
    with tempfile.TemporaryDirectory() as temp:
        dest = Path(temp) / 'case'
        shutil.copytree(HERE / 'data', dest / 'data')
        shutil.copytree(HERE / 'results', dest / 'results')
        mutator(dest)
        try:
            verify(dest)
        except (AssertionError, ValueError):
            return True
        return False

def edit_candidate(dest, key, value):
    p = dest / 'results/candidates.json'
    j = json.loads(p.read_text())
    j['candidates'][0][key] = value
    p.write_text(json.dumps(j))

def run():
    checks = {}
    checks['valid_nominal_certificates'] = verify()[0]['status'] == 'VERIFIED'
    checks['negative_rate_rejected'] = rejected_change(lambda p: edit_candidate(p, 'rate', '-1'))
    checks['wrong_offset_rejected'] = rejected_change(lambda p: edit_candidate(p, 'p0', '0.2'))
    checks['wrong_positive_rate_rejected'] = rejected_change(lambda p: edit_candidate(p, 'rate', '0.01'))
    def damage_data(dest):
        p = dest / 'data/nominal_data.csv'
        p.write_bytes(p.read_bytes().replace(b'0.13957002373923574', b'0.14957002373923574'))
    checks['changed_data_rejected'] = rejected_change(damage_data)
    def damage_source(dest):
        p = dest / 'data/original/Shintani_A.csv'
        p.write_bytes(p.read_bytes()+b'\n')
    checks['changed_original_rejected'] = rejected_change(damage_source)
    checks['zero_exponent_exact'] = exp_negative_interval(Q(0)) == (Q(1), Q(1), 0)
    # 1/e must be below 3/8 and above 1/3; deliberately coarse independent bounds.
    lo, hi, _ = exp_negative_interval(Q(1))
    checks['one_over_e_enclosure_sanity'] = Q(1,3) < lo <= hi < Q(3,8)
    checks['outward_negative_rounding'] = decimal_out(Q(-1,3), False, 3) == '-0.334' and decimal_out(Q(-1,3), True, 3) == '-0.333'
    if not all(checks.values()):
        raise RuntimeError(checks)
    (HERE / 'results/integrity_checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps({'checks_passed': len(checks), 'checks_total': len(checks)}))

if __name__ == '__main__':
    run()
