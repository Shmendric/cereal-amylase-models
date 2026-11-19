"""Verify the files in the delivered repository snapshot using SHA-256."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent
EXCLUDED = {'MANIFEST_SHA256.json', 'SHA256SUMS.txt'}

def tracked(path):
    rel = path.relative_to(ROOT)
    return (path.is_file() and not path.is_symlink()
            and rel.as_posix() not in EXCLUDED
            and not any(part in {'.git', '__pycache__', '.venv', 'venv'} for part in rel.parts)
            and rel.parts[0] != 'outputs'
            and not rel.as_posix().startswith('research/ergebnisse/laufprotokolle/'))

def main():
    expected = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text(encoding='utf-8'))
    actual = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(ROOT.rglob('*')) if tracked(p)}
    missing = sorted(expected.keys() - actual.keys())
    added = sorted(actual.keys() - expected.keys())
    changed = sorted(k for k in expected.keys() & actual.keys() if expected[k] != actual[k])
    sums_ok = (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf-8') == ''.join(v + '  ' + k + '\n' for k, v in expected.items())
    passed = not (missing or added or changed) and sums_ok
    print(json.dumps({'status': 'PASS' if passed else 'FAIL', 'files': len(expected),
                      'missing': missing, 'added': added, 'changed': changed,
                      'text_manifest_matches': sums_ok}, indent=2))
    return 0 if passed else 1

if __name__ == '__main__':
    raise SystemExit(main())
