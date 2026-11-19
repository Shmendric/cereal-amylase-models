"""Run the frozen research edition without changing its scientific files."""
from pathlib import Path
import argparse
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-optimization', action='store_true', help='Recheck stored nominal parameters instead of re-optimising them.')
    parser.add_argument('--verify-only', action='store_true', help='Check the frozen research files without running all calculations.')
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs' / 'latest')
    args = parser.parse_args()
    if not __debug__:
        parser.error('Run Python without -O or -OO; the checks use assertions.')
    if args.verify_only:
        command = [sys.executable, str(ROOT / 'research' / 'verify_release.py')]
    else:
        command = [sys.executable, str(ROOT / 'research' / 'reproduce.py'), '--output', str(args.output.resolve())]
        if args.skip_optimization:
            command.append('--skip-optimization')
    return subprocess.call(command)

if __name__ == '__main__':
    raise SystemExit(main())
