#!/usr/bin/env python3
"""Generate all physics verification figures."""
import subprocess
import sys
from pathlib import Path

SCRIPTS = [
    'verify_qwz_chern.py',
    'verify_bhz_z2.py',
    'verify_landau_levels.py',
    'sweep_rashba_bdg.py',
]

def main():
    results = []
    for script in SCRIPTS:
        result = subprocess.run([sys.executable, f'scripts/{script}'], capture_output=True)
        status = 'PASS' if result.returncode == 0 else 'FAIL'
        results.append((script, status))
        print(f'{script}: {status}')

    # Summary
    passed = sum(1 for _, s in results if s == 'PASS')
    print(f'\nPassed: {passed}/{len(results)}')
    return 0 if passed == len(results) else 1

if __name__ == '__main__':
    sys.exit(main())