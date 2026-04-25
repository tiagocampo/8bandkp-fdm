#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: verify_scattering_transition_sanity.py <scattering_rates.dat>", file=sys.stderr)
        return 2

    path = Path(sys.argv[1])
    if not path.exists():
        print(f"missing file: {path}", file=sys.stderr)
        return 1

    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 7:
        print("unexpected scattering_rates.dat format", file=sys.stderr)
        return 1

    dE_meV = data[:, 2]
    min_dE = float(np.min(dE_meV))

    # Intersubband scattering should not be reported between numerically
    # degenerate Kramers partners. A near-zero transition energy indicates
    # that spin partners were misclassified as distinct subbands.
    if min_dE < 1.0e-3:
        print(
            f"found near-zero intersubband transition energy: {min_dE:.6e} meV",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
