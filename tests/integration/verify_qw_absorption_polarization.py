#!/usr/bin/env python3
"""Verify qualitative TE/TM ordering for the simple QW absorption benchmark.

Usage:
    verify_qw_absorption_polarization.py <absorption_TE.dat> <absorption_TM.dat>
"""

import sys

import numpy as np


def summarize(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    energy = data[:, 0]
    alpha = data[:, 1]
    peak_idx = int(np.argmax(alpha))
    strong = np.where(alpha > 0.05 * alpha.max())[0]
    onset = float(energy[strong[0]]) if strong.size else float("nan")
    integral = float(np.trapezoid(alpha, energy))
    return {
        "peak_energy": float(energy[peak_idx]),
        "peak_alpha": float(alpha[peak_idx]),
        "onset": onset,
        "integral": integral,
    }


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <absorption_TE.dat> <absorption_TM.dat>")
        sys.exit(1)

    te = summarize(sys.argv[1])
    tm = summarize(sys.argv[2])

    print(f"TE summary: {te}")
    print(f"TM summary: {tm}")

    if not te["onset"] < tm["onset"]:
        print("FAIL: TE onset is not earlier than TM onset")
        sys.exit(1)

    if not te["integral"] > tm["integral"]:
        print("FAIL: TE integrated absorption is not larger than TM")
        sys.exit(1)

    if not te["peak_energy"] < tm["peak_energy"]:
        print("FAIL: TE peak does not occur below the TM peak energy")
        sys.exit(1)

    print("PASS: QW absorption preserves the expected qualitative TE/TM ordering")


if __name__ == "__main__":
    main()
