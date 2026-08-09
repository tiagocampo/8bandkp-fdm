"""TDD-red/green pin for P2: `lecture_13_topological.py` degeneracy detection.

Pins the four failure modes of the slim Pfaffian `bcrit_pfaffian` extraction
in `scripts/lecture_13_topological.py:section_wire_rung()`:

  1. file absent               -> bcrit_pfaffian is None  (file not found)
  2. file present, empty match -> bcrit_pfaffian is None  (file empty/unmatched)
  3. file saturated at floor   -> bcrit_pfaffian is None  (T6: max > 0 but
                                                          (max-min)/max < tol)
  4. file all-zero entries     -> bcrit_pfaffian is None  (P2: every line is 0;
                                                          pmax == 0 short-
                                                          circuit fix)

Test (4) is the P2 regression: pre-fix, when every |Pf| entry is 0, the
guard `pmax > 0 and (pmax-pmin)/pmax < tol` short-circuits to False, and
the else branch returns `min(pf_mags, key=lambda B: pf_mags[B])` — which
returns the **first** B key in dict insertion order, not a phase boundary.
Post-fix, an explicit `elif pmax == 0` branch sets `bcrit_pfaffian = None`
with the same WARN as the saturated case.

Strategy: import `lecture_13_topological.py`, monkey-patch its `REPO` to a
tmpdir containing a synthetic `wire_slim_pfaffian_witness.dat`, and stub
`_run_verifier` so `section_wire_rung()` returns from the verifier calls
cleanly and reaches the Pfaffian block. We assert `bcrit_pfaffian is None`
for the degenerate (cases 3, 4) and absent (case 1) regimes, and that
`bcrit_pfaffian` is a float for the varying (case 5) regime (added as a
regression guard against over-correction).
"""
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

# Make scripts/ importable
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import lecture_13_topological as L13  # noqa: E402


def _make_fake_verifier(stdout):
    """Return a _run_verifier stub that yields `stdout` and rc=0."""
    def _stub(_verifier_name, _exe, _config):
        return 0, stdout
    return _stub


_FAKE_1D_STDOUT = (
    "minigap(meV): Bx=0.0:2.907 Bx=2.8:0.019 Bx=5.0:3.799\n"
    "By=0.0 Bx=0.0 mu=0.6601\n"
)


class PfaffianDegeneracyDetection(unittest.TestCase):
    """Pin the four failure modes + one regression-vary case."""

    def _run_with_pf(self, pf_text):
        """Invoke L13.section_wire_rung against a tmpdir with `pf_text`
        written to output/wire_slim_pfaffian_witness.dat. Mocks the
        `_run_verifier` so the function reaches the Pfaffian block.
        Returns (bcrit_curve, bcrit_2d, bcrit_pfaffian).
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            (tmp_path / "output").mkdir(parents=True, exist_ok=True)
            if pf_text is not None:
                (tmp_path / "output" / "wire_slim_pfaffian_witness.dat").write_text(pf_text)

            with mock.patch.object(L13, "REPO", tmp_path), \
                 mock.patch.object(L13, "_run_verifier",
                                   _make_fake_verifier(_FAKE_1D_STDOUT)):
                buf = io.StringIO()
                with redirect_stdout(buf):
                    _, (bcrit_curve, bcrit_2d, bcrit_pfaffian) = L13.section_wire_rung("/dev/null/exe")
        return bcrit_curve, bcrit_2d, bcrit_pfaffian

    # --- Case 1: file absent ---
    def test_pfaffian_absent_returns_none(self):
        _, _, bpf = self._run_with_pf(pf_text=None)
        self.assertIsNone(bpf, "Pfaffian file absent must yield bcrit_pfaffian=None")

    # --- Case 2: file present but no matched lines ---
    def test_pfaffian_empty_match_returns_none(self):
        _, _, bpf = self._run_with_pf(pf_text="# only comments, no body\n")
        self.assertIsNone(bpf, "Unmatched Pfaffian file must yield bcrit_pfaffian=None")

    # --- Case 3: file saturated at FEST precision floor (T6) ---
    def test_pfaffian_saturated_floor_returns_none(self):
        # All values equal: pmax > 0 but (pmax-pmin)/pmax = 0 < tol
        pf_text = "\n".join(
            f"B={B:.6E} |Pf|=4.000E-08" for B in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        ) + "\n"
        _, _, bpf = self._run_with_pf(pf_text=pf_text)
        self.assertIsNone(
            bpf,
            "Pfaffian file saturated at floor must yield bcrit_pfaffian=None "
            "(T6 detection)"
        )

    # --- Case 4: file all-zero (P2 fix) ---
    def test_pfaffian_all_zero_returns_none(self):
        # Closure regime: every entry's |Pf|=0 (s2_sign=0 across the grid).
        pf_text = "\n".join(
            f"B={B:.6E} |Pf|=0.000E+00" for B in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        ) + "\n"
        _, _, bpf = self._run_with_pf(pf_text=pf_text)
        # PRE-FIX: returns the first B key (0.0), not None. POST-FIX: None.
        self.assertIsNone(
            bpf,
            "Pfaffian file with all-zero entries must yield bcrit_pfaffian=None "
            "(P2 fix; pre-fix short-circuit returned the first B key)"
        )

    # --- Case 5: regression guard — varying profile must still produce a float ---
    def test_pfaffian_varying_returns_minimum(self):
        pf_text = "\n".join(
            f"B={B:.6E} |Pf|={pf:.6E}"
            for B, pf in [(0.0, 5e-6), (1.0, 4e-6), (2.0, 1e-7),
                           (3.0, 3e-6), (4.0, 2e-6), (5.0, 6e-6)]
        ) + "\n"
        _, _, bpf = self._run_with_pf(pf_text=pf_text)
        self.assertIsNotNone(bpf, "Varying Pfaffian profile must yield a numeric B_crit")
        self.assertAlmostEqual(bpf, 2.0, places=4)


if __name__ == "__main__":
    # Useful when running this test directly:
    #   python3 tests/integration/test_pfaffian_degeneracy_detection.py
    unittest.main(verbosity=2)
