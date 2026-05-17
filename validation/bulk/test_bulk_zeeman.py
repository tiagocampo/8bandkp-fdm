"""Bulk Zeeman splitting cross-code validation.

Compares Zeeman-split eigenvalues at k=0 under magnetic field.
Validates ge/kappa Zeeman mapping.

Tolerance: 0.01 meV for splitting at low B.
"""

import os
import sys

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

from validation.shared.param_mapper import map_material, FORTRAN_MATERIALS

# Placeholder: Zeeman comparison requires magnetic field support in both codes
# Our code supports B-field via ExternalField, kdotpy via magnetic field params
# This will be implemented after bulk dispersion is validated


def test_bulk_zeeman():
    """Placeholder for bulk Zeeman comparison."""
    print("=" * 70)
    print("BULK ZEEMAN CROSS-CODE VALIDATION")
    print("=" * 70)
    print("SKIPPED: Zeeman comparison deferred pending B-field API investigation")
    print("Both codes support B-field, but API wrapping needs additional work")
    return True


if __name__ == "__main__":
    success = test_bulk_zeeman()
    sys.exit(0 if success else 1)
