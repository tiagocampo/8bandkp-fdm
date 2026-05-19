"""Parameter mapping from our Fortran code format to kdotpy format.

Hardcoded material parameters from parameters.f90, converted to kdotpy
conventions using the mapping formulas resolved in the plan.

Key mappings:
  F = (1/meff - 1) / 2        (CB diagonal equality)
  kappa = 0.0                  (eliminates interface discrepancy)
  P = sqrt(EP * HBARM0_KDOTPY * 1000)  (meV*nm)
  Ec = Eg * 1000               (meV, with Ev = 0)
  Ev = 0.0                     (energy reference)
  delta_so = deltaSO * 1000    (meV)
  strain_C1 = ac * 1000        (meV)
  strain_Dd = av * 1000        (meV)
  strain_Du = -1.5 * b_dp * 1000  (meV, convention)
  strain_Duprime = -0.5 * sqrt(3) * d_dp * 1000  (meV, convention)
  a = a0 / 10                  (nm)
"""

import math


# Materials from parameters.f90 (Vurgaftman 2001, Winkler 2003 for W-variants).
# Values in our native units: eV, Angstrom, kbar.
FORTRAN_MATERIALS = {
    "GaAs": {
        "meff": 0.067, "EP": 28.8, "Eg": 1.519, "deltaSO": 0.341,
        "gamma1": 6.98, "gamma2": 2.06, "gamma3": 2.93,
        "EV": -0.8, "EC": 0.719, "eps0": 12.90,
        "C11": 1221.0, "C12": 566.0, "C44": 599.0, "a0": 5.65325,
        "ac": -7.17, "av": 1.16, "b_dp": -2.0, "d_dp": -4.8,
    },
    "GaAsW": {
        "meff": 0.0665, "EP": 28.89, "Eg": 1.519, "deltaSO": 0.341,
        "gamma1": 6.85, "gamma2": 2.10, "gamma3": 2.90,
        "EV": -0.8, "EC": 0.719, "eps0": 12.90,
        "C11": 1221.0, "C12": 566.0, "C44": 599.0, "a0": 5.65325,
        "ac": -7.17, "av": 1.16, "b_dp": -2.0, "d_dp": -4.8,
    },
    "Al20Ga80As": {
        "meff": 0.0836, "EP": 27.26, "Eg": 1.835, "deltaSO": 0.3288,
        "gamma1": 6.336, "gamma2": 1.812, "gamma3": 2.628,
        "EV": -0.906, "EC": 0.929, "eps0": 12.332,
        "C11": 1221.0*0.8 + 1250.0*0.2, "C12": 566.0*0.8 + 534.0*0.2,
        "C44": 599.0*0.8 + 542.0*0.2, "a0": 5.65325*0.8 + 5.6611*0.2,
        "ac": -7.17*0.8 + (-5.64)*0.2, "av": 1.16*0.8 + 2.47*0.2,
        "b_dp": -2.0*0.8 + (-2.3)*0.2, "d_dp": -4.8*0.8 + (-3.4)*0.2,
    },
    "Al15Ga85As": {
        "meff": 0.0795, "EP": 27.645, "Eg": 1.756, "deltaSO": 0.3319,
        "gamma1": 6.497, "gamma2": 1.874, "gamma3": 2.7035,
        "EV": -0.8795, "EC": 0.8765, "eps0": 12.474,
        "C11": 1221.0*0.85 + 1250.0*0.15, "C12": 566.0*0.85 + 534.0*0.15,
        "C44": 599.0*0.85 + 542.0*0.15, "a0": 5.65325*0.85 + 5.6611*0.15,
        "ac": -7.17*0.85 + (-5.64)*0.15, "av": 1.16*0.85 + 2.47*0.15,
        "b_dp": -2.0*0.85 + (-2.3)*0.15, "d_dp": -4.8*0.85 + (-3.4)*0.15,
    },
    "Al30Ga70As": {
        "meff": 0.093, "EP": 26.32, "Eg": 1.977, "deltaSO": 0.353,
        "gamma1": 6.107, "gamma2": 1.773, "gamma3": 2.543,
        "EV": -0.959, "EC": 1.018, "eps0": 12.048,
        "C11": 1221.0*0.70 + 1250.0*0.30, "C12": 566.0*0.70 + 534.0*0.30,
        "C44": 599.0*0.70 + 542.0*0.30, "a0": 5.65325*0.70 + 5.6611*0.30,
        "ac": -7.17*0.70 + (-5.64)*0.30, "av": 1.16*0.70 + 2.47*0.30,
        "b_dp": -2.0*0.70 + (-2.3)*0.30, "d_dp": -4.8*0.70 + (-3.4)*0.30,
    },
    "InAs": {
        "meff": 0.026, "EP": 21.5, "Eg": 0.417, "deltaSO": 0.39,
        "gamma1": 20.0, "gamma2": 8.5, "gamma3": 9.2,
        "EV": -0.59, "EC": -0.173, "eps0": 15.15,
        "C11": 832.9, "C12": 452.6, "C44": 395.9, "a0": 6.0583,
        "ac": -5.08, "av": 1.00, "b_dp": -1.8, "d_dp": -3.6,
    },
    "InAsW": {
        "meff": 0.0229, "EP": 22.2, "Eg": 0.418, "deltaSO": 0.38,
        "gamma1": 20.4, "gamma2": 8.3, "gamma3": 9.1,
        "EV": -0.59, "EC": -0.172, "eps0": 15.15,
        "C11": 832.9, "C12": 452.6, "C44": 395.9, "a0": 6.0583,
        "ac": -5.08, "av": 1.00, "b_dp": -1.8, "d_dp": -3.6,
    },
    "AlAs": {
        "meff": 0.15, "EP": 21.1, "Eg": 3.099, "deltaSO": 0.28,
        "gamma1": 3.76, "gamma2": 0.82, "gamma3": 1.42,
        "EV": -1.33, "EC": 1.769, "eps0": 10.06,
        "C11": 1250.0, "C12": 534.0, "C44": 542.0, "a0": 5.6611,
        "ac": -5.64, "av": 2.47, "b_dp": -2.3, "d_dp": -3.4,
    },
    "GaSb": {
        "meff": 0.039, "EP": 27.0, "Eg": 0.812, "deltaSO": 0.76,
        "gamma1": 13.4, "gamma2": 4.7, "gamma3": 6.0,
        "EV": 0.0, "EC": 0.812, "eps0": 15.70,
        "C11": 884.2, "C12": 402.4, "C44": 432.3, "a0": 6.0959,
        "ac": -7.5, "av": 0.8, "b_dp": -2.0, "d_dp": -4.7,
    },
    "AlSb": {
        "meff": 0.14, "EP": 18.7, "Eg": 2.386, "deltaSO": 0.676,
        "gamma1": 5.18, "gamma2": 1.19, "gamma3": 1.97,
        "EV": -0.41, "EC": 1.976, "eps0": 12.04,
        "C11": 876.5, "C12": 434.1, "C44": 407.6, "a0": 6.1355,
        "ac": -4.5, "av": 1.4, "b_dp": -1.35, "d_dp": -4.3,
    },
    "InSb": {
        "meff": 0.0135, "EP": 23.3, "Eg": 0.235, "deltaSO": 0.81,
        "gamma1": 34.8, "gamma2": 15.5, "gamma3": 16.5,
        "EV": 0.0, "EC": 0.235, "eps0": 16.80,
        "C11": 684.7, "C12": 373.5, "C44": 311.1, "a0": 6.4794,
        "ac": -6.94, "av": 0.36, "b_dp": -2.0, "d_dp": -4.7,
    },
    "InP": {
        "meff": 0.0795, "EP": 20.7, "Eg": 1.4236, "deltaSO": 0.108,
        "gamma1": 5.08, "gamma2": 1.60, "gamma3": 2.10,
        "EV": -0.70, "EC": 0.7236, "eps0": 12.61,
        "C11": 1011.0, "C12": 561.0, "C44": 456.0, "a0": 5.8687,
        "ac": -6.0, "av": 0.6, "b_dp": -2.0, "d_dp": -5.0,
    },
    "InPW": {
        "meff": 0.0803, "EP": 20.56, "Eg": 1.423, "deltaSO": 0.11,
        "gamma1": 4.95, "gamma2": 1.65, "gamma3": 2.35,
        "EV": -0.70, "EC": 0.723, "eps0": 12.61,
        "C11": 1011.0, "C12": 561.0, "C44": 456.0, "a0": 5.8687,
        "ac": -6.35, "av": 1.70, "b_dp": -2.0, "d_dp": -5.6,
    },
    "InSbW": {
        "meff": 0.0139, "EP": 24.4, "Eg": 0.237, "deltaSO": 0.81,
        "gamma1": 37.1, "gamma2": 16.5, "gamma3": 17.7,
        "EV": 0.0, "EC": 0.237, "eps0": 16.80,
        "C11": 684.7, "C12": 373.5, "C44": 311.1, "a0": 6.4794,
        "ac": -6.94, "av": 0.36, "b_dp": -2.0, "d_dp": -4.7,
    },
    "AlAsW": {
        "meff": 0.15, "EP": 21.12, "Eg": 3.13, "deltaSO": 0.3,
        "gamma1": 3.25, "gamma2": 0.65, "gamma3": 1.21,
        "EV": -1.33, "EC": 1.800, "eps0": 10.06,
        "C11": 1250.0, "C12": 534.0, "C44": 542.0, "a0": 5.6611,
        "ac": -5.64, "av": 2.47, "b_dp": -2.3, "d_dp": -3.4,
    },
    "GaSbW": {
        "meff": 0.041, "EP": 22.37, "Eg": 0.812, "deltaSO": 0.76,
        "gamma1": 13.4, "gamma2": 4.7, "gamma3": 5.7,
        "EV": -0.03, "EC": 0.782, "eps0": 15.70,
        "C11": 884.2, "C12": 402.4, "C44": 432.3, "a0": 6.0959,
        "ac": -7.5, "av": 0.8, "b_dp": -2.0, "d_dp": -4.7,
    },
    "AlSbW": {
        "meff": 0.12, "EP": 18.8, "Eg": 2.384, "deltaSO": 0.673,
        "gamma1": 4.15, "gamma2": 1.01, "gamma3": 1.71,
        "EV": -0.41, "EC": 1.974, "eps0": 12.04,
        "C11": 876.5, "C12": 434.1, "C44": 407.6, "a0": 6.1355,
        "ac": -4.5, "av": 1.4, "b_dp": -1.35, "d_dp": -4.3,
    },
    "Ga47In53AsW": {
        "meff": 0.038, "EP": 25.26, "Eg": 0.8166, "deltaSO": 0.362,
        "gamma1": 11.97, "gamma2": 4.36, "gamma3": 5.15,
        "EV": -0.689, "EC": 0.128, "eps0": 14.0925,
        "C11": 1221.0*0.47 + 832.9*0.53, "C12": 566.0*0.47 + 452.6*0.53,
        "C44": 599.0*0.47 + 395.9*0.53, "a0": 5.65325*0.47 + 6.0583*0.53,
        "ac": -7.17*0.47 + (-5.08)*0.53, "av": 1.16*0.47 + 1.00*0.53,
        "b_dp": -2.0*0.47 + (-1.8)*0.53, "d_dp": -4.8*0.47 + (-3.6)*0.53,
    },
    "Al47In53AsW": {
        "meff": 0.0779, "EP": 21.7, "Eg": 1.693, "deltaSO": 0.342,
        "gamma1": 6.17, "gamma2": 1.62, "gamma3": 2.31,
        "EV": -0.938, "EC": 0.755, "eps0": 12.7577,
        "C11": 1250.0*0.47 + 832.9*0.53, "C12": 534.0*0.47 + 452.6*0.53,
        "C44": 542.0*0.47 + 395.9*0.53, "a0": 5.6611*0.47 + 6.0583*0.53,
        "ac": -5.64*0.47 + (-5.08)*0.53, "av": 2.47*0.47 + 1.00*0.53,
        "b_dp": -2.3*0.47 + (-1.8)*0.53, "d_dp": -3.4*0.47 + (-3.6)*0.53,
    },
    # II-VI materials (Pfeuffer-Jeschke PhD thesis, 2000; Novik et al., PRB 72, 2005)
    # HgTe: deltaSO=1.003 eV from Pfeuffer-Jeschke; C11/C12/C44 from Landolt-Bornstein
    "HgTe": {
        "meff": 1.0, "EP": 18.8, "Eg": -0.303, "deltaSO": 1.003,
        "gamma1": 4.1, "gamma2": 0.5, "gamma3": 1.3,
        "EV": 0.0, "EC": -0.303, "eps0": 15.0,
        "C11": 532.0, "C12": 368.0, "C44": 201.0, "a0": 6.461,
        "ac": -0.48, "av": -0.35, "b_dp": 0.44, "d_dp": -1.7,
    },
    # CdTe: deltaSO=0.91 eV from Pfeuffer-Jeschke; C11/C12/C44 from Landolt-Bornstein
    "CdTe": {
        "meff": 1.2195, "EP": 18.8, "Eg": 1.606, "deltaSO": 0.91,
        "gamma1": 5.0, "gamma2": 1.3, "gamma3": 2.1,
        "EV": 0.0, "EC": 1.606, "eps0": 10.0,
        "C11": 532.0, "C12": 368.0, "C44": 201.0, "a0": 6.481,
        "ac": -0.38, "av": -0.17, "b_dp": 0.30, "d_dp": -4.9,
    },
}


def map_material(mat_name, qw_mode=False):
    """Map our material parameters to a kdotpy Material object.

    Args:
        mat_name: material name (must be in FORTRAN_MATERIALS)
        qw_mode: if True, use actual EC/EV for band offsets (needed for QW);
                 if False, use Ec=Eg*1000, Ev=0 (works for bulk k=0)

    Returns:
        kdotpy Material object with converted parameters

    Raises:
        KeyError: if material name is not in our database
    """
    if mat_name not in FORTRAN_MATERIALS:
        raise KeyError(
            f"Unknown material '{mat_name}'. "
            f"Available: {sorted(FORTRAN_MATERIALS.keys())}"
        )

    from kdotpy.materials import Material
    from kdotpy.physconst import hbarm0 as HBARM0_KDOTPY

    m = FORTRAN_MATERIALS[mat_name]

    F = (1.0 / m["meff"] - 1.0) / 2.0
    P = math.sqrt(m["EP"] * 1000.0 * HBARM0_KDOTPY)
    delta_so = m["deltaSO"] * 1000.0  # meV
    a_nm = m["a0"] / 10.0  # Angstrom -> nm

    if qw_mode:
        # Use actual EC/EV from parameters.f90 for correct band offsets
        Ec = m["EC"] * 1000.0  # eV -> meV
        Ev = m["EV"] * 1000.0  # eV -> meV
    else:
        # Bulk k=0 reference: Ev=0, Ec=Eg (gap is correct, offsets don't matter)
        Ec = m["Eg"] * 1000.0  # meV, with Ev=0
        Ev = 0.0
    param = {
        "F": F,
        "gamma1": m["gamma1"],
        "gamma2": m["gamma2"],
        "gamma3": m["gamma3"],
        "Ec": Ec,
        "Ev": Ev,
        "delta_so": delta_so,
        "P": P,
        "kappa": 0.0,
        "q": 0.0,
        "ge": 2.0,
        "a": a_nm,
        "elasticity_c11": m["C11"],
        "elasticity_c12": m["C12"],
        "elasticity_c44": m["C44"],
    }

    if "ac" in m:
        param["strain_C1"] = m["ac"] * 1000.0
        param["strain_Dd"] = m["av"] * 1000.0
        param["strain_Du"] = -0.75 * m["b_dp"] * 1000.0
        param["strain_Duprime"] = -0.5 * math.sqrt(3.0) * m["d_dp"] * 1000.0

    return Material(mat_name, param=param)


def get_fortran_params(mat_name):
    """Get our Fortran material parameters as a dict.

    Args:
        mat_name: material name

    Returns:
        dict with our native parameters (eV, Angstrom units)
    """
    if mat_name not in FORTRAN_MATERIALS:
        raise KeyError(f"Unknown material '{mat_name}'")
    return dict(FORTRAN_MATERIALS[mat_name])


def list_materials():
    """Return sorted list of available material names."""
    return sorted(FORTRAN_MATERIALS.keys())
