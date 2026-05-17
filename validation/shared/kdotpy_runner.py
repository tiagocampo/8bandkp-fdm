"""kdotpy library wrapper for cross-code validation.

Wraps kdotpy's Python API for programmatic eigenvalue computation.
Uses kdotpy as a library (not CLI) for clean structured access.
"""

import sys
import os

from kdotpy.config import initialize_config
from kdotpy.physparams import PhysParams
from kdotpy.vector import Vector
from kdotpy.diagonalization.diagonalization import hbulk, hz


_initialized = False


def init_kdotpy():
    """Initialize kdotpy configuration. Call once before any calculations."""
    global _initialized
    if not _initialized:
        initialize_config()
        _initialized = True


def run_bulk(material, k_points):
    """Run bulk 8-band k.p calculation for a list of k-points.

    Args:
        material: kdotpy Material object
        k_points: list of (kx, ky, kz) tuples in 1/nm

    Returns:
        list of numpy arrays, each with 8 eigenvalues in meV (ascending)
    """
    init_kdotpy()
    pp = PhysParams(
        kdim=3,
        m_layers=[material],
        l_layers=[1.0],
        norbitals=8,
        rel_strain=0.0,
    )

    results = []
    for kx, ky, kz in k_points:
        k = Vector(kx, ky, kz)
        dd = hbulk(k, pp)
        results.append(sorted(dd.eival))
    return results


def run_bulk_single(material, kx=0.0, ky=0.0, kz=0.0):
    """Run bulk calculation for a single k-point.

    Args:
        material: kdotpy Material object
        kx, ky, kz: wavevector components in 1/nm

    Returns:
        numpy array with 8 eigenvalues in meV (ascending)
    """
    return run_bulk(material, [(kx, ky, kz)])[0]


def run_qw(barrier_material, well_material, l_well_nm, l_barrier_nm,
           k_points, zres=0.25, neig=50, energy=0.0):
    """Run QW 8-band k.p calculation for a list of in-plane k-points.

    Args:
        barrier_material: kdotpy Material object for barrier
        well_material: kdotpy Material object for well
        l_well_nm: well width in nm
        l_barrier_nm: barrier width in nm (one side)
        k_points: list of (kx, ky) tuples in 1/nm (in-plane k)
        zres: grid resolution in nm
        neig: number of eigenvalues to compute
        energy: target energy for shift-and-invert (meV)

    Returns:
        list of numpy arrays, each with eigenvalues in meV (ascending)
    """
    init_kdotpy()

    total_length = 2 * l_barrier_nm + l_well_nm
    pp = PhysParams(
        kdim=2,
        m_layers=[barrier_material, well_material, barrier_material],
        l_layers=[l_barrier_nm, l_well_nm, l_barrier_nm],
        norbitals=8,
        zres=zres,
        rel_strain=0.0,
    )

    results = []
    for kx, ky in k_points:
        k = Vector(kx, ky, 0.0)
        dd = hz(k, pp, neig=neig, energy=energy)
        results.append(sorted(dd.eival))
    return results
