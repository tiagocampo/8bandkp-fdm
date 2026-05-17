"""Fortran executable runner for cross-code validation.

Wraps our Fortran executables via subprocess, reusing the star_helpers
pattern for running executables and parsing output.
"""

import os
import shutil
import subprocess
import tempfile


def run_bulk(build_dir, material, k_points, work_dir=None, timeout=120):
    """Run bandStructure executable for bulk calculation.

    Generates input.cfg, runs the executable, parses eigenvalues.

    Args:
        build_dir: path to build/ directory containing src/bandStructure
        material: material name (e.g. 'GaAs')
        k_points: list of (kx, ky, kz) tuples in 1/Angstrom
        work_dir: working directory (created if None)
        timeout: execution timeout in seconds

    Returns:
        list of (|k|, [eigenvalues]) tuples, eigenvalues in eV ascending
    """
    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="fortran_bulk_")

    exe_path = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(f"Executable not found: {exe_path}")

    # Build k-point list as waveVector/waveVectorMax/waveVectorStep
    k_mags = []
    for kx, ky, kz in k_points:
        k_mag = (kx**2 + ky**2 + kz**2) ** 0.5
        k_mags.append(k_mag)

    # Generate input.cfg for bulk
    cfg_content = _build_bulk_config(material, k_points)
    cfg_path = os.path.join(work_dir, "input.cfg")
    with open(cfg_path, "w") as f:
        f.write(cfg_content)

    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    result = subprocess.run(
        [exe_path],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"bandStructure failed (rc={result.returncode}):\n"
            f"stdout: {result.stdout[-500:]}\n"
            f"stderr: {result.stderr[-500:]}"
        )

    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        raise FileNotFoundError(f"No eigenvalues.dat in {output_dir}")

    return _parse_eigenvalues(eig_path)


def run_qw(build_dir, barrier_material, well_material,
           l_well_ang, l_barrier_ang, total_length_ang,
           k_points, fdstep=201, fdorder=2, numcb=4, numvb=8,
           work_dir=None, timeout=120):
    """Run bandStructure executable for QW calculation.

    Uses the 2-layer "last-layer-wins" pattern:
    Layer 1: barrier covers full domain
    Layer 2: well overwrites center

    Args:
        build_dir: path to build/ directory
        barrier_material: material name for barrier
        well_material: material name for well
        l_well_ang: well width in Angstrom
        l_barrier_ang: barrier width in Angstrom (one side)
        total_length_ang: total domain length in Angstrom
        k_points: list of (kx, ky) tuples in 1/Angstrom (in-plane k)
        fdstep: number of FD grid points
        fdorder: FD order (2, 4, 6, 8)
        numcb: number of CB bands to output
        numvb: number of VB bands to output
        work_dir: working directory
        timeout: execution timeout

    Returns:
        list of (|k|, [eigenvalues]) tuples, eigenvalues in eV ascending
    """
    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="fortran_qw_")

    exe_path = os.path.join(build_dir, "src", "bandStructure")
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(f"Executable not found: {exe_path}")

    # Build k-point list
    if len(k_points) == 1:
        kx, ky = k_points[0]
        k_mag = (kx**2 + ky**2)**0.5
        if k_mag < 1e-12:
            direction = "k0"
            k_max = 0.0
            n_steps = 1
        else:
            # Single non-zero k-point: use n_steps=2 to avoid div-by-zero in Fortran
            # (Fortran computes k_i = (i-1)*kmax/(nsteps-1), needs nstep >= 2)
            if abs(ky) < 1e-12:
                direction = "kx"
                k_max = abs(kx)
            elif abs(kx) < 1e-12:
                direction = "ky"
                k_max = abs(ky)
            else:
                raise ValueError(
                    f"Non-axis-aligned single k-point ({kx}, {ky}) not supported for QW."
                )
            n_steps = 2
    else:
        k_mags = [(kx**2 + ky**2)**0.5 for kx, ky in k_points]
        k_max = max(k_mags)
        n_steps = len(k_points)

        nonzero = [(kx, ky) for kx, ky in k_points if (kx**2 + ky**2) > 1e-12]
        if nonzero and all(abs(p[1]) < 1e-12 for p in nonzero):
            direction = "kx"
        elif nonzero and all(abs(p[0]) < 1e-12 for p in nonzero):
            direction = "ky"
        else:
            direction = "k0"

    # Last-layer-wins: barrier covers full domain, well overwrites center
    # Use centered coordinates (regression configs use symmetric -L/2 to L/2)
    half_total = total_length_ang / 2.0
    z_min = -half_total
    z_max = half_total
    well_start = -l_well_ang / 2.0
    well_end = l_well_ang / 2.0

    cfg_content = "\n".join([
        f"waveVector: {direction}",
        f"waveVectorMax: {k_max:.10f}",
        f"waveVectorStep: {n_steps}",
        "confinement: 1",
        f"FDstep: {fdstep}",
        f"FDorder: {fdorder}",
        "numLayers: 2",
        f"material1: {barrier_material} {z_min:.1f} {z_max:.1f} 0",
        f"material2: {well_material} {well_start:.1f} {well_end:.1f} 0",
        f"numcb: {numcb}",
        f"numvb: {numvb}",
        "ExternalField: 0  EF",
        "EFParams: 0.0",
    ]) + "\n"

    cfg_path = os.path.join(work_dir, "input.cfg")
    with open(cfg_path, "w") as f:
        f.write(cfg_content)

    output_dir = os.path.join(work_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    result = subprocess.run(
        [exe_path],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"bandStructure failed (rc={result.returncode}):\n"
            f"stdout: {result.stdout[-500:]}\n"
            f"stderr: {result.stderr[-500:]}"
        )

    eig_path = os.path.join(output_dir, "eigenvalues.dat")
    if not os.path.isfile(eig_path):
        raise FileNotFoundError(f"No eigenvalues.dat in {output_dir}")

    return _parse_eigenvalues(eig_path)


def _build_bulk_config(material, k_points):
    """Build input.cfg content for bulk band structure.

    Uses our input.cfg format: waveVector: <direction>, waveVectorMax: <kmax>,
    waveVectorStep: <nsteps>. kmax in 1/Angstrom, nsteps is integer count.
    """
    if len(k_points) == 1:
        kx, ky, kz = k_points[0]
        k_mag = (kx**2 + ky**2 + kz**2) ** 0.5

        if k_mag < 1e-12:
            direction = "k0"
            k_max = 0.0
            n_steps = 1
        elif abs(ky) < 1e-12 and abs(kz) < 1e-12:
            direction = "kx"
            k_max = abs(kx)
            n_steps = 1
        elif abs(kx) < 1e-12 and abs(kz) < 1e-12:
            direction = "ky"
            k_max = abs(ky)
            n_steps = 1
        elif abs(kx) < 1e-12 and abs(ky) < 1e-12:
            direction = "kz"
            k_max = abs(kz)
            n_steps = 1
        else:
            raise ValueError(
                f"Non-axis-aligned single k-point ({kx}, {ky}, {kz}) not supported. "
                f"Use axis-aligned k-points (only one non-zero component)."
            )
    else:
        k_mags = [(kx**2 + ky**2 + kz**2)**0.5 for kx, ky, kz in k_points]
        k_max = max(k_mags)
        n_steps = len(k_points)

        nonzero = [(kx, ky, kz) for kx, ky, kz in k_points
                    if (kx**2 + ky**2 + kz**2) > 1e-12]
        if nonzero and all(abs(p[1]) < 1e-12 and abs(p[2]) < 1e-12 for p in nonzero):
            direction = "kx"
        elif nonzero and all(abs(p[0]) < 1e-12 and abs(p[2]) < 1e-12 for p in nonzero):
            direction = "ky"
        elif nonzero and all(abs(p[0]) < 1e-12 and abs(p[1]) < 1e-12 for p in nonzero):
            direction = "kz"
        else:
            direction = "k0"

    lines = [
        f"waveVector: {direction}",
        f"waveVectorMax: {k_max:.10f}",
        f"waveVectorStep: {n_steps}",
        "confinement: 0",
        "FDstep: 101",
        "FDorder: 2",
        "numLayers: 1",
        f"material1: {material}",
        "numcb: 2",
        "numvb: 6",
        "ExternalField: 0  EF",
        "EFParams: 0.0",
    ]

    return "\n".join(lines) + "\n"


def _parse_eigenvalues(filepath):
    """Parse eigenvalues.dat, returning list of (|k|, [eigenvalues]).

    Each line: |k| eval_1 eval_2 ... eval_N
    Comment lines start with '#'.

    Returns:
        list of (float, list[float])
    """
    results = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = [float(x) for x in line.split()]
            if len(vals) >= 2:
                k = vals[0]
                evals = sorted(vals[1:])
                results.append((k, evals))
    return results
