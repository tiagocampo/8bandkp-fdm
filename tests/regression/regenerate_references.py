#!/usr/bin/env python3
"""Regenerate regression reference data after const convention fix.

Runs each config through the appropriate executable and copies output files
to the reference data directory. Existing reference files are overwritten.

Usage: python3 regenerate_references.py <build_dir>
  build_dir -- path to build/ (contains src/ executables)
"""
import os
import sys
import shutil
import subprocess
import tempfile

# Maps test name -> (executable, config, ref_data_dir, output_files)
REGENERATION_MAP = [
    # Bulk band structure
    ("bandStructure", "bulk_gaas_kx.toml", "bulk_gaas_kx", ["eigenvalues.dat"]),
    ("bandStructure", "bulk_inas_kx.toml", "bulk_inas_kx", ["eigenvalues.dat"]),
    # QW band structure
    ("bandStructure", "qw_alsbw_gasbw_inasw.toml", "qw_alsb_gasb_inas", ["eigenvalues.dat"]),
    # Gfactor (compares stdout, captured as output.txt)
    ("gfactorCalculation", "gfactor_qw_cb.toml", "gfactor_qw_cb", ["output.txt"]),
    ("gfactorCalculation", "gfactor_bulk_gaas_vb.toml", "gfactor_bulk_gaas_vb", ["output.txt"]),
    ("gfactorCalculation", "gfactor_bulk_gaasw_cb.toml", "gfactor_bulk_gaasw_cb", ["output.txt"]),
    ("gfactorCalculation", "gfactor_qw_vb.toml", "gfactor_qw_vb", ["output.txt"]),
    # Landau
    ("bandStructure", "landau_bulk_GaAs.toml", "landau_bulk_GaAs", ["eigenvalues.dat"]),
    ("bandStructure", "landau_bulk_InAs.toml", "landau_bulk_InAs", ["eigenvalues.dat"]),
    ("bandStructure", "landau_bulk_InAs_Bsweep.toml", "landau_bulk_InAs_Bsweep", ["eigenvalues.dat"]),
    # SC loop
    ("bandStructure", "sc_gaas_alas_qw.toml", "sc_gaas_alas_qw", ["eigenvalues.dat"]),
    ("bandStructure", "sc_qw_inas_alsb.toml", "sc_qw_inas_alsb", ["eigenvalues.dat"]),
    ("bandStructure", "sc_bulk_gaas_doped.toml", "sc_bulk_gaas_doped",
     ["eigenvalues.dat", "sc_charge.dat", "sc_potential_profile.dat"]),
    ("bandStructure", "sc_mod_doped_gaas_algaas.toml", "sc_mod_doped_gaas_algaas", ["eigenvalues.dat"]),
    ("bandStructure", "sc_qcse_gaas_algaas.toml", "sc_qcse_gaas_algaas", ["eigenvalues.dat"]),
    ("bandStructure", "sc_qcse_gaas_algaas_ef.toml", "sc_qcse_gaas_algaas_ef", ["eigenvalues.dat"]),
    # Wire
    ("bandStructure", "wire_gaas_rectangle.toml", "wire_gaas_rectangle", ["eigenvalues.dat"]),
    ("bandStructure", "wire_inas_gaas_strain.toml", "wire_inas_gaas_core_shell", ["eigenvalues.dat"]),
    ("bandStructure", "wire_gaas_hexagon.toml", "wire_gaas_hexagon",
     ["eigenvalues.dat", "output.txt"]),
    # SC wire
    ("bandStructure", "sc_wire_gaas.toml", "sc_wire_gaas",
     ["eigenvalues.dat", "output.txt"]),
    # Wire gfactor
    ("gfactorCalculation", "wire_insb_gfactor.toml", "wire_insb_gfactor", ["output.txt"]),
    # Optics
    ("opticalProperties", "qw_optics_commutator.toml", "qw_optics_commutator",
     ["absorption_TE.dat", "absorption_TM.dat"]),
]


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <build_dir>")
        sys.exit(1)

    build_dir = os.path.abspath(sys.argv[1])
    script_dir = os.path.dirname(os.path.abspath(__file__))
    configs_dir = os.path.join(script_dir, "configs")
    data_dir = os.path.join(script_dir, "data")

    # Flush output for each print
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, line_buffering=True)

    n_ok = 0
    n_fail = 0

    for exe_name, config, ref_dir, output_files in REGENERATION_MAP:
        exe_path = os.path.join(build_dir, "src", exe_name)
        config_path = os.path.join(configs_dir, config)
        ref_path = os.path.join(data_dir, ref_dir)

        label = f"{exe_name} + {config} -> {ref_dir}"
        print(f"Regenerating: {label}")

        if not os.path.isfile(exe_path):
            print(f"  SKIP: executable not found: {exe_path}")
            continue

        if not os.path.isfile(config_path):
            print(f"  SKIP: config not found: {config_path}")
            continue

        work_dir = tempfile.mkdtemp(prefix="regen_")
        try:
            shutil.copy(config_path, os.path.join(work_dir, "input.toml"))
            os.makedirs(os.path.join(work_dir, "output"), exist_ok=True)

            result = subprocess.run(
                [exe_path],
                cwd=work_dir,
                capture_output=True,
                text=True,
                timeout=600,
            )

            if result.returncode != 0:
                print(f"  FAIL: {exe_name} returned {result.returncode}")
                if result.stderr:
                    print(f"    stderr: {result.stderr[:200]}")
                n_fail += 1
                continue

            os.makedirs(ref_path, exist_ok=True)

            # Save stdout as output.txt if needed (gfactor tests compare stdout)
            if "output.txt" in output_files:
                stdout_path = os.path.join(work_dir, "output", "output.txt")
                with open(stdout_path, "w") as f:
                    f.write(result.stdout)

            copied = 0
            for fname in output_files:
                src = os.path.join(work_dir, "output", fname)
                dst = os.path.join(ref_path, fname)
                if os.path.isfile(src):
                    shutil.copy2(src, dst)
                    copied += 1
                else:
                    print(f"  WARN: {fname} not produced")

            if copied > 0:
                print(f"  OK: {copied}/{len(output_files)} files regenerated")
                n_ok += 1
            else:
                print(f"  FAIL: no output files produced")
                n_fail += 1

        except subprocess.TimeoutExpired:
            print(f"  FAIL: timeout")
            n_fail += 1
        except Exception as e:
            print(f"  FAIL: {e}")
            n_fail += 1
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)

    print(f"\nSummary: {n_ok} regenerated, {n_fail} failed")


if __name__ == "__main__":
    main()
