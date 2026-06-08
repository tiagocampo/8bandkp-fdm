#!/usr/bin/env python3
"""Convert old label:value config files to TOML format for the 8-band k.p solver.

Reads a .cfg file in the old positional label:value format and writes a
corresponding .toml file that the toml-f based input parser accepts.

Usage:
    python3 convert_cfg_to_toml.py path/to/config.cfg [path/to/output.toml]
    python3 convert_cfg_to_toml.py --batch tests/regression/configs/
"""

import sys
import os
import re
from pathlib import Path


def parse_cfg(filepath):
    """Parse a .cfg file into a list of (label, value) tuples."""
    lines = []
    with open(filepath, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('!'):
                continue
            # Split on first colon
            if ':' in line:
                label, _, value = line.partition(':')
                lines.append((label.strip(), value.strip()))
    return lines


def to_toml_string(s):
    """Escape a string for TOML."""
    if '"' in s:
        return '"' + s.replace('\\', '\\\\').replace('"', '\\"') + '"'
    return '"' + s + '"'


def to_toml_value(v):
    """Format a value for TOML output."""
    if v.lower() in ('true', 't'):
        return 'true'
    if v.lower() in ('false', 'f'):
        return 'false'
    # Try int
    try:
        iv = int(v)
        return str(iv)
    except ValueError:
        pass
    # Try float
    try:
        fv = float(v)
        return repr(fv)
    except ValueError:
        pass
    # String
    return to_toml_string(v)


def convert_cfg_to_toml(cfg_path, toml_path=None):
    """Convert a .cfg file to .toml format."""
    lines = parse_cfg(cfg_path)

    if toml_path is None:
        toml_path = str(Path(cfg_path).with_suffix('.toml'))

    # Collect into sections
    confinement = None
    fdorder = None
    fdstep = None

    # Wave vector
    wv_mode = None
    wv_max = None
    wv_nsteps = None

    # Bands
    numcb = None
    numvb = None

    # Materials (bulk/QW)
    materials = []  # list of dicts with 'name', optional 'z_min', 'z_max'

    # Wire
    wire_nx = None
    wire_ny = None
    wire_dx = None
    wire_dy = None
    wire_shape = None
    wire_radius = None
    wire_width = None
    wire_height = None
    wire_polygon_verts = []

    # Wire regions
    regions = []

    # Landau
    landau_nx = None
    landau_width = None
    landau_sweep = None

    # External field
    ef_enabled = None  # the integer before "EF"
    ef_type = None
    ef_value = None

    # B field
    b_field_components = None
    b_field_g_factor = None
    b_sweep = None

    # SC
    sc_enabled = None
    sc_max_iter = None
    sc_tolerance = None
    sc_mixing_alpha = None
    sc_diis = None
    sc_temperature = None
    sc_fermi_mode = None
    sc_fermi_level = None
    sc_num_kpar = None
    sc_kpar_max = None
    sc_bc_type = None
    sc_bc_left = None
    sc_bc_right = None

    # Doping (uniform: doping<N>: ND NA; delta: delta<N>: NS fwhm pos)
    uniform_doping = []  # list of (ND, NA)
    delta_doping = []    # list of (NS, fwhm, pos)

    # Strain
    strain_enabled = None
    strain_ref = None
    strain_solver = None
    piezo = None
    strain_substrate_val = None  # numeric substrate value (bulk)

    # Topology
    topo_enabled = None
    topo_mode = None
    topo_compute_chern = None
    topo_compute_hall = None
    topo_qwz_u = None
    topo_compute_z2 = None
    topo_z2_method = None
    topo_extract_edge_states = None
    topo_edge_E_window = None
    topo_compute_ldos = None
    topo_ldos_eta = None
    topo_ldos_E_range = []  # may have two entries
    topo_ldos_num_E = None

    # Gap sweep
    topo_compute_gap_sweep = None
    gap_sweep_B = None
    gap_sweep_mu = None
    sweep_model = None

    # Conductance
    topo_compute_conductance = None
    conductance_method = None
    berry_nk = None
    landauer_energy = None

    # Spectral
    topo_compute_spectral = None
    spectral_k_grid = None
    spectral_E_grid = None

    # BdG
    bdg_enabled = None
    bdg_mu = None
    bdg_delta_0 = None
    bdg_g_factor = None
    bdg_B_vec = None
    bdg_gauge = None
    bdg_kz = None
    bdg_self_consistent = None

    # Optics
    optics_enabled = None
    opt_linewidth_lorentzian = None
    opt_linewidth_gaussian = None
    opt_refractive_index = None
    opt_Emin = None
    opt_Emax = None
    opt_NEnergyPoints = None
    opt_temperature = None
    opt_carrier_density = None
    opt_gain = None
    opt_gain_carrier_density = None
    opt_ISBT = None
    opt_spontaneous = None
    opt_spin_resolved = None

    # Exciton
    exciton_enabled = None
    exciton_method = None

    # Scattering
    scattering_enabled = None
    phonon_energy = None
    eps_inf = None
    eps_0 = None

    # Solver (converted from legacy FEAST)
    feast_emin = None
    feast_emax = None
    feast_m0 = None

    # G-factor
    which_band = None
    band_idx = None

    # Parse all lines
    i = 0
    while i < len(lines):
        label, value = lines[i]

        if label == 'waveVector':
            wv_mode = value
        elif label == 'waveVectorMax':
            wv_max = value
        elif label == 'waveVectorStep':
            # Old format: waveVectorStep = number of k-points (integer count)
            # New TOML: nsteps = integer count, step = real step size
            wv_nsteps = value
        elif label == 'waveVectorNsteps':
            wv_nsteps = value
        elif label == 'confinement':
            confinement = int(value)
        elif label == 'FDstep':
            fdstep = value
        elif label == 'FDorder':
            fdorder = value
        elif label == 'numcb':
            numcb = value
        elif label == 'numvb':
            numvb = value
        elif label.startswith('material'):
            # material<N>: Name [z_min z_max [lattice_mismatch]]
            parts = value.split()
            mat = {'name': parts[0]}
            if len(parts) >= 3:
                mat['z_min'] = parts[1]
                mat['z_max'] = parts[2]
            if len(parts) >= 4:
                mat['lattice_mismatch'] = parts[3]
            materials.append(mat)
        elif label == 'ExternalField':
            # ExternalField: <0|1> EF
            parts = value.split()
            ef_enabled = int(parts[0])
            if len(parts) >= 2:
                ef_type = parts[1]
        elif label == 'EFParams':
            ef_value = value
        elif label == 'wire_nx':
            wire_nx = value
        elif label == 'wire_ny':
            wire_ny = value
        elif label == 'wire_dx':
            wire_dx = value
        elif label == 'wire_dy':
            wire_dy = value
        elif label == 'wire_shape':
            wire_shape = value
        elif label == 'wire_radius':
            wire_radius = value
        elif label == 'wire_width':
            wire_width = value
        elif label == 'wire_height':
            wire_height = value
        elif label == 'wire_polygon':
            # wire_polygon: <nverts> x1 y1 x2 y2 ...
            pass  # handled by vertex lines
        elif label == 'vertex':
            wire_polygon_verts.append(value.split())
        elif label == 'numRegions':
            pass  # derived from region count
        elif label == 'region':
            # region: Material inner outer
            parts = value.split()
            regions.append({
                'material': parts[0],
                'inner': parts[1],
                'outer': parts[2]
            })
        elif label == 'landau_nx':
            landau_nx = value
        elif label == 'landau_width':
            landau_width = value
        elif label == 'landau_sweep':
            landau_sweep = value
        elif label == 'b_field':
            b_field_components = value.split()
        elif label == 'g_factor':
            if bdg_enabled == 'T':
                bdg_g_factor = value
            else:
                b_field_g_factor = value
        elif label == 'b_sweep':
            b_sweep = value.split()
        elif label == 'SC':
            sc_enabled = value
        elif label in ('max_iter', 'SC_max_iter'):
            sc_max_iter = value
        elif label in ('tolerance', 'SC_tolerance'):
            sc_tolerance = value
        elif label in ('mixing_alpha', 'SC_mixing'):
            sc_mixing_alpha = value
        elif label in ('diis_history', 'SC_diis'):
            sc_diis = value
        elif label in ('temperature', 'SC_temperature'):
            sc_temperature = value
        elif label in ('fermi_mode', 'SC_fermi_mode'):
            sc_fermi_mode = value
        elif label in ('fermi_level', 'SC_fermi_level'):
            sc_fermi_level = value
        elif label in ('num_kpar', 'SC_num_kpar'):
            sc_num_kpar = value
        elif label in ('kpar_max', 'SC_kpar_max'):
            sc_kpar_max = value
        elif label in ('bc_type', 'SC_bc'):
            sc_bc_type = value
        elif label in ('bc_left', 'SC_bc_left'):
            sc_bc_left = value
        elif label in ('bc_right', 'SC_bc_right'):
            sc_bc_right = value
        elif re.match(r'^doping\d+$', label):
            # doping<N>: ND NA
            parts = value.split()
            nd = parts[0] if len(parts) > 0 else '0.0'
            na = parts[1] if len(parts) > 1 else '0.0'
            uniform_doping.append((nd, na))
        elif re.match(r'^delta\d+$', label):
            # delta<N>: NS fwhm pos
            parts = value.split()
            ns = parts[0] if len(parts) > 0 else '0.0'
            fwhm = parts[1] if len(parts) > 1 else '10.0'
            pos = parts[2] if len(parts) > 2 else '0.0'
            delta_doping.append((ns, fwhm, pos))
        elif label in ('strain', 'Strain'):
            strain_enabled = value
        elif label in ('strain_ref', 'strainSubstrate'):
            if value.replace('.', '').replace('-', '').replace('e', '').replace('E', '').replace('+', '').isdigit():
                strain_substrate_val = value
                if strain_enabled is None:
                    strain_enabled = 'T'  # strainSubstrate implicitly enables strain
            else:
                strain_ref = value
                if strain_enabled is None:
                    strain_enabled = 'T'
        elif label == 'strain_solver':
            strain_solver = value
        elif label == 'piezo':
            piezo = value
        elif label in ('topology', 'Topology'):
            topo_enabled = value
        elif label == 'mode':
            topo_mode = value
        elif label == 'compute_chern':
            topo_compute_chern = value
        elif label == 'compute_hall':
            topo_compute_hall = value
        elif label == 'qwz_u':
            topo_qwz_u = value
        elif label == 'compute_z2':
            topo_compute_z2 = value
        elif label == 'z2_method':
            topo_z2_method = value
        elif label == 'extract_edge_states':
            topo_extract_edge_states = value
        elif label == 'edge_E_window':
            topo_edge_E_window = value
        elif label == 'compute_ldos':
            topo_compute_ldos = value
        elif label == 'ldos_eta':
            topo_ldos_eta = value
        elif label == 'ldos_E_range':
            topo_ldos_E_range.append(value)
        elif label == 'ldos_num_E':
            topo_ldos_num_E = value
        elif label == 'compute_gap_sweep':
            topo_compute_gap_sweep = value
        elif label == 'gap_sweep_B':
            gap_sweep_B = value.split()
        elif label == 'gap_sweep_mu':
            gap_sweep_mu = value.split()
        elif label == 'sweep_model':
            sweep_model = value
        elif label == 'compute_conductance':
            topo_compute_conductance = value
        elif label == 'conductance_method':
            conductance_method = value
        elif label == 'berry_nk':
            berry_nk = value
        elif label == 'landauer_energy':
            landauer_energy = value
        elif label == 'compute_spectral':
            topo_compute_spectral = value
        elif label == 'spectral_k_grid':
            spectral_k_grid = value.split()
        elif label == 'spectral_E_grid':
            spectral_E_grid = value.split()
        elif label in ('bdg', 'Bdg'):
            bdg_enabled = value
        elif label == 'mu':
            bdg_mu = value
        elif label == 'delta_0':
            bdg_delta_0 = value
        elif label == 'B_vec':
            bdg_B_vec = value.split()
        elif label == 'gauge':
            bdg_gauge = value
        elif label == 'kz':
            bdg_kz = value
        elif label == 'self_consistent':
            bdg_self_consistent = value
        elif label in ('Optics', 'optics'):
            optics_enabled = value
        elif label == 'LinewidthLorentzian':
            opt_linewidth_lorentzian = value
        elif label == 'LinewidthGaussian':
            opt_linewidth_gaussian = value
        elif label == 'RefractiveIndex':
            opt_refractive_index = value
        elif label == 'Emin':
            opt_Emin = value
        elif label == 'Emax':
            opt_Emax = value
        elif label == 'NEnergyPoints':
            opt_NEnergyPoints = value
        elif label == 'Temperature':
            opt_temperature = value
        elif label == 'CarrierDensity':
            opt_carrier_density = value
        elif label == 'Gain':
            opt_gain = value
        elif label == 'GainCarrierDensity':
            opt_gain_carrier_density = value
        elif label == 'ISBT':
            opt_ISBT = value
        elif label == 'SpontaneousEnabled':
            opt_spontaneous = value
        elif label == 'SpinResolved':
            opt_spin_resolved = value
        elif label in ('Exciton', 'exciton'):
            exciton_enabled = value
        elif label == 'ExcitonMethod':
            exciton_method = value
        elif label == 'method':
            # Context-dependent: exciton method or topology mode
            # Use recently-seen-section heuristic
            if exciton_enabled is not None and exciton_enabled.upper() in ('T', 'TRUE', '1'):
                exciton_method = value
            elif topo_enabled is not None and topo_enabled.upper() in ('T', 'TRUE', '1'):
                # If topology mode already set by a 'mode' line, this is exciton
                if topo_mode is not None:
                    exciton_method = value
                else:
                    topo_mode = value
            else:
                exciton_method = value
        elif label in ('Scattering', 'scattering'):
            scattering_enabled = value
        elif label == 'PhononEnergy':
            phonon_energy = value
        elif label == 'EpsInf':
            eps_inf = value
        elif label == 'Eps0':
            eps_0 = value
        elif label == 'feast_emin':
            feast_emin = value
        elif label == 'feast_emax':
            feast_emax = value
        elif label == 'feast_m0':
            feast_m0 = value
        elif label == 'whichBand':
            which_band = value
        elif label == 'bandIdx':
            band_idx = value
        elif label == 'numLayers':
            pass  # derived from material count

        i += 1

    # ---- Build TOML output ----
    out = []

    # Confinement
    conf_map = {0: 'bulk', 1: 'qw', 2: 'wire', 3: 'landau'}
    conf_str = conf_map.get(confinement, 'bulk')
    out.append(f'confinement = {to_toml_string(conf_str)}')

    # FDorder
    if fdorder is not None:
        out.append(f'FDorder = {fdorder}')

    # fd_step
    if fdstep is not None:
        out.append(f'fd_step = {fdstep}')

    # G-factor (top-level — must come BEFORE any [section] or [[array]] headers)
    if which_band is not None:
        out.append(f'which_band = {which_band}')
    if band_idx is not None:
        out.append(f'band_idx = {band_idx}')

    # Wave vector
    out.append('')
    out.append('[wave_vector]')
    if wv_mode is not None:
        out.append(f'mode = {to_toml_string(wv_mode)}')
    if wv_max is not None:
        out.append(f'max = {wv_max}')
    if wv_nsteps is not None:
        out.append(f'nsteps = {wv_nsteps}')

    # Bands
    out.append('')
    out.append('[bands]')
    if numcb is not None:
        out.append(f'num_cb = {numcb}')
    if numvb is not None:
        out.append(f'num_vb = {numvb}')

    # Materials
    out.append('')
    for mat in materials:
        out.append('[[material]]')
        out.append(f'name = {to_toml_string(mat["name"])}')
        if 'z_min' in mat:
            out.append(f'z_min = {mat["z_min"]}')
        if 'z_max' in mat:
            out.append(f'z_max = {mat["z_max"]}')
        out.append('')

    # Wire section
    if conf_str == 'wire':
        out.append('[wire]')
        if wire_nx is not None:
            out.append(f'nx = {wire_nx}')
        if wire_ny is not None:
            out.append(f'ny = {wire_ny}')
        if wire_dx is not None:
            out.append(f'dx = {wire_dx}')
        if wire_dy is not None:
            out.append(f'dy = {wire_dy}')

        # Wire geometry
        if wire_shape is not None:
            out.append('')
            out.append('[wire.geometry]')
            out.append(f'shape = {to_toml_string(wire_shape)}')
            if wire_shape in ('circle', 'hexagon') and wire_radius is not None:
                out.append(f'radius = {wire_radius}')
            elif wire_shape == 'rectangle':
                if wire_width is not None:
                    out.append(f'width = {wire_width}')
                if wire_height is not None:
                    out.append(f'height = {wire_height}')
            elif wire_shape == 'polygon' and wire_polygon_verts:
                verts_str = ', '.join(f'[{v[0]}, {v[1]}]' for v in wire_polygon_verts)
                out.append(f'vertices = [{verts_str}]')

        # Regions
        if regions:
            out.append('')
            for reg in regions:
                out.append('[[region]]')
                out.append(f'material = {to_toml_string(reg["material"])}')
                out.append(f'inner = {reg["inner"]}')
                out.append(f'outer = {reg["outer"]}')
                out.append('')

    # Landau section
    if conf_str == 'landau':
        out.append('[landau]')
        if landau_nx is not None:
            out.append(f'nx = {landau_nx}')
        if landau_width is not None:
            out.append(f'width = {landau_width}')
        if landau_sweep is not None:
            out.append(f'sweep = {to_toml_string(landau_sweep)}')

    # External field
    if ef_enabled is not None and ef_enabled != 0:
        out.append('')
        out.append('[external_field]')
        if ef_type is not None:
            out.append(f'type = {to_toml_string(ef_type)}')
        if ef_value is not None:
            out.append(f'value = {ef_value}')

    # B field
    if b_field_components is not None:
        out.append('')
        out.append('[b_field]')
        comps = ', '.join(b_field_components)
        out.append(f'components = [{comps}]')
        if b_field_g_factor is not None:
            out.append(f'g_factor = {b_field_g_factor}')

    # B sweep is handled in the BdG section output below

    # Strain
    if strain_enabled is not None and strain_enabled.upper() in ('T', 'TRUE', '1'):
        out.append('')
        out.append('[strain]')
        if strain_ref is not None:
            out.append(f'reference = {to_toml_string(strain_ref)}')
        if strain_solver is not None:
            out.append(f'solver = {to_toml_string(strain_solver)}')
        if piezo is not None:
            out.append(f'piezoelectric = {to_toml_value(piezo)}')
        if strain_substrate_val is not None:
            # numeric substrate value (used for bulk strained)
            out.append(f'strain_substrate = {strain_substrate_val}')

    # SC section
    sc_enabled_flag = False
    if sc_enabled is not None and sc_enabled.upper() in ('T', 'TRUE', '1'):
        sc_enabled_flag = True
    if sc_enabled is not None and sc_enabled not in ('0', 'F', 'f', 'False', 'false'):
        sc_enabled_flag = True
    if sc_enabled_flag or sc_max_iter is not None or sc_tolerance is not None:
        out.append('')
        out.append('[sc]')
        if sc_max_iter is not None:
            out.append(f'max_iterations = {sc_max_iter}')
        if sc_tolerance is not None:
            out.append(f'tolerance = {sc_tolerance}')
        if sc_mixing_alpha is not None:
            out.append(f'mixing_alpha = {sc_mixing_alpha}')
        if sc_diis is not None:
            out.append(f'diis_history = {sc_diis}')
        if sc_temperature is not None:
            out.append(f'temperature = {sc_temperature}')
        if sc_fermi_mode is not None:
            fm = sc_fermi_mode
            # Could be integer or string
            try:
                fm_int = int(fm)
                fm_str = 'charge_neutrality' if fm_int == 0 else 'fixed'
            except ValueError:
                fm_str = fm
            out.append(f'fermi_mode = {to_toml_string(fm_str)}')
        if sc_fermi_level is not None:
            out.append(f'fermi_level = {sc_fermi_level}')
        if sc_num_kpar is not None:
            out.append(f'num_kpar = {sc_num_kpar}')
        if sc_kpar_max is not None:
            out.append(f'kpar_max = {sc_kpar_max}')
        if sc_bc_type is not None:
            out.append(f'bc_type = {to_toml_string(sc_bc_type)}')
        if sc_bc_left is not None:
            out.append(f'bc_left = {sc_bc_left}')
        if sc_bc_right is not None:
            out.append(f'bc_right = {sc_bc_right}')

    # Doping
    for nd, na in uniform_doping:
        out.append('')
        out.append('[[doping]]')
        out.append(f'ND = {nd}')
        out.append(f'NA = {na}')

    for ns, fwhm, pos in delta_doping:
        out.append('')
        out.append('[[doping]]')
        out.append(f'type = "delta"')
        out.append(f'NS = {ns}')
        out.append(f'fwhm = {fwhm}')
        out.append(f'pos = {pos}')

    # Topology
    topo_enabled_flag = False
    if topo_enabled is not None and topo_enabled.upper() in ('T', 'TRUE', '1'):
        topo_enabled_flag = True
    if topo_enabled_flag or topo_mode is not None:
        out.append('')
        out.append('[topology]')
        if topo_mode is not None:
            out.append(f'mode = {to_toml_string(topo_mode)}')
        if topo_compute_chern is not None:
            out.append(f'compute_chern = {to_toml_value(topo_compute_chern)}')
        if topo_compute_hall is not None:
            out.append(f'compute_hall = {to_toml_value(topo_compute_hall)}')
        if topo_qwz_u is not None:
            out.append(f'qwz_u = {topo_qwz_u}')
        if topo_compute_z2 is not None:
            out.append(f'compute_z2 = {to_toml_value(topo_compute_z2)}')
        if topo_z2_method is not None:
            out.append(f'z2_method = {to_toml_string(topo_z2_method)}')
        if topo_extract_edge_states is not None:
            out.append(f'extract_edge_states = {to_toml_value(topo_extract_edge_states)}')
        if topo_edge_E_window is not None:
            out.append(f'edge_E_window = {topo_edge_E_window}')
        if topo_compute_ldos is not None:
            out.append(f'compute_ldos = {to_toml_value(topo_compute_ldos)}')
        if topo_ldos_eta is not None:
            out.append(f'ldos_eta = {topo_ldos_eta}')
        if topo_ldos_E_range:
            vals = ', '.join(topo_ldos_E_range)
            out.append(f'ldos_E_range = [{vals}]')
        if topo_ldos_num_E is not None:
            out.append(f'ldos_num_E = {topo_ldos_num_E}')

        # Gap sweep
        if topo_compute_gap_sweep is not None:
            out.append(f'compute_gap_sweep = {to_toml_value(topo_compute_gap_sweep)}')
        if gap_sweep_B is not None:
            out.append(f'gap_sweep_B_min = {gap_sweep_B[0]}')
            if len(gap_sweep_B) > 1:
                out.append(f'gap_sweep_B_max = {gap_sweep_B[1]}')
            if len(gap_sweep_B) > 2:
                out.append(f'gap_sweep_nB = {gap_sweep_B[2]}')
        if gap_sweep_mu is not None:
            out.append(f'gap_sweep_mu_min = {gap_sweep_mu[0]}')
            if len(gap_sweep_mu) > 1:
                out.append(f'gap_sweep_mu_max = {gap_sweep_mu[1]}')
            if len(gap_sweep_mu) > 2:
                out.append(f'gap_sweep_nMu = {gap_sweep_mu[2]}')
        if sweep_model is not None:
            out.append(f'sweep_model = {to_toml_string(sweep_model)}')

        # Conductance
        if topo_compute_conductance is not None:
            out.append(f'compute_conductance = {to_toml_value(topo_compute_conductance)}')
        if conductance_method is not None:
            out.append(f'conductance_method = {to_toml_string(conductance_method)}')
        if berry_nk is not None:
            out.append(f'berry_nk = {berry_nk}')
        if landauer_energy is not None:
            out.append(f'landauer_energy = {landauer_energy}')

        # Spectral
        if topo_compute_spectral is not None:
            out.append(f'compute_spectral = {to_toml_value(topo_compute_spectral)}')
        if spectral_k_grid is not None:
            out.append(f'spectral_k_min = {spectral_k_grid[0]}')
            if len(spectral_k_grid) > 1:
                out.append(f'spectral_k_max = {spectral_k_grid[1]}')
            if len(spectral_k_grid) > 2:
                out.append(f'spectral_nk = {spectral_k_grid[2]}')
        if spectral_E_grid is not None:
            out.append(f'spectral_E_min = {spectral_E_grid[0]}')
            if len(spectral_E_grid) > 1:
                out.append(f'spectral_E_max = {spectral_E_grid[1]}')
            if len(spectral_E_grid) > 2:
                out.append(f'spectral_nE = {spectral_E_grid[2]}')
            if len(spectral_E_grid) > 3:
                out.append(f'spectral_eta = {spectral_E_grid[3]}')

    # BdG
    bdg_enabled_flag = False
    if bdg_enabled is not None and bdg_enabled.upper() in ('T', 'TRUE', '1'):
        bdg_enabled_flag = True
    if bdg_enabled_flag or bdg_mu is not None or bdg_delta_0 is not None or b_sweep is not None:
        out.append('')
        out.append('[bdg]')
        if bdg_mu is not None:
            out.append(f'mu = {bdg_mu}')
        if bdg_delta_0 is not None:
            out.append(f'delta_0 = {bdg_delta_0}')
        if bdg_g_factor is not None:
            out.append(f'g_factor = {bdg_g_factor}')
        if bdg_B_vec is not None:
            comps = ', '.join(bdg_B_vec)
            out.append(f'B_vec = [{comps}]')
        if bdg_gauge is not None:
            out.append(f'gauge = {to_toml_string(bdg_gauge)}')
        if bdg_kz is not None:
            out.append(f'kz = {bdg_kz}')
        if b_sweep is not None:
            comps = ', '.join(b_sweep)
            out.append(f'B_sweep = [{comps}]')

    # Optics
    optics_enabled_flag = False
    if optics_enabled is not None and optics_enabled.upper() in ('T', 'TRUE', '1'):
        optics_enabled_flag = True
    if optics_enabled_flag:
        out.append('')
        out.append('[optics]')
        if opt_linewidth_lorentzian is not None:
            out.append(f'linewidth_lorentzian = {opt_linewidth_lorentzian}')
        if opt_linewidth_gaussian is not None:
            out.append(f'linewidth_gaussian = {opt_linewidth_gaussian}')
        if opt_refractive_index is not None:
            out.append(f'refractive_index = {opt_refractive_index}')
        if opt_Emin is not None:
            out.append(f'E_min = {opt_Emin}')
        if opt_Emax is not None:
            out.append(f'E_max = {opt_Emax}')
        if opt_NEnergyPoints is not None:
            out.append(f'num_energy_points = {opt_NEnergyPoints}')
        if opt_temperature is not None:
            out.append(f'temperature = {opt_temperature}')
        if opt_carrier_density is not None:
            out.append(f'carrier_density = {opt_carrier_density}')
        if opt_gain is not None:
            out.append(f'gain_enabled = {to_toml_value(opt_gain)}')
        if opt_gain_carrier_density is not None:
            out.append(f'gain_carrier_density = {opt_gain_carrier_density}')
        if opt_ISBT is not None:
            out.append(f'ISBT = {to_toml_value(opt_ISBT)}')
        if opt_spontaneous is not None:
            out.append(f'spontaneous = {to_toml_value(opt_spontaneous)}')
        if opt_spin_resolved is not None:
            out.append(f'spin_resolved = {to_toml_value(opt_spin_resolved)}')

    # Exciton
    exciton_enabled_flag = False
    if exciton_enabled is not None and exciton_enabled.upper() in ('T', 'TRUE', '1'):
        exciton_enabled_flag = True
    if exciton_enabled_flag:
        out.append('')
        out.append('[exciton]')
        if exciton_method is not None:
            out.append(f'method = {to_toml_string(exciton_method)}')

    # Scattering
    scattering_enabled_flag = False
    if scattering_enabled is not None and scattering_enabled.upper() in ('T', 'TRUE', '1'):
        scattering_enabled_flag = True
    if scattering_enabled_flag:
        out.append('')
        out.append('[scattering]')
        if phonon_energy is not None:
            out.append(f'phonon_energy = {phonon_energy}')
        if eps_inf is not None:
            out.append(f'eps_inf = {eps_inf}')
        if eps_0 is not None:
            out.append(f'eps_0 = {eps_0}')

    # Solver (migrated from legacy FEAST)
    if feast_emin is not None or feast_emax is not None or feast_m0 is not None:
        out.append('')
        out.append('[solver]')
        # m0 = -1 means dense fallback
        if feast_m0 is not None and feast_m0 == '-1':
            out.append('method = "DENSE"')
        else:
            out.append('method = "FEAST"')
        out.append('mode = "ENERGY"')
        if feast_emin is not None:
            out.append(f'emin = {feast_emin}')
        if feast_emax is not None:
            out.append(f'emax = {feast_emax}')
        if feast_m0 is not None and feast_m0 != '-1':
            out.append(f'm0 = {feast_m0}')

    # Clean up trailing blank lines
    while out and out[-1] == '':
        out.pop()
    out.append('')  # single trailing newline

    with open(toml_path, 'w') as f:
        f.write('\n'.join(out))

    return toml_path


def main():
    if len(sys.argv) < 2:
        print("Usage: convert_cfg_to_toml.py <config.cfg> [output.toml]")
        print("       convert_cfg_to_toml.py --batch <directory>")
        sys.exit(1)

    if sys.argv[1] == '--batch':
        if len(sys.argv) < 3:
            print("Usage: convert_cfg_to_toml.py --batch <directory>")
            sys.exit(1)
        directory = sys.argv[2]
        count = 0
        for cfg_file in sorted(Path(directory).glob('*.cfg')):
            toml_file = cfg_file.with_suffix('.toml')
            convert_cfg_to_toml(str(cfg_file), str(toml_file))
            print(f"Converted: {cfg_file.name} -> {toml_file.name}")
            count += 1
        print(f"\nConverted {count} files.")
    else:
        cfg_path = sys.argv[1]
        toml_path = sys.argv[2] if len(sys.argv) > 2 else None
        result = convert_cfg_to_toml(cfg_path, toml_path)
        print(f"Converted: {cfg_path} -> {result}")


if __name__ == '__main__':
    main()
