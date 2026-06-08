# Input Parameter Reference (`input.toml`)

Complete reference for every parameter accepted by the `input.toml` file, parsed by `src/io/input_parser.f90` using the `toml-f` library. Parameters are organized into TOML sections; sections are order-independent within the file.

**Config format:** Standard TOML (Tom's Obvious Minimal Language). Lines beginning with `#` are comments. Section presence enables optional physics blocks (no separate enable flags).

**Modes:** B = bulk (`confinement = "bulk"`), Q = quantum well (`confinement = "qw"`), W = wire (`confinement = "wire"`), L = Landau (`confinement = "landau"`), G = g-factor, S = self-consistent.

See `tests/regression/configs/` for canonical examples.

---

## 1. Top-Level Parameters

Required parameters that control the simulation mode and discretization.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `confinement` | string | (required) | `"bulk"`, `"qw"`, `"wire"`, `"landau"` | all | Simulation mode: `"bulk"` = 8x8 Hamiltonian, `"qw"` = 8NxN (1D confinement along z), `"wire"` = 2D confinement in x-y, `"landau"` = Landau levels with magnetic field. |
| `FDorder` | integer | `2` | 2, 4, 6, 8, 10 | Q W | Finite-difference stencil order. Higher orders give better accuracy but require more grid points. Grid must have >= `FDorder + 1` points. |
| `fd_step` | integer | `1` | >= 1 (bulk), >= 3 (QW) | B Q | Number of finite-difference grid points along the confinement direction. For bulk, set to any value (internally forced to 1). For wire/Landau, this key is not used; the grid is derived from the `[wire]` or `[landau]` section. |
| `which_band` | integer | `0` | 0, 1 | G | Band edge for g-factor calculation: 0 = conduction band, 1 = valence band. |
| `band_idx` | integer | `1` | >= 1 | G | Subband doublet index: 1 = ground state, 2 = first excited, etc. |

---

## 2. Wave Vector (`[wave_vector]`)

Controls the k-space sweep direction and range.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `mode` | string | `"k0"` | `"k0"`, `"kx"`, `"ky"`, `"kz"`, `"kxky"`, `"kpar"` | B Q W G S | Wave-vector direction for the k-sweep. `"k0"` means single-point calculation at Gamma. |
| `max` | float | `0.0` | >= 0 (1/A) | B Q W G S | Maximum wave-vector magnitude in inverse Angstroms. 0 for single-point. |
| `nsteps` | integer | `100` | >= 0 | B Q W G S | Number of steps in the k-sweep. 0 or 1 for single-point. |
| `step` | float | `0.01` | > 0 (1/A) | B Q W G S | Wave-vector step size in inverse Angstroms. |

Example:
```toml
[wave_vector]
mode = "kx"
max = 0.1
nsteps = 21
```

---

## 3. Band Counts (`[bands]`)

Number of eigenvalues to compute in each band group.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `num_cb` | integer | `2` | >= 1 | B Q W G S | Number of conduction-band eigenvalues to compute. Must be even (each state is a spin doublet). Bulk: max 2. QW/wire: depends on system size. |
| `num_vb` | integer | `6` | >= 1 | B Q W G S | Number of valence-band eigenvalues to compute. Must be even. Bulk: max 6. QW/wire: depends on system size. |

The total eigenvalue count is `evnum = num_cb + num_vb`. For bulk, the Hamiltonian is 8x8 so the maximum is `num_cb=2`, `num_vb=6`.

Example:
```toml
[bands]
num_cb = 4
num_vb = 8
```

---

## 4. Material Layers (`[[material]]`) -- Bulk and QW

Material layers for bulk and QW simulations. Each layer is a TOML array-of-tables entry. The number of layers is inferred from the count of `[[material]]` entries (no explicit `numLayers` field).

For **bulk**, exactly one `[[material]]` entry is required (no `z_min`/`z_max`).

For **QW**, two or more entries define the layered structure.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `name` | string | (required) | any valid material name | B Q | Material identifier from the parameter database (e.g. `GaAs`, `InAsW`, `Al30Ga70As`). |
| `z_min` | float | (required for QW) | any (A) | Q | Layer lower boundary in Angstroms. |
| `z_max` | float | (required for QW) | any (A) | Q | Layer upper boundary in Angstroms. |

**Last-layer-wins** convention: later material entries overwrite earlier ones at overlapping z-positions. Use the 2-layer pattern (barrier covers full domain, well overwrites center).

Example -- bulk:
```toml
[[material]]
name = "GaAs"
```

Example -- QW (GaAs well in AlGaAs barriers):
```toml
[[material]]
name = "Al30Ga70As"
z_min = -200
z_max = 200

[[material]]
name = "GaAs"
z_min = -50
z_max = 50
```

Example -- broken-gap InAs/GaSb (3-layer):
```toml
[[material]]
name = "AlSbW"
z_min = -250
z_max = 250

[[material]]
name = "GaSbW"
z_min = -135
z_max = 135

[[material]]
name = "InAsW"
z_min = -35
z_max = 35
```

---

## 5. Wire Geometry (`[wire]`) -- confinement = "wire"

Wire cross-section parameters. Only used when `confinement = "wire"`.

### `[wire]` grid parameters

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `nx` | integer | (required) | >= 3 | W | Number of grid points in the x direction. |
| `ny` | integer | (required) | >= 3 | W | Number of grid points in the y direction. |
| `dx` | float | (required) | > 0 (A) | W | Grid spacing in x (Angstroms). |
| `dy` | float | (required) | > 0 (A) | W | Grid spacing in y (Angstroms). |

### `[wire.geometry]` cross-section shape

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `shape` | string | (required) | `"rectangle"`, `"circle"`, `"hexagon"`, `"polygon"` | W | Cross-section shape. Determines which dimensional parameters are needed. |
| `width` | float | -- | > 0 (A) | W | Full width in x for `"rectangle"` shape. |
| `height` | float | -- | > 0 (A) | W | Full height in y for `"rectangle"` shape. |
| `radius` | float | -- | > 0 (A) | W | Radius for `"circle"` and `"hexagon"` shapes. |
| `vertices` | array of [x,y] pairs | -- | >= 3 vertices | W | Polygon vertices for `"polygon"` shape. Each vertex is `[x, y]` in Angstroms. |

### `[[region]]` material regions

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `material` | string | (required) | any valid material name | W | Material for this region. |
| `inner` | float | (required) | >= 0 (A) | W | Inner radius/distance of the region in Angstroms. |
| `outer` | float | (required) | > `inner` (A) | W | Outer radius/distance of the region in Angstroms. |

Example -- rectangular wire with single region:
```toml
[wire]
nx = 21
ny = 21
dx = 3.0
dy = 3.0

[wire.geometry]
shape = "rectangle"
width = 63.0
height = 63.0

[[region]]
material = "GaAs"
inner = 0.0
outer = 100.0
```

Example -- circular core-shell wire:
```toml
[wire]
nx = 41
ny = 41
dx = 1.0
dy = 1.0

[wire.geometry]
shape = "circle"
radius = 20.0

[[region]]
material = "InAs"
inner = 0.0
outer = 10.0

[[region]]
material = "GaAs"
inner = 10.0
outer = 20.0
```

Example -- hexagonal wire:
```toml
[wire]
nx = 31
ny = 31
dx = 2.0
dy = 2.0

[wire.geometry]
shape = "hexagon"
radius = 30.0

[[region]]
material = "GaAs"
inner = 0.0
outer = 100.0
```

---

## 6. Landau Levels (`[landau]`) -- confinement = "landau"

Landau level calculation parameters. Only used when `confinement = "landau"`.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `nx` | integer | `100` | >= 3 | L | Number of grid points in the confinement direction. |
| `width` | float | `2000.0` | > 0 (A) | L | Domain width in Angstroms. |
| `sweep` | string | `"ky"` | `"ky"`, `"kz"`, `"B"` | L | Sweep variable: `"ky"` = in-plane k, `"kz"` = out-of-plane k, `"B"` = magnetic field magnitude. |

Example:
```toml
[landau]
nx = 100
width = 2000.0
sweep = "B"
```

---

## 7. External Field (`[external_field]`)

Optional. Electric field applied to the structure. If the section is absent, no external field is applied.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `type` | string | `"EF"` | `"EF"` | B Q W G S | Field type. Currently only `"EF"` (electric field) is supported. |
| `value` | float | `0.0` | any (V/A) | B Q W G S | Electric field strength in V/Angstrom. |

Example:
```toml
[external_field]
type = "EF"
value = 0.070
```

---

## 8. Magnetic Field (`[b_field]`)

Optional. Magnetic field for Landau levels and Zeeman splitting. If the section is absent, B = 0.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `components` | 3 floats | `[0.0, 0.0, 0.0]` | any (T) | L W | Magnetic field components `[Bx, By, Bz]` in Tesla. |
| `g_factor` | float | `2.0` | > 0 | L | Lande g-factor for Landau level calculation. |

Example:
```toml
[b_field]
components = [0.0, 0.0, 5.0]
g_factor = 2.0
```

---

## 9. Strain (`[strain]`)

Optional. Enables strain calculation. If the section is absent, strain is not computed.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `reference` | string | -- | material name | Q W | Reference material for strain. Uses the lattice constant of this material as the unstrained reference. |
| `strain_substrate` | float | -- | > 0 (A) | B | Bulk strain: explicit substrate lattice constant in Angstroms. Must be specified inside the `[strain]` section. |
| `solver` | string | `"pardiso"` | `"pardiso"` | Q W | Strain solver backend. Currently only PARDISO direct solver is supported. |
| `piezoelectric` | boolean | `false` | `true`/`false` | Q W | Include piezoelectric polarization in the strain calculation. |

Example -- QW strain with substrate reference:
```toml
[strain]
reference = "GaAs"
```

Example -- bulk strain with explicit lattice constant:
```toml
[strain]
strain_substrate = 5.869
```

---

## 10. Self-Consistent Schrodinger-Poisson (`[sc]`)

Optional. Enables the self-consistent loop. If the section is absent, no SC calculation is performed.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `max_iterations` | integer | `100` | >= 1 | Q S | Maximum number of SC iterations. |
| `tolerance` | float | `1.0e-6` | > 0 (eV) | S | Convergence threshold: maximum absolute change in electrostatic potential (infinity norm) in eV. |
| `mixing_alpha` | float | `0.3` | (0, 1] | S | Linear mixing parameter for the initial warm-up phase before DIIS kicks in. Smaller values are more stable but slower. |
| `diis_history` | integer | `7` | >= 1 | S | Number of previous iterations stored for DIIS/Pulay acceleration. |
| `temperature` | float | `300.0` | > 0 (K) | S | Temperature for Fermi-Dirac statistics. |
| `fermi_mode` | string | `"charge_neutrality"` | `"charge_neutrality"`, `"fixed"` | S | Fermi level determination: `"charge_neutrality"` finds the Fermi level automatically, `"fixed"` uses the `fermi_level` value. |
| `fermi_level` | float | `0.0` | any (eV) | S | Fixed Fermi level in eV. Only used when `fermi_mode = "fixed"`. |
| `num_kpar` | integer | `201` | >= 1, odd preferred | S | Number of in-plane wave-vector sampling points for charge density integration. Odd numbers enable Simpson's rule. |
| `kpar_max` | float | `0.0` | >= 0 (1/A) | S | Maximum in-plane wave-vector. 0 = automatic determination. |
| `bc_type` | string | `"DD"` | `"DD"`, `"DN"`, `"ND"`, `"NN"` | S | Poisson boundary condition type: D = Dirichlet, N = Neumann. |
| `bc_left` | float | `0.0` | any (eV) | S | Left boundary potential in eV (Dirichlet). |
| `bc_right` | float | `0.0` | any (eV) | S | Right boundary potential in eV (Dirichlet for `"DD"`). |

Example:
```toml
[sc]
max_iterations = 50
tolerance = 1.0e-6
mixing_alpha = 0.3
diis_history = 7
temperature = 300.0
fermi_mode = "charge_neutrality"
num_kpar = 41
kpar_max = 0.2
bc_type = "DD"
bc_left = 0.0
bc_right = 0.0
```

---

## 11. Doping (`[[doping]]`)

Optional. Per-layer doping specifications for self-consistent calculations. Each entry is a TOML array-of-tables. Used with `[sc]`.

### Uniform doping

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `ND` | float | `0.0` | >= 0 (cm^-3) | S | Donor concentration (n-type). |
| `NA` | float | `0.0` | >= 0 (cm^-3) | S | Acceptor concentration (p-type). |

### Delta doping

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `type` | string | `"uniform"` | `"uniform"`, `"delta"` | S | Doping profile type. `"delta"` activates delta-function doping. |
| `NS` | float | -- | >= 0 (cm^-2) | S | Sheet doping density for delta doping. |
| `fwhm` | float | -- | > 0 (A) | S | Full width at half maximum of the Gaussian smearing for delta doping. |
| `pos` | float | -- | any (A) | S | Position of the delta doping peak in Angstroms. |

Example -- uniform per-layer doping:
```toml
[[doping]]
ND = 0.0
NA = 0.0

[[doping]]
ND = 5.0e18
NA = 0.0
```

Example -- delta doping:
```toml
[[doping]]
type = "delta"
NS = 5.0
fwhm = 10.0
pos = 0.0
```

---

## 12. Topological Analysis (`[topology]`)

Optional. Enables topological invariant computation. If the section is absent, no topological analysis is performed.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `mode` | string | `"qhe"` | `"qhe"`, `"qshe"`, `"bdg"`, `"spectral"`, `"conductance"`, `"sweep"` | all | Topological mode: `"qhe"` = Chern number, `"qshe"` = Z2 invariant, `"bdg"` = Majorana modes, `"spectral"` = spectral function, `"conductance"` = conductance, `"sweep"` = phase diagram. |
| `compute_chern` | boolean | `false` | `true`/`false` | all | Compute Chern number via Berry curvature integration. |
| `compute_hall` | boolean | `false` | `true`/`false` | all | Compute Hall conductance. |
| `qwz_u` | float | `0.0` | any (eV) | QHE | QWZ model mass parameter for Chern number computation. |
| `bhz_M` | float | `10.0` | any (eV) | all | BHZ model mass parameter for Z2/gap sweep computation. |
| `compute_z2` | boolean | `false` | `true`/`false` | all | Compute Z2 invariant. |
| `z2_method` | string | `"auto"` | `"auto"`, `"gap"`, `"parity"` | all | Z2 computation method. `"auto"` selects based on geometry. |
| `extract_edge_states` | boolean | `false` | `true`/`false` | W | Extract edge state energies and localization length. |
| `edge_E_window` | float | `0.01` | > 0 (eV) | W | Energy window for edge state detection. |
| `compute_ldos` | boolean | `false` | `true`/`false` | W | Compute local density of states via complex PARDISO. |
| `ldos_eta` | float | `0.001` | > 0 (eV) | W | Lorentzian broadening for LDOS. |
| `ldos_E_range` | 2 floats | `[-0.1, 0.1]` | any (eV) | W | Energy range for LDOS computation. |
| `ldos_num_E` | integer | `200` | >= 1 | W | Number of energy points for LDOS. |

### Gap sweep sub-section

Nested under `[topology]` when `mode = "sweep"`:

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `compute_gap_sweep` | boolean | `false` | `true`/`false` | B/Q/W | Enable gap sweep phase diagram. |
| `gap_sweep_B_min` | float | `0.0` | >= 0 (T) | W/Q | Minimum magnetic field for gap sweep. |
| `gap_sweep_B_max` | float | `1.0` | >= `gap_sweep_B_min` (T) | W/Q | Maximum magnetic field for gap sweep. |
| `gap_sweep_nB` | integer | `20` | >= 1 | W/Q | Number of B-field points in gap sweep. |
| `gap_sweep_mu` | 3 values | `[0.0, 0.01, 20]` | `[min, max, npts]` (eV) | W/Q | Chemical-potential sweep grid. |
| `sweep_model` | string | `"bhz_analytic"` | `"bhz_analytic"`, `"wire_bdg"`, `"qw_fukane"` | B/Q/W | Gap sweep evaluator. |

### Spectral sub-section

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `compute_spectral` | boolean | `false` | `true`/`false` | B/Q/W | Compute spectral function. |
| `spectral_k_grid` | 3 values | `[-0.1, 0.1, 100]` | `[k_min, k_max, nk]` (1/A) | B/Q/W | Wave-vector grid for spectral function. |
| `spectral_E_grid` | 4 values | `[-0.05, 0.05, 200, 0.001]` | `[E_min, E_max, nE, eta]` (eV) | B/Q/W | Energy grid and Lorentzian width. |

### Conductance sub-section

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `compute_conductance` | boolean | `false` | `true`/`false` | B/Q/W | Compute conductance. |
| `conductance_method` | string | `"kubo_chern"` | `"kubo_chern"`, `"kubo_berry"`, `"landauer"` | B/Q/W | Conductance method. |
| `berry_nk` | integer | `50` | >= 2 | B/Q | Kubo Berry grid size. |
| `landauer_energy` | float | `0.0` | any (eV) | W | Energy for Landauer helper. |

### Output files

| File | Condition | Columns | Description |
|---|---|---|---|
| `topology_result.dat` | `[topology]` present | header + data | Topological invariants, conductance, min gap, Majorana count. |
| `spectral_function.dat` | `mode = "spectral"` | `k(1/AA) E(eV) A(k,E)` | Spectral function heatmap. |
| `z2_phase_diagram.dat` | `mode = "sweep"` | `B(T) mu(eV) z2 gap(eV)` | Z2/gap sweep results. |
| `z2_transitions.dat` | `mode = "sweep"` | `B(T) mu(eV)` | Detected transition midpoints. |

Example -- QHE Chern number:
```toml
[topology]
mode = "qhe"
compute_chern = true
compute_hall = true
```

Example -- BdG Majorana modes:
```toml
[topology]
mode = "bdg"
extract_edge_states = true
edge_E_window = 0.01
```

---

## 13. BdG Parameters (`[bdg]`)

Optional. Bogoliubov-de Gennes parameters for topological superconductor / Majorana analysis. Used with `topology.mode = "bdg"`.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `mu` | float | `0.0` | any (eV) | W | Chemical potential in eV. |
| `delta_0` | float | `0.0` | >= 0 (eV) | W | s-wave superconducting gap amplitude in eV. |
| `g_factor` | float | `2.0` | > 0 | W | Lande g-factor for Zeeman splitting. |
| `B_vec` | 3 floats | `[0.0, 0.0, 0.0]` | any (T) | W | Magnetic field `[Bx, By, Bz]` for Zeeman splitting in Tesla. |
| `gauge` | string | `"landau_x"` | `"landau_x"`, `"landau_z"`, `"zeeman"` | W | Gauge choice for magnetic field coupling. |
| `kz` | float | `0.0` | any (1/A) | W | Out-of-plane wave vector. |
| `self_consistent` | boolean | `false` | `true`/`false` | W | Enable self-consistent gap computation (future). |
| `B_sweep` | 3 floats | -- | `[min, max, step]` (T) | W | B-field sweep parameters for phase diagram. |

Example:
```toml
[bdg]
mu = 0.0
delta_0 = 0.001
g_factor = 2.0
B_vec = [0.0, 0.0, 0.0]
gauge = "landau_x"
B_sweep = [0.5, 10.0, 0.5]
```

---

## 14. Optical Spectra (`[optics]`)

Optional. Enables optical spectra calculation. If the section is absent, no optical calculation is performed.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `linewidth_lorentzian` | float | `0.030` | >= 0 (eV) | Q W | Lorentzian FWHM linewidth for homogeneous broadening. |
| `linewidth_gaussian` | float | `0.005` | >= 0 (eV) | Q W | Gaussian FWHM linewidth for inhomogeneous broadening. |
| `refractive_index` | float | `3.3` | > 0 | Q W | Background refractive index used in the prefactor. |
| `E_min` | float | `0.5` | > 0 (eV) | Q W | Minimum photon energy in the spectral grid. |
| `E_max` | float | `2.0` | > `E_min` (eV) | Q W | Maximum photon energy in the spectral grid. |
| `num_energy_points` | integer | `200` | >= 1 | Q W | Number of energy grid points. |
| `temperature` | float | `300.0` | > 0 (K) | Q W | Temperature for Fermi-Dirac occupation factors. |
| `carrier_density` | float | `0.0` | >= 0 (cm^-2) | Q W | 2D carrier density for equilibrium absorption. 0 means intrinsic. |
| `gain_enabled` | boolean | `false` | `true`/`false` | Q W | Enable gain spectrum calculation. |
| `gain_carrier_density` | float | `3.0e12` | > 0 (cm^-2) | Q W | 2D carrier density for gain calculation. |
| `ISBT` | boolean | `false` | `true`/`false` | Q W | Enable intersubband transition (ISBT) absorption. |
| `spontaneous` | boolean | `false` | `true`/`false` | Q W | Enable spontaneous emission spectrum. |
| `spin_resolved` | boolean | `false` | `true`/`false` | Q W | Enable spin-up/down decomposition of absorption spectra. |

### Output files

| File | Condition | Columns | Description |
|---|---|---|---|
| `absorption_TE.dat` | `[optics]` present | E(eV), alpha(cm^-1) | Interband TE absorption |
| `absorption_TM.dat` | `[optics]` present | E(eV), alpha(cm^-1) | Interband TM absorption |
| `absorption_ISBT.dat` | `ISBT = true` | E(eV), alpha(cm^-1) | ISBT absorption |
| `gain_TE.dat` | `gain_enabled = true` | E(eV), gain(cm^-1) | Interband TE gain |
| `gain_TM.dat` | `gain_enabled = true` | E(eV), gain(cm^-1) | Interband TM gain |
| `spontaneous_TE.dat` | `spontaneous = true` | E(eV), rate(cm^-1) | Spontaneous emission TE |
| `spontaneous_TM.dat` | `spontaneous = true` | E(eV), rate(cm^-1) | Spontaneous emission TM |
| `absorption_TE_up.dat` | `spin_resolved = true` | E(eV), alpha(cm^-1) | Spin-up TE absorption |
| `absorption_TE_dw.dat` | `spin_resolved = true` | E(eV), alpha(cm^-1) | Spin-down TE absorption |

Example:
```toml
[optics]
linewidth_lorentzian = 0.030
linewidth_gaussian = 0.005
refractive_index = 3.3
E_min = 1.3
E_max = 2.0
num_energy_points = 300
temperature = 300.0
carrier_density = 0.0
gain_enabled = false
ISBT = false
spontaneous = false
spin_resolved = false
```

---

## 15. Exciton (`[exciton]`)

Optional. Enables exciton binding energy calculation.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `method` | string | `"variational"` | `"variational"` | Q | Exciton solver method. |

Example:
```toml
[exciton]
method = "variational"
```

---

## 16. Scattering (`[scattering]`)

Optional. Phonon scattering parameters.

| Key | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `phonon_energy` | float | -- | > 0 (eV) | Q | LO phonon energy in eV. |
| `eps_inf` | float | -- | > 0 | Q | High-frequency dielectric constant. |
| `eps_0` | float | -- | > 0 | Q | Static dielectric constant. |

Example:
```toml
[scattering]
phonon_energy = 0.036
eps_inf = 10.9
eps_0 = 12.9
```

---

## 17. Eigensolver (`[solver]`)

Optional. Configures the eigensolver method and parameters for all confinement
modes (bulk, QW, wire, Landau). When omitted, the solver uses smart defaults
based on the confinement mode.

| Key | Type | Default | Valid range | Description |
|---|---|---|---|---|
| `method` | string | `"AUTO"` | `AUTO`, `DENSE`, `FEAST` | Eigensolver algorithm. `AUTO` selects based on confinement mode and matrix size. |
| `mode` | string | `"AUTO"` | `AUTO`, `FULL`, `INDEX`, `ENERGY` | Eigenvalue selection mode. `FULL` = all eigenvalues, `INDEX` = range il:iu, `ENERGY` = range [emin, emax]. |
| `emin` | float | `0.0` | any (eV) | Lower bound of the energy search interval (ENERGY mode). 0 = automatic. |
| `emax` | float | `0.0` | any (eV) | Upper bound of the energy search interval (ENERGY mode). 0 = automatic. |
| `m0` | integer | `0` | >= 0 | Initial subspace dimension for FEAST. 0 = automatic (2 * number of eigenvalues sought). |

### Smart defaults (AUTO)

When `method = "AUTO"` (the default), the code selects:

| Confinement | Matrix size | Method chosen |
|---|---|---|
| bulk | 8 x 8 | `DENSE` (zheevx) |
| QW | <= threshold | `DENSE` (zheevx) |
| QW | > threshold | `FEAST` (sparse CSR) |
| wire | any | `FEAST` (sparse CSR) |
| Landau | any | `DENSE` (zheevx) |

When `mode = "AUTO"` (the default), the code selects:

| Confinement | Mode chosen |
|---|---|
| bulk | `FULL` |
| QW | `ENERGY` (uses emin/emax or auto-computed window) |
| wire | `ENERGY` (uses emin/emax or auto-computed window) |
| Landau | `FULL` |

### Valid combinations

| method | mode | Notes |
|---|---|---|
| `DENSE` | `FULL` | Bulk/Landau default; all eigenvalues via zheevx |
| `DENSE` | `ENERGY` | Dense solve with energy window filtering |
| `FEAST` | `ENERGY` | Sparse FEAST with contour [emin, emax]; requires emin < emax |
| `FEAST` | `FULL` | FEAST with auto-computed full energy window |
| `FEAST` | `INDEX` | **Rejected** by validation (FEAST does not support index-based selection) |

### Examples

FEAST with explicit subspace dimension:
```toml
[solver]
method = "FEAST"
mode = "ENERGY"
emin = -1.5
emax = 2.0
m0 = 128
```

Dense LAPACK with energy window:
```toml
[solver]
method = "DENSE"
mode = "ENERGY"
emin = -1.5
emax = 2.0
```

Minimal (uses smart defaults):
```toml
[solver]
emin = -1.5
emax = 2.0
```

### Migration from `[feast]`

The legacy `[feast]` section is no longer supported. Migrate as follows:

| Old `[feast]` | New `[solver]` |
|---|---|
| `m0 = -1` | `method = "DENSE"`, `mode = "ENERGY"` |
| `m0 = 0` or absent | `method = "FEAST"`, `mode = "ENERGY"` |
| `m0 > 0` | `method = "FEAST"`, `mode = "ENERGY"`, keep `m0` |

---

## Complete Examples

### Bulk GaAs, k along x

```toml
confinement = "bulk"
FDorder = 2
fd_step = 101

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 11

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "GaAs"
```

### Quantum well (GaAs/AlGaAs), k-sweep

```toml
confinement = "qw"
FDorder = 2
fd_step = 201

[wave_vector]
mode = "kx"
max = 0.1
nsteps = 21

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "Al30Ga70As"
z_min = -200
z_max = 200

[[material]]
name = "GaAs"
z_min = -50
z_max = 50
```

### g-Factor (bulk GaAs conduction band)

```toml
confinement = "bulk"
FDorder = 2
fd_step = 1
which_band = 0
band_idx = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "GaAs"
```

### Quantum wire (GaAs rectangular cross section)

```toml
confinement = "wire"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "kz"
max = 0.1
nsteps = 3

[bands]
num_cb = 8
num_vb = 16

[wire]
nx = 21
ny = 21
dx = 3.0
dy = 3.0

[wire.geometry]
shape = "rectangle"
width = 63.0
height = 63.0

[[region]]
material = "GaAs"
inner = 0.0
outer = 100.0

[solver]
method = "FEAST"
mode = "ENERGY"
emin = -1.5
emax = 2.0
m0 = 128
```

### Landau levels (GaAs, B = 5 T)

```toml
confinement = "landau"
FDorder = 2
fd_step = 100

[wave_vector]
mode = "ky"
max = 0.05
nsteps = 10

[bands]
num_cb = 4
num_vb = 4

[[material]]
name = "GaAs"

[landau]
nx = 100
width = 2000.0

[b_field]
components = [0.0, 0.0, 5.0]
```

### Self-consistent SP (GaAs/AlAs QW, n-doped well)

```toml
confinement = "qw"
FDorder = 2
fd_step = 101
which_band = 0
band_idx = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "AlAs"
z_min = -150
z_max = 150

[[material]]
name = "GaAs"
z_min = -50
z_max = 50

[sc]
max_iterations = 50
tolerance = 1.0e-6
mixing_alpha = 0.3
diis_history = 7
temperature = 300.0
fermi_mode = "charge_neutrality"
num_kpar = 41
kpar_max = 0.2
bc_type = "DD"
bc_left = 0.0
bc_right = 0.0

[[doping]]
ND = 0.0
NA = 0.0

[[doping]]
ND = 5.0e18
NA = 0.0
```

### Optical properties (GaAs/AlGaAs QW absorption)

```toml
confinement = "qw"
FDorder = 4
fd_step = 101
which_band = 0
band_idx = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 0

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "Al30Ga70As"
z_min = -200
z_max = 200

[[material]]
name = "GaAs"
z_min = -50
z_max = 50

[optics]
linewidth_lorentzian = 0.030
linewidth_gaussian = 0.005
refractive_index = 3.3
E_min = 1.3
E_max = 2.0
num_energy_points = 300
temperature = 300.0
carrier_density = 0.0
gain_enabled = false
gain_carrier_density = 0.0
ISBT = false
spontaneous = false
spin_resolved = false
```

### Strained QW (InAs/GaAs)

```toml
confinement = "qw"
FDorder = 2
fd_step = 201

[wave_vector]
mode = "k0"
max = 0.1
nsteps = 1

[bands]
num_cb = 4
num_vb = 8

[[material]]
name = "GaAs"
z_min = -60
z_max = 60

[[material]]
name = "InAs"
z_min = -10
z_max = 10

[strain]
reference = "GaAs"
```

### Topological analysis (QHE Chern number)

```toml
confinement = "bulk"
FDorder = 2
fd_step = 1

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "GaAs"

[topology]
mode = "qhe"
compute_chern = true
compute_hall = true
```

### BdG / Majorana modes

```toml
confinement = "qw"
FDorder = 2
fd_step = 7

[wave_vector]
mode = "k0"
max = 0.0
nsteps = 1

[bands]
num_cb = 2
num_vb = 6

[[material]]
name = "GaAs"
z_min = -30
z_max = 30

[topology]
mode = "bdg"
extract_edge_states = true
edge_E_window = 0.01

[bdg]
mu = 0.0
delta_0 = 0.001
g_factor = 2.0
B_vec = [0.0, 0.0, 0.0]
gauge = "landau_x"
```

---

## Notes

- All TOML keys use `snake_case` convention.
- Sections are order-independent within the file.
- Optional sections are enabled by their presence in the file (no separate enable flags).
- Comments use `#` prefix (standard TOML).
- Material names are case-sensitive and must match the database in `src/core/parameters.f90` (e.g. `GaAs`, `InAsW`, `Al20Ga80As`).
- Alloy materials like `Al20Ga80As` use linear interpolation between endpoint binaries with the indicated composition fraction.
- The `fd_step` field is the user-specified grid size for bulk and QW modes. Internally, the grid system (`cfg%grid`) manages the actual point count, accessible via `cfg%grid%npoints()`: bulk returns 1, QW returns `fd_step`, wire returns `nx * ny`, Landau returns `landau.nx`.
- For bulk mode, the `[[material]]` entry does not need `z_min`/`z_max`.
- For wire mode, `fd_step` is not used; the grid is specified via the `[wire]` section (`nx`, `ny`, `dx`, `dy`).
- The `topologicalAnalysis` executable uses the same `input.toml` format with `[topology]` and `[bdg]` sections.
