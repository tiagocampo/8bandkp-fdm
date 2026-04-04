# Input Parameter Reference (`input.cfg`)

Complete reference for every parameter accepted by the `input.cfg` file, parsed by `src/io/input_parser.f90`. Parameters are read as label-value pairs; lines beginning with `!` are comments.

Parameter order in the file is **positional** (the parser reads sequentially). The groupings below reflect logical purpose, not file order. See the example configs (`bulk.example`, `quantumwell.example`, `gfactor.example`) for correct ordering.

**Modes:** B = bulk (`confinement=0`), Q = quantum well (`confinement=1`), W = wire (`confinement=2`), G = g-factor, S = self-consistent.

---

## 1. General

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `waveVector` | string(2) | `k0` | `k0`, `kx`, `ky`, `kz` | B Q W G S | Wave-vector direction for the k-sweep. `k0` means single-point calculation at Gamma. |
| `waveVectorMax` | float | `0.0` | >= 0 (1/A) | B Q W G S | Maximum wave-vector magnitude. 0 for single-point. |
| `waveVectorStep` | integer | `100` | >= 0 | B Q W G S | Number of steps in the k-sweep. 0 or 1 for single-point. |

---

## 2. Structure

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `confinement` | integer | `0` | 0, 1, 2 | all | Simulation mode: 0 = bulk (8x8), 1 = quantum well (8NxN, 1D confinement along z), 2 = wire (2D confinement in x-y). |
| `FDstep` | integer | `1` | >= 1 (bulk), >= 3 (QW) | B Q | Number of finite-difference grid points along the confinement direction. Bulk mode forces this to 1. Wire mode sets this equal to `wire_ny`. Must be >= `FDorder + 1` for QW. |
| `FDorder` | integer | `2` | 2, 4, 6, 8, 10 | Q W | Finite-difference stencil order. Higher orders give better accuracy but require more grid points and wider stencils. Grid must have >= `FDorder + 1` points. |
| `numLayers` | integer | `1` | >= 1 | B Q W | Number of material layers. Bulk must be 1. For wire mode, overridden to match `numRegions`. |

---

## 3. Structure -- Materials (QW mode)

Each layer is specified on a separate line. The label uses a 1-based index suffix.

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `materialN` | string + 2 floats | (required) | -- | Q | Layer specification: `materialN: <name> <startPos> <endPos>`. `name` is a material identifier from the parameter database (e.g. `GaAs`, `InAsW`). `startPos` and `endPos` are the layer boundaries in Angstroms. First layer's `startPos` to `endPos` defines the total structure extent. |

Example:
```
material1: AlSb -250 250
material2: GaSb -135 135
material3: InAs  -35  35
```

---

## 4. Structure -- Materials (Bulk mode)

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `material1` | string | (required) | any valid material name | B | Single material name for bulk calculation (e.g. `GaAs`, `InAsW`). |

---

## 5. Structure -- Wire Geometry (confinement=2)

These parameters appear after `numLayers` only in wire mode.

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `wire_nx` | integer | `0` | >= 3, >= `FDorder + 1` | W | Number of grid points in the x (confinement) direction. |
| `wire_ny` | integer | `0` | >= 3, >= `FDorder + 1` | W | Number of grid points in the y (confinement) direction. |
| `wire_dx` | float | `0.0` | > 0 (A) | W | Grid spacing in x (Angstroms). |
| `wire_dy` | float | `0.0` | > 0 (A) | W | Grid spacing in y (Angstroms). |
| `wire_shape` | string(16) | `rectangle` | `rectangle`, `circle`, `hexagon`, `polygon` | W | Cross-section shape of the wire. Determines which dimensional parameters follow. |
| `wire_radius` | float | `0.0` | > 0 (A) | W | Radius for `circle` and `hexagon` shapes. |
| `wire_width` | float | `0.0` | > 0 (A) | W | Full width in x for `rectangle` shape. |
| `wire_height` | float | `0.0` | > 0 (A) | W | Full height in y for `rectangle` shape. |
| `wire_polygon` | integer + N pairs | -- | nverts >= 3 | W | For `polygon` shape: first the vertex count, then one line per vertex with x y coordinates. |
| `wire_vertexN` | 2 floats | -- | -- | W | Individual polygon vertex coordinates (x, y in Angstroms). Read `nverts` times after `wire_polygon`. |
| `numRegions` | integer | `0` | >= 1 | W | Number of material regions in the wire cross-section. |
| `region` | string + 2 floats | (required) | -- | W | Region specification: `region: <material> <inner> <outer>`. `inner` and `outer` define the radial (or distance-based) extent of the region in Angstroms. Read `numRegions` times. |

Example (rectangle wire with single region):
```
wire_nx: 11
wire_ny: 11
wire_dx: 2.0
wire_dy: 2.0
wire_shape: rectangle
wire_width: 22.0
wire_height: 22.0
numRegions: 1
region: GaAs  0.0  100.0
```

Example (core-shell circle wire):
```
wire_nx: 41
wire_ny: 41
wire_dx: 1.0
wire_dy: 1.0
wire_shape: circle
wire_radius: 20.0
numRegions: 2
region: InAs  0.0  10.0
region: GaAs  10.0  20.0
```

---

## 6. Eigenvalues

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `numcb` | integer | `2` | >= 1 | B Q W G S | Number of conduction-band eigenvalues to compute. Must be even (each state is a spin doublet). Bulk: max 2. QW/wire: depends on system size. |
| `numvb` | integer | `6` | >= 1 | B Q W G S | Number of valence-band eigenvalues to compute. Must be even. Bulk: max 6. QW/wire: depends on system size. |

The total eigenvalue count is `evnum = numcb + numvb`. For bulk, the Hamiltonian is 8x8 so the maximum is `numcb=2`, `numvb=6`. For QW/wire, the Hamiltonian is 8NxN where N is the number of grid points.

---

## 7. External Field

These two lines always appear together.

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `ExternalField` | integer + string(2) | `0 EF` | 0 or 1; `EF` | B Q W G S | First value enables/disables external field (0 = off, 1 = on). Second value is the field type: currently only `EF` (electric field) is supported. |
| `EFParams` | float | `0.0005` | any (V/A) | B Q W G S | Electric field strength in V/A. Only read when `ExternalField = 1 EF`. For QW mode, the grid must not start at z=0 when the field is active. |

---

## 8. g-Factor Parameters

Optional. These are read after the external field block. If not present, the defaults below are used and the parser proceeds to SC parameters.

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `whichBand` | integer | `0` | 0, 1 | G | Band edge for g-factor calculation: 0 = conduction band, 1 = valence band. |
| `bandIdx` | integer | `1` | >= 1 | G | Subband doublet index: 1 = ground state, 2 = first excited, etc. |

---

## 9. Self-Consistent Schrodinger-Poisson Parameters

Optional. Triggered by `SC: 1`. All subsequent SC parameters use their defaults if missing or if the file ends early.

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `SC` | integer | `0` | 0, 1 | Q S | Enable self-consistent loop: 0 = off, 1 = on. All parameters below are only read when `SC = 1`. |
| `max_iter` | integer | `100` | >= 1 | S | Maximum number of SC iterations before giving up. |
| `tolerance` | float | `1.0e-6` | > 0 (eV) | S | Convergence threshold: maximum absolute change in electrostatic potential (infinity norm) in eV. |
| `mixing_alpha` | float | `0.3` | (0, 1] | S | Linear mixing parameter for the initial warm-up phase before DIIS kicks in. Smaller values are more stable but slower. |
| `diis_history` | integer | `7` | >= 1 | S | Number of previous iterations stored for DIIS/Pulay acceleration. Larger values can accelerate convergence but increase memory and may become unstable. |
| `temperature` | float | `300.0` | > 0 (K) | S | Temperature for Fermi-Dirac statistics in the charge density calculation. |
| `fermi_mode` | integer | `0` | 0, 1 | S | Fermi level determination: 0 = charge neutrality (automatically finds the Fermi level), 1 = fixed Fermi level (uses `fermi_level` value). |
| `fermi_level` | float | `0.0` | any (eV) | S | Fixed Fermi level in eV. Only used when `fermi_mode = 1`. |
| `num_kpar` | integer | `201` | >= 1, odd preferred | S | Number of in-plane wave-vector (k_parallel) sampling points for the charge density integration. Odd numbers enable Simpson's rule. |
| `kpar_max` | float | `0.0` | >= 0 (1/A) | S | Maximum in-plane wave-vector for charge density integration. 0 = automatic determination. |
| `bc_type` | string(2) | `DD` | `DD`, `DN` | S | Poisson boundary condition type: `DD` = Dirichlet-Dirichlet (potential fixed at both ends), `DN` = Dirichlet-Neumann (potential fixed at left, field=0 at right). |
| `bc_left` | float | `0.0` | any (eV) | S | Left boundary electrostatic potential in eV (Dirichlet). |
| `bc_right` | float | `0.0` | any (eV) | S | Right boundary electrostatic potential in eV (Dirichlet for `DD`, ignored for `DN`). |

---

## 10. Doping

Per-layer doping specifications, read only when `SC = 1`. One line per layer, 1-indexed.

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `dopingN` | 2 floats | `0.0 0.0` | >= 0 (cm^-3) | S | Doping for layer N: `dopingN: <ND> <NA>`. `ND` = donor concentration (n-type), `NA` = acceptor concentration (p-type). Read `numLayers` times. Defaults to `(0, 0)` if the line is missing. |

Example:
```
doping1: 0.0 0.0
doping2: 1.0e18 0.0
doping3: 0.0 0.0
```

---

## 11. Strain Parameters

Optional. Triggered by reading a strain-enabled line after the SC block. Uses defaults if missing.

| Name | Type | Default | Valid range | Modes | Description |
|---|---|---|---|---|---|
| `strain` | logical | `.false.` | `T`/`.true.`, `F`/`.false.` | Q W | Enable strain calculation. When enabled, the remaining strain parameters are read. |
| `strain_reference` | string(20) | `substrate` | `substrate`, material name | Q W | Reference lattice for strain calculation. `substrate` uses the first layer as reference. |
| `strain_solver` | string(20) | `pardiso` | `pardiso` | Q W | Strain solver backend. Currently only PARDISO direct solver is supported. |
| `strain_piezo` | logical | `.false.` | `T`/`.true.`, `F`/`.false.` | Q W | Include piezoelectric polarization in the strain calculation. |

---

## File Order Summary

Parameters must appear in the following sequence in `input.cfg`:

```
waveVector: <k0|kx|ky|kz>
waveVectorMax: <float>
waveVectorStep: <int>
confinement: <0|1|2>
FDstep: <int>
FDorder: <2|4|6|8|10>
numLayers: <int>
```

Then, **mode-dependent** block:

- **Bulk** (`confinement=0`): one `material1: <name>` line.
- **QW** (`confinement=1`): `numLayers` lines of `materialN: <name> <startPos> <endPos>`.
- **Wire** (`confinement=2`): grid parameters, shape, dimensions, `numRegions`, then `numRegions` region lines.

Then **common** block (always present):

```
numcb: <int>
numvb: <int>
ExternalField: <0|1> EF
EFParams: <float>
```

Then **optional** blocks, in order:

```
whichBand: <0|1>          ! g-factor (optional)
bandIdx: <int>            ! g-factor (optional)
SC: <0|1>                 ! self-consistent (optional)
  max_iter: <int>           }
  tolerance: <float>        }
  mixing_alpha: <float>     } SC parameters (only if SC=1)
  diis_history: <int>       }
  temperature: <float>      }
  fermi_mode: <0|1>         }
  fermi_level: <float>      }
  num_kpar: <int>           }
  kpar_max: <float>         }
  bc_type: <DD|DN>          }
  bc_left: <float>          }
  bc_right: <float>         }
  doping1: <ND> <NA>        }
  doping2: <ND> <NA>        } per-layer doping (only if SC=1)
  ...                       }
strain: <T|F>              ! strain (optional)
  strain_reference: <str>   }
  strain_solver: <str>      } strain parameters (only if strain=T)
  strain_piezo: <T|F>       }
```

---

## Notes

- The parser is **sequential**: lines must appear in the exact order shown above.
- Optional blocks (g-factor, SC, strain) can be omitted entirely; defaults from `defs.f90` are used.
- Comments use `!` prefix but are only reliable in positions the parser does not try to read.
- Material names are case-sensitive and must match the database in `src/core/parameters.f90` (e.g. `GaAs`, `InAsW`, `Al20Ga80As`).
- The `FDstep` field is named `fdStep` internally. For bulk it is forced to 1 regardless of input.
- For wire mode, `numLayers` is overridden to equal `numRegions` after parsing.
- Alloy materials like `Al20Ga80As` use linear interpolation between endpoint binaries (e.g. AlAs and GaAs with 20/80 composition).
