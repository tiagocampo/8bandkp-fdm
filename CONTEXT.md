# 8-band k.p FDM solver

Fortran solver for the 8-band zinc-blende k·p Hamiltonian across four confinement
modes (bulk, quantum well, wire, Landau). This glossary fixes terms that are
overloaded across those modes — it is a glossary only, not a spec.

> **Companion file.** `UBIQUITOUS_LANGUAGE.md` is the **domain-language
> backbone** (bands, confinement structures, k·p model, observables) this
> ambiguity layer sits on top of. The two are a deliberate split: that file
> defines the terms; this one disambiguates them where a word is overloaded
> across modes.

## Language

**Requested band window**:
The subset of the band spectrum a calculation targets: the `num_cb` lowest
conduction bands plus the `num_vb` highest valence bands
(`evnum = num_cb + num_vb`).
Its meaning depends on confinement. For QW, wire, and Landau it is a **compute
directive** — it selects the `il:iu` eigenvalue window actually extracted. For
bulk it is a **display filter only**, because bulk is always fully diagonalized
(all 8 bands are computed regardless; `evnum` only limits how many are written).
_Avoid_: "number of eigenvalues", "nev" — `nev` is the solver-internal count
*returned*, which may exceed or fall short of the requested window and is
clamped at the call site.

**AUTO** (method or mode = `"AUTO"`):
"Resolve to a valid default for the resolved method/confinement pair" — not
"the confinement's default, independently." Method resolves first (confinement
default unless overridden); mode then resolves to a default *compatible with
that method* (FEAST → ENERGY; DENSE → the confinement's native mode). The
contract is that `AUTO` never produces an invalid method×mode combination.
_Avoid_: "default" unqualified — a method default and a mode default are two
different resolutions and must not be picked independently.

**Eigensolver dispatch** (format vs backend):
At a solve call site, the caller selects the entry point by the **format of the
matrix it holds** — a dense array or a CSR (sparse) matrix — never by which
backend will run. The **backend** (dense LAPACK vs FEAST) is fixed when the
solver is constructed (by `[solver].method`, resolved through `AUTO`) and is
invisible at the call site; each backend accepts both formats, converting
internally. Format and backend are **independent axes**: dense/CSR is what the
caller holds, DENSE-LAPACK/FEAST is who does the work.
_Avoid_: "the solve path" or "the eigensolver path" unqualified, and "solve is
the legacy alias" — naming one axis without saying which has caused repeated
confusion (the most-used entry point got mislabeled *legacy* because the two
axes were conflated).

**Block** (always qualified — the bare word is overloaded): three senses that
must not be conflated. A **block-table entry** is one row of the k·p / strain /
Zeeman topology tables — a static descriptor of how a single band-pair couples
(the single-source-of-truth *data*). A **spatial block** is an N×N region of a
built QW/wire Hamiltonian (the 8N×8N matrix seen as an 8×8 grid of N×N
sub-matrices), with its scalar analog one (band-i, band-j) entry of the bulk
8×8. A **Nambu block** is one of the four 8N×8N quarters of the BdG 16N×16N
matrix, unrelated to the topology tables. The k·p-term *formula* (Q−T,
½(Q+T)) is yet another descriptor — the k·p-term descriptor — not a block.
_Avoid_: bare "block" for a matrix region when a table entry could be meant;
"block-formula" for the recipe (say k·p-term descriptor).

**Geometry**: reserved for the **wire cross-section shape** (rectangle / circle /
hexagon / polygon in the x–y plane) — the `geometry` module, the
`[wire.geometry]` section, and the `wire_geometry` type are all wire-specific
despite the general-sounding module name. For the spatial layout implied by
*any* confinement (bulk's single point, QW's z-grid, wire's x–y grid, Landau's
x-grid), the word is **confinement** / **confinement structure**, not "geometry".
_Avoid_: "geometry" for the QW z-grid or any non-wire spatial layout.

**Confinement direction** (`conf_direction`): the single FDM-discretized
confinement axis — meaningful only for the 1D-confinement modes: QW → `'z'`
(growth axis), Landau → `'x'`, bulk → `'n'` (none). The wire is 2D-confined in
the x–y plane, so it has *no* single confinement axis; `conf_direction('wire')`
returns the sentinel `'w'`, which no live path consumes (every wire path
branches to wire-specific code). Note the wire's free/propagation axis is z —
the opposite role from QW's `'z'`.
_Avoid_: treating `'w'` as an axis or as `'n'` (the wire IS confined); a generic
`conf_direction == 'z'` branch without a prior wire guard (it would misread the
wire as a QW and flip the axis). "Growth direction" (strain's term, = z) and
"propagation direction" (the wire's free z-axis) name different roles than the
QW confinement direction, even when they coincide on z.

**Mode**: overloaded across three namespaces — disambiguate by context. (1)
**Confinement mode** — an informal synonym for a `confinement` value (bulk / qw
/ wire / landau); it is *not* a code field (the field is `confinement`). (2)
**Eigensolver mode** — `solver%mode`, the eigenvalue selection (`FULL` / `INDEX`
/ `ENERGY` / `AUTO`). (3) **Topology mode** — `topo%mode`, the invariant class
(`qhe` / `qshe` / `bdg` / ...). "QW mode" / "wire mode" mean confinement mode;
"FEAST mode" / "energy mode" mean eigensolver mode.
_Avoid_: bare "mode" in comments near solver or topology code without a
qualifier — three different fields answer to the same word. "confinement_mode"
as if it were a field (it is not; the field is `confinement`).

**Velocity operator**: the matrix v_α = −i[r_α, H]/ħ, built two *equivalent*
ways — the **commutator** −i[r_α, H] element-wise on CSR (for FDM-discretized
directions, `build_velocity_matrices`) and the **k-derivative** ∂H/∂k_α via the
`g='g'` / `g='g3'` modes (for translationally free directions). These are the
same operator — the k·p identity ∂H/∂k = −i[r, H]; the source is chosen by which
direction is discretized. The raw built object carries one factor of ħ: it is
M_α = ħ v_α = (ħ/m₀) p_α (units eV·Å), i.e. *proportional to* velocity, not
equal to it. **Momentum matrix** p_α = m₀ v_α is the same object rescaled. Each consumer
cancels the raw matrix's ħ² at a different stage — gfactor inside the
oscillator-strength formula (f = (2/m₀)|p|²/ΔE via `hbar2O2m0`), optics inside
the absorption prefactor (`optics_apply_prefactor`: 2πe²/(n_r·c·ε₀·ħ²·E)). Both
yield correct absolute results.
_Avoid_: treating the commutator and k-derivative sources as different operators
(they are identical); "g1" / "g2" modes (absent from the source — only `g='g'`
and `g='g3'` exist); "current operator" (not used in this codebase).
