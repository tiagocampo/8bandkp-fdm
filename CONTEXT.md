# 8-band k.p FDM solver

Fortran solver for the 8-band zinc-blende k·p Hamiltonian across four confinement
modes (bulk, quantum well, wire, Landau). This glossary fixes terms that are
overloaded across those modes — it is a glossary only, not a spec.

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
