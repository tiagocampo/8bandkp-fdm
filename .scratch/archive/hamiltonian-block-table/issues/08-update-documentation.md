## Parent

`.scratch/hamiltonian-block-table/PRD.md`

## What to build

Update all documentation that references the Hamiltonian module structure, strain/Zeeman tables, or module dependency graph to reflect the C4 refactoring. This is the final documentation sweep after all code changes are validated.

Files to update (in priority order):

**CLAUDE.md:**
- Directory listing (`src/physics/`): add `hamiltonian_blocks.f90`
- Module dependency graph: add `hamiltonian_blocks.f90` between `finitedifferences` and `hamiltonianConstructor`/`hamiltonian_wire`
- Key design concepts: update Hamiltonian description to mention table-driven block structure
- Boundaries section: add `get_strain_table()` and `get_zeeman_table()` as single-source-of-truth exports
- Known issues: update wire g-factor velocity operator reference if module changed

**README.md:**
- Architecture diagram: add `hamiltonian_blocks` to the dependency chain

**docs/lecture/12-extending-the-code.md:**
- Module dependency graph: add `hamiltonian_blocks.f90` and update edges
- Code examples referencing `hamiltonianConstructor.f90` for new couplings: update guidance (k.p terms go to table)

**docs/lecture/04-strain.md:**
- API table for `strain_solver`: add `strain_entry`, `get_strain_table`, `zeeman_entry`, `get_zeeman_table`

**docs/lecture/01-bulk-band-structure.md:**
- `ZB8bandBulk` references: update if location or interface changes

**docs/lecture/02-quantum-well.md:**
- `ZB8bandQW` and `confinementInitialization` references: update if interface changes

**docs/lecture/06-optical-properties.md:**
- Velocity matrix and Bir-Pikus references in `hamiltonianConstructor.f90`: update locations

**docs/lecture/08-quantum-wire.md:**
- `confinementInitialization` kpterms_2d reference: update if affected

**docs/lecture/10-qcse.md:**
- `externalFieldSetup_electricField` reference: update if location changes

**docs/lecture/13-topological-superconductivity.md:**
- `add_peierls_coo` reference: update if location changes

**docs/adr/0001-simulation-setup-fat-type.md:**
- Update C4 candidate status from "future" to "done"
- Update dispatch seam description if `setup_build_H` changed

**docs/solutions/logic-errors/bir-pikus-qeps-sign-bug-2026-05-18.md:**
- Update six-file coordination list (add `hamiltonian_blocks.f90`)
- Update stale comment references and line numbers
- Update single-source-of-truth to include `get_strain_table()`

**docs/solutions/best-practices/8band-verification-ladder-2026-05-09.md:**
- Update module name if file references changed

**docs/plans/archive/2026-04-26-hamiltonian-performance-refactor-plan.md:**
- Update line count expectations for both Hamiltonian modules

## Acceptance criteria

- [ ] CLAUDE.md dependency graph includes `hamiltonian_blocks.f90`
- [ ] CLAUDE.md directory listing includes `hamiltonian_blocks.f90`
- [ ] CLAUDE.md boundaries section lists `get_strain_table()` and `get_zeeman_table()` as single-source-of-truth
- [ ] README.md architecture diagram updated
- [ ] docs/lecture/12-extending-the-code.md module graph updated
- [ ] docs/lecture/04-strain.md API table includes new strain_solver exports
- [ ] docs/adr/0001 C4 status updated from "future" to "done"
- [ ] docs/solutions/logic-errors/bir-pikus-qeps-sign-bug-2026-05-18.md coordination list updated
- [ ] All other lecture docs referencing affected modules are updated
- [ ] No stale references to removed subroutines (`add_bp_strain_dense`, `apply_bp_strain_inline`)

## Blocked by

- `.scratch/hamiltonian-block-table/issues/07-validate-full-refactoring.md`
