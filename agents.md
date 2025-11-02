# Agents

## Overview
- Tracks the people and automation that help operate the `8bandkp-fdm` project.
- “Agent” here means any entity—human or software—that can run workflows, modify code, or provide support.
- Use this document to know who to contact and how to collaborate or extend the automation surface.

## Active Agents

### Human Maintainer
- **Tiago de Campos**  
  - Role: Project author and primary reviewer.  
  - Scope: Reviews Fortran changes, curates material parameter sets, maintains scientific references.  
  - Contact: Follow the channels defined in the parent project (email or issue tracker).

### LLM Support
- **Codex CLI (GPT-based)**  
  - Role: Drafts documentation, assists with Fortran refactors, suggests Makefile tweaks, and prepares input examples.  
  - Invocation: Use the Codex CLI session from the project root; provide clear prompts and mention expected outputs (code, docs, analysis).  
  - Guardrails: Keeps modifications within the repository workspace; relies on maintainers for domain validation and publication decisions.

### Manual Operations
- **Build Operators**  
  - Role: Developers or researchers who compile the solvers locally or on HPC resources.  
  - Scope: Run `make all`, execute `./bandStructure` and `./gfactorCalculation`, archive raw outputs.  
  - Notes: Ensure required libraries (FFTW, BLAS/LAPACK or MKL) are available on the target machine.

## Interaction Playbooks
- Code updates: open a branch from the target commit, run `make all`, validate results with representative inputs (`input.cfg`, `quantumwell.example`).  
- Documentation refresh: coordinate with Tiago de Campos for scientific accuracy; Codex CLI can draft sections before review.  
- Data runs: document input files and solver versions alongside generated datasets.

## Onboarding New Agents
1. Define the agent type (human, automation, LLM) and intended responsibilities.  
2. Record access requirements (repository rights, compute resources).  
3. Update this file with contact info, scope, and guardrails.  
4. Walk through compilation and execution steps to ensure the agent can reproduce baseline results.

## Maintenance
- Review this document whenever contributors change or automation is added/retired.  
- At project milestones (e.g., publication, release), confirm that agent responsibilities and playbooks still match current workflows.

