# Triage Labels

The skills speak in terms of five canonical triage roles. This file maps those roles to the actual label strings used in this repo's issue tracker.

| Label in mattpocock/skills | Label in our tracker | Meaning                                  |
| -------------------------- | -------------------- | ---------------------------------------- |
| `needs-triage`             | `needs-triage`       | Maintainer needs to evaluate this issue  |
| `needs-info`               | `needs-info`         | Waiting on reporter for more information |
| `ready-for-agent`          | `ready-for-agent`    | Fully specified, ready for an AFK agent  |
| `ready-for-human`          | `ready-for-human`    | Requires human implementation            |
| `wontfix`                  | `wontfix`            | Will not be actioned                     |
| (closed work, not in skill vocab) | `closed`       | Work shipped; spec/issue kept as documentation; not awaiting triage |

When a skill mentions a role (e.g. "apply the AFK-ready triage label"), use the corresponding label string from this table. Note: `closed` is a project-local extension for specs / tickets whose work has shipped; the upstream skills never reference it, so applying it requires explicit context (the as-built spec preamble flags it).

Edit the right-hand column to match whatever vocabulary you actually use.
