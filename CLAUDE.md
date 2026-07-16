# Project Instructions for AI Agents

This file provides instructions and context for AI coding agents working on this project.

This project uses **bd** (beads) for issue tracking. Run `bd onboard` to get started.

This checkout uses the Beads 1.x embedded-Dolt format. Before mutating issues,
run `bd version`; do not use a legacy 0.x binary against this repository. If
`type -a bd` reports more than one installation, explicitly select the 1.x
binary before running any Beads command.

On a fresh clone, restore the issue database from the pushed Dolt history before
claiming work:

```bash
chmod 700 .beads
bd bootstrap --yes
bd prime
```

If bootstrap cannot authenticate, follow the forwarded-SSH-agent procedure
below and retry. The tracked `.beads/issues.jsonl` is a portable checkpoint and
interchange copy; `bd bootstrap` is the preferred history-preserving restore.

## Branch and durability policy

- Do implementation work on a dedicated topic branch, never directly on `main` or
  `master`.
- Make reviewable commits and push the topic branch periodically after coherent,
  tested milestones, as well as at session completion, so work is not stranded on
  one machine.
- If an SSH push fails, inspect the currently forwarded agent and its identities
  (`SSH_AUTH_SOCK` and `ssh-add -l`), obtain or select the latest forwarded key,
  and retry the push.  Never place a private key in the repository or logs.
- A current explicit user or orchestrator instruction not to push still takes
  precedence over this standing policy.

## Landing the Plane (Session Completion)

When ending a work session, complete all of these steps. Work is not complete
until the required push succeeds.

1. File Beads issues for remaining work.
2. Run relevant builds, tests, linters, and other quality gates.
3. Close finished issues and accurately update unfinished claims.
4. Commit coherent changes, push Beads Dolt history, and push the topic branch:

   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status
   ```

5. Verify that the tree is clean and the branch is up to date with its remote.
6. Leave a durable handoff with exact validation and remaining blockers.

Never stop at "ready to push" when current instructions authorize publishing.
If a push fails, resolve authentication or transport and retry; do not leave
work stranded locally.

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:6cd5cc61 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->


## Build & Test

_Add your build and test commands here_

```bash
# Example:
# npm install
# npm test
```

## Architecture Overview

_Add a brief overview of your project architecture_

## Conventions & Patterns

_Add your project-specific conventions here_
