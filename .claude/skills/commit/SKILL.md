# Commit Workflow

## Step 0 — Detect Branch

Run `git branch --show-current` to determine the current branch, then follow the corresponding workflow.

## Branch: `development` (Full Workflow)

1. Run `pytest` and ensure all tests pass (skip only if the user explicitly says so)
2. Update any `context4agent/` documentation that references changed modules
3. Stage all changes with `git add -A`
4. Write a conventional commit message summarizing changes
5. Push to `development`

## Branch: `main` (Public-Only Workflow)

On `main`, only user-facing public content should be committed. Developer-only private content must NOT be staged or committed.

### Public (user-facing) — DO commit
- `src/`
- `test/`
- `README.md`
- `pyproject.toml`
- `data_example/`
- `.gitignore`
- `.gitattributes`

### Private (developer-only) — DO NOT commit
- `CLAUDE.md`
- `.claude/`
- `context4agent/`
- `.vscode/`

### Steps
1. Identify which changed files are public vs private (use `git status`)
2. If there are no public changes, inform the user and stop
3. Stage ONLY the public files by name (never use `git add -A` on main)
4. Write a conventional commit message summarizing the public changes
5. Push to `main`
6. If private files were skipped, list them and inform the user
