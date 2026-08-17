# Issue tracker: GitHub

Issues and specs for this repo live in the `adamrtalbot/hap.py` GitHub repository. Use the `gh` CLI for all operations and pass `--repo adamrtalbot/hap.py` explicitly because this clone also has an upstream remote.

## Conventions

- **Create an issue**: `gh issue create --repo adamrtalbot/hap.py --title "..." --body "..."`. Use a heredoc for multi-line bodies.
- **Read an issue**: `gh issue view <number> --repo adamrtalbot/hap.py --comments`, filtering comments by `jq` and also fetching labels.
- **List issues**: `gh issue list --repo adamrtalbot/hap.py --state open --json number,title,body,labels,comments --jq '[.[] | {number, title, body, labels: [.labels[].name], comments: [.comments[].body]}]'` with appropriate `--label` and `--state` filters.
- **Comment on an issue**: `gh issue comment <number> --repo adamrtalbot/hap.py --body "..."`
- **Apply or remove labels**: `gh issue edit <number> --repo adamrtalbot/hap.py --add-label "..."` / `--remove-label "..."`
- **Close**: `gh issue close <number> --repo adamrtalbot/hap.py --comment "..."`

## Pull requests as a triage surface

**PRs as a request surface: no.** _(Set to `yes` if this repo treats external PRs as feature requests; `/triage` reads this flag.)_

When set to `yes`, PRs run through the same labels and states as issues, using the `gh pr` equivalents:

- **Read a PR**: `gh pr view <number> --repo adamrtalbot/hap.py --comments` and `gh pr diff <number> --repo adamrtalbot/hap.py`.
- **List external PRs for triage**: `gh pr list --repo adamrtalbot/hap.py --state open --json number,title,body,labels,author,authorAssociation,comments`, retaining only `CONTRIBUTOR`, `FIRST_TIME_CONTRIBUTOR`, or `NONE`.
- **Comment, label, or close**: use `gh pr comment`, `gh pr edit`, or `gh pr close` with `--repo adamrtalbot/hap.py`.

GitHub shares one number space across issues and PRs. Resolve an ambiguous `#42` with `gh pr view 42 --repo adamrtalbot/hap.py`, then fall back to `gh issue view 42 --repo adamrtalbot/hap.py`.

## When a skill says “publish to the issue tracker”

Create an issue in `adamrtalbot/hap.py`.

## When a skill says “fetch the relevant ticket”

Run `gh issue view <number> --repo adamrtalbot/hap.py --comments`.

## Wayfinding operations

Used by `/wayfinder`. The map is a single issue with child issues as tickets.

- **Map**: an issue labelled `wayfinder:map`, holding the Notes, Decisions-so-far, and Fog body.
- **Child ticket**: an issue linked to the map as a GitHub sub-issue. Where sub-issues are unavailable, add it to a task list in the map and put `Part of #<map>` at the top of the child body. Use a `wayfinder:<type>` label.
- **Blocking**: use GitHub’s native issue dependencies. Where unavailable, use a `Blocked by: #<n>` line at the top of the child body.
- **Frontier query**: list the map’s open children, remove assigned or blocked issues, and select the first in map order.
- **Claim**: `gh issue edit <n> --repo adamrtalbot/hap.py --add-assignee @me`.
- **Resolve**: comment with the answer, close the child, then append a context pointer to the map’s Decisions-so-far.
