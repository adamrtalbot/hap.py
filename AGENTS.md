## Agent skills

### Issue tracker

Issues are tracked in the `adamrtalbot/hap.py` GitHub repository. See `docs/agents/issue-tracker.md`.

### Triage labels

The repository uses the five default triage labels. See `docs/agents/triage-labels.md`.

### Domain docs

This is a single-context repository. See `docs/agents/domain.md`.

## Working rules

Never grill a question whose answer is measurable. Run the gate, record the
measurement in `verification/README.md`, then decide once. See
`docs/adr/0009-close-the-map-and-measure-what-is-left.md`.

An ADR states a decision and the rule it creates. Measurements belong in
`verification/README.md`, because an ADR that carries them becomes the input to
the next session.

The wayfinder map is closed. The programme is
[spec: the hap-rs release validation programme](https://github.com/adamrtalbot/hap.py/issues/37);
work runs from there, not from a decision map.
