# Test Directory Layout

This directory contains test scaffolding for the current code under
`scripts/structure/utils/`.

Recommended structure:

```text
test/
  README.md
  conftest.py
  unit/
    structure/
      utils/
        test_layer_parser.py
        test_water_parser.py
  integration/
    structure/
      utils/
        test_water_layer_pipeline.py
  fixtures/
    atoms/
      README.md
    expected/
      README.md
```

Notes:

- `unit/`: small, deterministic tests for single functions.
- `integration/`: frame-level workflows that combine layer + water analysis.
- `fixtures/atoms/`: minimal ASE frame fixtures (tiny synthetic systems).
- `fixtures/expected/`: expected arrays/profiles used for regression checks.
- Keep z-profile tests explicit in naming (for example `*_z_distribution_*`)
  so they remain distinct from future "z-window statistic" tests.
