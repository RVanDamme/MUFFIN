<!-- Adapted for the Nextflow plugin: generic host tools and OAuth MCP. -->
<!-- Purpose: Migration-focused notes from official Snakemake testing docs. -->

# snakemake-testing-unit-tests

Official pattern for building small, reproducible tests from real workflow runs.

## Why this source helps

Snakemake can auto-generate Pytest rule tests from a successful run:

```bash
snakemake --generate-unit-tests
```

It generates per-rule tests and expected outputs captured from representative jobs; inspect the generated files in the location selected by the tool.

## Migration use

- Use generated tests to establish source behavior before conversion.
- Keep fixtures tiny and deterministic (`--notemp` + small data) so tests are fast and versionable.
- Port intent to `nf-test` as:
  - spec tests (expected behavior)
  - regression tests (bug reproductions)
- Reuse output comparison patterns when exact byte parity is required.

## Practical flow

1. Run Snakemake on small fixture data.
2. Generate Snakemake unit tests.
3. Run pytest to validate baseline.
4. Port equivalent checks to `nf-test` during DSL2 migration.

## Sources

- ReadTheDocs: https://snakemake.readthedocs.io/en/stable/snakefiles/testing.html
- Upstream source (authoritative): https://github.com/snakemake/snakemake/blob/main/docs/snakefiles/testing.rst
