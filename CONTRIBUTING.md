# Contributing

## Development setup

1. Clone the repository and start R from the repository root.
2. Install hard dependencies from `DESCRIPTION`.
3. Install only the optional `Suggests` needed for the module being changed.
4. Keep raw data and generated outputs outside Git.

## Pull requests

- Keep changes scoped and document scientific assumptions.
- Add or update a regression test for every behavior change.
- Preserve `pixel_id`, subject, section, coordinate-unit and provenance fields.
- Do not describe pixels, tiles or sections as independent biological subjects.
- Do not describe approximate-mass matches as confirmed structures.
- Label method-inspired approximations explicitly; do not imply exact upstream
  software reproduction.
- Run `Rscript scripts/run_tests.R` before opening a pull request.

## Data and privacy

Do not commit raw MSI files, downloaded annotation databases, private paths,
credentials, patient identifiers or controlled-access data. Use environment
variables and document public accessions and licenses instead.

## Style

- Use base R-compatible code unless a dependency is justified in `DESCRIPTION`.
- Validate dimensions, identifiers, units and exchangeability blocks at public
  API boundaries.
- Return explicit QC/status fields instead of silently replacing failed fits.
