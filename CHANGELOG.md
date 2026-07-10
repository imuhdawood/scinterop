# Changelog

## v0.2.0 (2026-07-10)

### Added
- Support for `.qs` (qs serialization) files via the RDS adapter --
  `si.read("file.qs")` and `si.rds.write_rds(obj, "out.qs")` now work
  with the `qs` R package
- `r-qs` added to `environment.yml` for automatic installation
- `CITATION.cff` for GitHub citation metadata
- `LICENSE` (MIT)
- CI/CD workflows: tests (push/PR), docs (GitHub Pages), publish (PyPI)
- MkDocs documentation site at https://imuhdawood.github.io/scinterop/
- Google-style docstrings on all 13 modules

### Changed
- README restructured: conda env creation from GitHub URL, merged
  full-env section into quick start
- README polished: removed conversational phrasing and em-dashes,
  added Citation section
- R bridge scripts detect `.qs` extension and use `qread()`/`qsave()`
  transparently

## v0.1.0 (2026-07-10)

Initial release.

- CanonicalObject data model with validation
- Format adapters: H5AD, 10X MTX, RDS (Seurat via R bridge)
- Auto-detection of file formats
- Subprocess runners for R and Python
- Scratch manager for temporary files
- Provenance logging to JSONL
- CLI tool (`scinterop detect`, `convert`, `run`)
