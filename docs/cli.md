# CLI Reference

## `scinterop detect`

Detect the format of a single-cell data file or directory.

```bash
scinterop detect path/to/file
```

## `scinterop convert`

Convert between single-cell formats.

```bash
scinterop convert input.h5ad output.rds
scinterop convert input.h5ad output.rds --seurat
scinterop convert input.rds output.h5ad --r-exe /path/to/Rscript
```

### Options

- `--r-exe` — Path to Rscript executable.
- `--python-exe` — Path to Python executable.
- `--seurat` — Create a Seurat object when writing RDS.
- `--debug` — Keep temporary files for debugging.

## `scinterop run`

Execute an external R or Python script.

```bash
scinterop run script.R --executor r
scinterop run script.py --executor python
scinterop run script.py --executor python --conda-env myenv
```
