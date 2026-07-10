# Format Detection

The `detect()` function identifies data formats without reading the
full file. It uses file extensions and directory contents.

## Resolution order

1. **Directory** → check for 10X MTX trio files.
2. **Existing file** → check extension.
3. **Non-existent path with extension** → infer from extension.
4. Otherwise → raise `DetectionError`.

## Supported extensions

| Extension | Format |
|-----------|--------|
| `.h5ad` / `.h5` | H5AD (AnnData) |
| `.rds` | RDS (Seurat / R list) |
| `.mtx` / `.mtx.gz` | 10X MTX |
