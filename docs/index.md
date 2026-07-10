# scinterop

**Single-cell data interoperability between R and Python formats.**

Convert seamlessly between H5AD (AnnData), RDS (Seurat), and 10X MTX
while preserving expression data, cell/gene metadata, dimensional
reductions, and expression layers.

## Quick start

```bash
pip install scinterop

# CLI usage
scinterop detect data/pbmc3k.h5ad
scinterop convert data/pbmc3k.h5ad output.rds
```

## Python API

```python
from scinterop import read, convert

obj = read("data/pbmc3k.h5ad")
print(f"Dataset: {obj.n_cells} cells x {obj.n_genes} genes")

convert("data/pbmc3k.h5ad", "output/seurat_obj.rds", seurat=True)
```

## Supported formats

| Format | Extension | Reader | Writer | Engine |
|--------|-----------|--------|--------|--------|
| AnnData / H5AD | `.h5ad` | ✅ | ✅ | Python (anndata) |
| Seurat / RDS | `.rds` | ✅ | ✅ | R (Seurat + anndata) |
| 10X MTX | `.mtx` / directory | ✅ | ✅ | Python (scipy.io) |
