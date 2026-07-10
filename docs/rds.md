# RDS Adapter

Read and write RDS files via an R subprocess bridge.

## Requirements

- R >= 4.3 with the ``anndata`` R package
- ``Seurat`` R package (for Seurat objects)

## Reading

```python
from scinterop.rds import read_rds

obj = read_rds("data/seurat_obj.rds")
```

The adapter writes an intermediate H5AD using R's `anndata` package,
then reads it into a CanonicalObject.

## Writing

```python
from scinterop.rds import write_rds

# Plain R list
write_rds(obj, "output.rds")

# Seurat object
write_rds(obj, "output.rds", seurat=True)
```
