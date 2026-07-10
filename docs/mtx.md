# MTX Adapter

Read and write 10X Genomics Market Exchange format.

## Reading

```python
from scinterop.mtx import read_mtx

# From a 10X directory
obj = read_mtx("data/filtered_feature_bc_matrix/")

# From a single .mtx file
obj = read_mtx("data/matrix.mtx.gz")
```

## Writing

```python
from scinterop.mtx import write_mtx

write_mtx(obj, "output/mtx_dir/")
```

Creates a directory with `matrix.mtx`, `barcodes.tsv`, and `features.tsv`.
