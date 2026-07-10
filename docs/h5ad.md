# H5AD Adapter

Read and write AnnData H5AD files.

## Reading

```python
from scinterop.h5ad import read_h5ad

obj = read_h5ad("data/pbmc3k.h5ad")
```

All slots are mapped: `X`, `obs`, `var`, `obsm`, `layers`, `uns`, `raw`.

## Writing

```python
from scinterop.h5ad import write_h5ad

write_h5ad(obj, "output.h5ad")
```

The ``.h5ad`` extension is automatically appended if missing.
