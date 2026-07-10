# scinterop: Single-Cell Interoperability

Bidirectional conversion between single-cell data formats (H5AD, 10X MTX, RDS/Seurat) with explicit environment paths, provenance logging, and scratch management.

**Core design principles:**
- No implicit global state: all paths are explicit
- Separate format adapters from execution runners
- No `rpy2` dependency: R bridge uses subprocess + intermediate H5AD
- All failures localize to the adapter that failed
- Every conversion logs a JSON provenance record

---

## Table of Contents

- [Install](#install)
- [Architecture overview](#architecture-overview)
- [CanonicalObject, the core data structure](#canonicalobject--the-core-data-structure)
- [Format detection](#format-detection)
- [Format adapters](#format-adapters)
  - [H5AD adapter (`scinterop.h5ad`)](#h5ad-adapter-scinteroph5ad)
  - [10X MTX adapter (`scinterop.mtx`)](#10x-mtx-adapter-scinteromtx)
  - [RDS adapter (`scinterop.rds`)](#rds-adapter-scinterords)
- [External script runners](#external-script-runners)
  - [R runner (`scinterop.r_runner`)](#r-runner-scinteropr_runner)
  - [Python runner (`scinterop.python_runner`)](#python-runner-scinteroppython_runner)
- [Scratch manager (`scinterop.cache`)](#scratch-manager-scinteropcache)
- [Provenance logging](#provenance-logging)
- [Validation](#validation)
- [Error hierarchy](#error-hierarchy)
- [CLI reference](#cli-reference)
- [Environment variables](#environment-variables)
- [Smoke test](#smoke-test)
- [Extending with a new format](#extending-with-a-new-format)

---

## Install

### Quick start (conda)

Create the environment directly from the GitHub-hosted `environment.yml`, then install scinterop:

```bash
conda env create -n scinterop -f https://raw.githubusercontent.com/imuhdawood/scinterop/main/environment.yml
conda activate scinterop
```

This gives you Python + R + Seurat + anndata (R package) in one shot -- everything needed for all format adapters.

### Minimal (pip only)

If you already have a Python environment with numpy, scipy, and pandas:

```bash
# Install from GitHub
pip install git+https://github.com/imuhdawood/scinterop.git

# Editable install (for development)
pip install -e /path/to/scinterop

# Optional: add H5AD support
pip install anndata
```

**Runtime deps:** `numpy`, `scipy`, `pandas`
**Optional:** `anndata` (for H5AD); `R` + `Seurat` + `anndata` R package (for RDS)

---

## Architecture overview

```
┌────────────────────────────────────────────────────────────────────┐
│                      scinterop.read() / convert()                  │
│                        (auto-detect + dispatch)                    │
├──────────┬──────────┬──────────┬──────────┬──────────┬────────────┤
│  detect  │ validate │   h5ad   │   mtx    │   rds    │ provenance │
│  .py     │  .py     │   .py    │   .py    │   .py    │   .py      │
├──────────┴──────────┴──────────┴──────────┴──────────┴────────────┤
│  r_runner.py  │  python_runner.py  │  cache.py  │  schema.py     │
│  (subprocess) │  (subprocess/conda) │  (scratch) │  (dataclass)   │
└───────────────┴────────────────────┴────────────┴────────────────┘
```

**Three API layers (all consistent):**

| Layer | Example | Purpose |
|-------|---------|---------|
| Object methods | `obj.to_seurat("out.rds")` | Convenience on CanonicalObject |
| Top-level auto-detect | `si.convert("in.h5ad", "out.rds")` | One-shot format conversion |
| Per-format explicit | `si.rds.write(obj, "out.rds")` | Fine control when you know the format |

---

## CanonicalObject, the core data structure

A plain Python dataclass with no external dependencies beyond `numpy`, `scipy`, `pandas`. This is the neutral representation that all adapters convert to and from.

```python
from scinterop import CanonicalObject
import numpy as np
import pandas as pd
from scipy import sparse as sp

obj = CanonicalObject(
    X=sp.csr_matrix(np.random.poisson(0.5, size=(100, 2000))),
    obs=pd.DataFrame({"cell_type": ["T cell"] * 100}),
    var=pd.DataFrame(index=[f"gene_{i}" for i in range(2000)]),
    obsm={"X_pca": np.random.randn(100, 10)},
    layers={"counts": sp.csr_matrix(np.random.poisson(0.5, size=(100, 2000)))},
    uns={"genome": "hg38", "params": {"n_pcs": 10}},
    raw=CanonicalObject(X=np.random.poisson(0.5, size=(100, 2000))),
)
```

### Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `X` | `np.ndarray` or `sp.spmatrix` | Yes | Expression matrix (cells × genes) |
| `obs` | `pd.DataFrame` | No | Cell-level metadata, one row per cell |
| `var` | `pd.DataFrame` | No | Gene-level metadata, one row per gene |
| `obsm` | `dict[str, np.ndarray]` | No | Multi-dimensional cell embeddings (X_pca, X_umap, etc.) |
| `layers` | `dict[str, np.ndarray or sp.spmatrix]` | No | Additional expression matrices, same shape as X |
| `uns` | `dict` | No | Unstructured metadata (parameters, colors, etc.) |
| `raw` | `CanonicalObject or None` | No | Raw counts reference |

### Properties and methods

```python
obj.shape       # -> (n_cells, n_genes)
obj.n_cells     # -> int
obj.n_genes     # -> int
obj.copy()      # -> deep copy

# Export convenience methods (delegate to format adapters)
obj.to_anndata("export.h5ad")
obj.to_seurat("export.rds")       # Seurat object in RDS
obj.to_rds("export.rds")          # Plain R list in RDS
obj.to_rds("export.rds", seurat=True)  # Same as to_seurat()
obj.to_mtx("export_mtx/")         # 10X-style MTX directory
```

### Auto-validation

The constructor validates that `obs` and `var` lengths match `X.shape`. See [Validation](#validation) for explicit checks.

---

## Format detection

The `detect()` function identifies the format of a file or directory without reading the full data.

```python
from scinterop import detect

result = detect("data.h5ad")
result.fmt     # -> "h5ad"
result.path    # -> Path("data.h5ad")
result.details # -> {"extension": ".h5ad"}
```

### Detection logic

1. **If path is a directory**: check for 10X MTX trio files (`matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, `features.tsv[.gz]`)
2. **If file exists**: check extension, optionally peek HDF5 contents for anndata vs Seurat signatures
3. **If file doesn't exist but has extension**: infer format from extension alone (useful for planning conversions)
4. **Anything else**: raise `DetectionError` with supported formats listed

### Supported formats

| Extension | Format constant | Also matches |
|-----------|----------------|--------------|
| `.h5ad` | `"h5ad"` | `.h5` (with anndata signature) |
| `.rds` | `"rds"` | (none) |
| `.mtx` | `"mtx"` | `.mtx.gz` |
| Directory | `"mtx"` | Contains 10X trio files |

### Detection errors

```python
try:
    result = detect("file.xyz")
except DetectionError as e:
    print(e)  # "Cannot determine format from path: file.xyz. Supported extensions: .h5ad, .rds, .mtx"
```

---

## Format adapters

Each format adapter is a standalone module with `read()` and `write()` functions.
They are callable directly or via the top-level `read()`/`convert()` API.

### H5AD adapter (`scinterop.h5ad`)

Read and write AnnData H5AD files. Requires `anndata` pip package.

```python
import scinterop as si

# Read H5AD
obj = si.h5ad.read_h5ad("data.h5ad")

# Write H5AD (auto-adds .h5ad extension if missing)
si.h5ad.write_h5ad(obj, "output.h5ad")
```

**What gets mapped:**
| AnnData slot | CanonicalObject field |
|--------------|----------------------|
| `X` | `X` |
| `obs` | `obs` |
| `var` | `var` |
| `obsm` | `obsm` |
| `layers` | `layers` |
| `uns` | `uns` |
| `raw` | `raw` (recursively) |

If `obs` or `var` are empty DataFrames, default indices (`cell_0`, `gene_0`, ...) are auto-generated for writing.

**Error example:**
```
H5adAdapterError: Failed to read H5AD file 'corrupt.h5ad': Unable to open file (File signature not found)
```

---

### 10X MTX adapter (`scinterop.mtx`)

Read and write 10X Genomics MTX format (directory with `matrix.mtx`, `barcodes.tsv`, `features.tsv`).

```python
# Read from 10X output directory
obj = si.mtx.read_mtx("cellranger_output/")

# Read single .mtx file (no barcodes/features; auto-generates names)
obj = si.mtx.read_mtx("matrix.mtx.gz")

# Write to 10X format
si.mtx.write_mtx(obj, "output_mtx/")
```

**Convention handling:**
- Files use the 10X convention **genes × cells** (features in rows, cells in columns)
- CanonicalObject stores **cells × genes**
- The adapter transposes automatically on both read and write

**Output directory structure:**
```
output_mtx/
├── matrix.mtx       # Sparse matrix (genes × cells)
├── barcodes.tsv     # Cell barcodes
└── features.tsv     # Gene identifiers (tab-delimited: id\tname\t"Gene Expression")
```

**Error example:**
```
MtxAdapterError: 10X MTX directory 'empty_dir/' is missing: matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz]
```

---

### RDS adapter (`scinterop.rds`)

Read and write RDS files via a subprocess R bridge. No `rpy2` dependency.

**How it works:**
1. The adapter writes a temporary R script as a template string
2. Runs it via `scinterop.r_runner.run_r()` using `Rscript`
3. The R script uses the `anndata` R package to convert between RDS and H5AD
4. The intermediate H5AD is read/written by the H5AD adapter
5. Temp files are cleaned up via the scratch manager

```
write_rds:
  CanonicalObject → [temp.h5ad] → R script reads H5AD, builds Seurat/list → output.rds

read_rds:
  input.rds → R script reads object, saves as H5AD → [temp.h5ad] → CanonicalObject
```

```python
# Write as Seurat object
si.rds.write_rds(obj, "seurat.rds", seurat=True)
obj.to_seurat("seurat.rds")           # same

# Write as plain R list (components: X, obs, var, obsm, layers, uns)
si.rds.write_rds(obj, "data.rds", seurat=False)
obj.to_rds("data.rds")                # same (default: seurat=False)

# Read RDS (auto-detects Seurat vs list)
obj = si.rds.read_rds("seurat_object.rds")
```

**What the R bridge does for Seurat output:**
1. Transposes expression matrix (cells×genes → genes×cells)
2. Converts to `dgCMatrix` (column-compressed sparse)
3. Creates `SeuratObject` with counts and metadata
4. Adds reductions from `obsm` (maps `X_pca` → `pca` DimReduc, etc.)

**What the R bridge does for R list output:**
1. Saves `X`, `obs`, `var`, `obsm`, `layers`, `uns` as a named R list
2. Access from R with `result$X`, `result$obs`, etc.

**Required R packages:**
- Reading: `anndata` (for both Seurat and list)
- Writing (Seurat): `anndata`, `Seurat`, `Matrix`
- Writing (list): `anndata`

**Installing R dependencies (conda, recommended):**
```bash
conda install -n scinterop r-base r-seurat r-matrix r-anndata \
  -c conda-forge -c bioconda -y
```

This installs R 4.5 + Seurat 5.5 + Matrix + anndata R package in your scinterop env, no separate CRAN install or compilation needed.

**Verification:**
```bash
conda run -n scinterop R -e 'library(anndata); library(Seurat); cat("OK\n")'
```

The R script fails with a clear message if packages are missing:
```
RdsAdapterError: Failed to write RDS file 'out.rds': R script failed with code 1.
  stderr:
    Error: R package 'Seurat' is required. Install with: install.packages('Seurat')
```

**Seurat v5 compatibility:**
The adapter auto-detects the SeuratObject version at runtime and uses the correct API: `layer="counts"` for Seurat v5+, `slot="counts"` for Seurat v4.

---

## External script runners

Run arbitrary R or Python scripts with explicit environment paths.
All output is captured to strings and optionally to log files.

### R runner (`scinterop.r_runner`)

```python
from scinterop import run_r

result = run_r(
    "library(Seurat); obj <- readRDS('data.rds'); print(dim(obj))",
    r_exe="/path/to/Rscript",
    log_path="/tmp/r_log.txt",
    timeout=300,
)
print(result.stdout)
print(result.stderr)
```

The `r_exe` parameter is resolved in this order:
1. Explicit argument: `r_exe="/opt/R/4.3/bin/Rscript"`
2. Environment variable: `SCINTEROP_R_EXE`
3. Default: `Rscript` (must be on PATH)

### Python runner (`scinterop.python_runner`)

```python
from scinterop import run_python

# Direct python
result = run_python("script.py", python_exe="/opt/python/3.11/bin/python")

# Conda environment
result = run_python(
    "import scanpy as sc; adata = sc.read_h5ad('data.h5ad'); print(adata)",
    conda_env="sc-env",
)

# With arguments
result = run_python("script.py", args=["--input", "data.h5ad", "--output", "out.rds"])
```

The `python_exe` / `conda_env` resolution:
1. If `conda_env` is set → use `conda run -n <env> python`
2. If `python_exe` is set → use it directly
3. Fallback: `SCINTEROP_PYTHON_EXE` env var → `python`

### Error handling for both runners

```python
try:
    run_r("nonexistent_command()")
except RExecError as e:
    print(e)  # "R script failed with code 1. stderr: ..."
```

Both runners raise `ExecError` on:
- Non-zero return code (with stderr in the error message)
- Timeout (configurable, default 600s)
- Executable not found (with hint to set env var)

---

## Scratch manager (`scinterop.cache`)

Manages temporary directories for intermediate files during conversion.

```python
from scinterop.cache import ScratchManager

mgr = ScratchManager()  # Uses SCINTEROP_SCRATCH or /tmp/scinterop_USER/

ctx = mgr.tempdir("my_conversion")
# ctx.path -> PosixPath('/tmp/scinterop_user/my_conversion_a1b2c3d4e5f6')

# Write a temp file
script_path = mgr.write_script(ctx, "print('hello')", "myscript")

# Clean up (removes directory)
mgr.cleanup(ctx)

# Or keep for debugging (controlled by SCINTEROP_DEBUG=1)
mgr.cleanup(ctx, keep=True)
```

**SCINTEROP_DEBUG behaviour:**
- `SCINTEROP_DEBUG=1` → all temp files are preserved after conversion
- Default (unset) → temp files are removed on success

---

## Provenance logging

Every conversion creates a JSON record with timestamps, versions, and system info.

```python
from scinterop.provenance import write_conversion_record, read_provenance_log

# Write a record (appends to provenance.jsonl in log_dir)
write_conversion_record(
    input_path="/data/input.h5ad",
    input_format="h5ad",
    output_path="/data/output.rds",
    output_format="rds",
    success=True,
    runtime_s=12.5,
    log_dir="/data/logs/",
)

# Read back all records
records = read_provenance_log("/data/logs/provenance.jsonl")
for r in records:
    print(r["timestamp"], r["input_format"], "->", r["output_format"])
```

**Record format:**
```json
{
  "timestamp": "2026-07-09T14:30:00+00:00",
  "package": "scinterop",
  "version": "0.1.0",
  "input_path": "/data/input.h5ad",
  "input_format": "h5ad",
  "output_path": "/data/output.rds",
  "output_format": "rds",
  "success": true,
  "runtime_seconds": 12.5,
  "system": {
    "platform": "Linux",
    "release": "5.14.0-1.el9.x86_64",
    "python": "3.12.12",
    "hostname": "compute-node-01"
  }
}
```

The `read()` and `convert()` functions at the top level automatically write provenance records.

---

## Validation

Structural validation of a `CanonicalObject` without loading format-specific libraries.

```python
from scinterop import validate, assert_valid, ValidationError

result = validate(obj)
result.valid    # True / False
result.errors   # list of error strings

# Raises ValidationError on first issue
assert_valid(obj)
```

**Checks performed:**

| Check | Condition | Error message |
|-------|-----------|---------------|
| Type | Is it a `CanonicalObject`? | `Expected CanonicalObject, got <type>` |
| X type | `np.ndarray` or `sp.spmatrix`? | `X must be numpy array or sparse matrix` |
| X dims | 2D? | `X must be 2-dimensional, got shape ...` |
| obs length | Matches X rows? | `obs has N rows but X has M cells` |
| var length | Matches X columns? | `var has N rows but X has M genes` |
| obsm | All 2D with same N rows? | `obsm['X_pca'] has N rows but X has M cells` |
| layers | Same shape as X? | `layers['counts'] shape ... != X shape ...` |
| raw | Valid CanonicalObject? | `raw object is invalid: ...` |

---

## Error hierarchy

```
ScinteropError (base)
├── DetectionError          # Format detection failure
├── FormatAdapterError      # Base for all adapter errors
│   ├── MtxAdapterError     # MTX read/write failure
│   ├── H5adAdapterError    # H5AD read/write failure
│   └── RdsAdapterError     # RDS read/write failure
├── ValidationError         # Object validation failure
├── ExecError               # Base for execution errors
│   ├── RExecError          # R script failure
│   └── PythonExecError     # Python script failure
├── CacheError              # Scratch directory failure
└── ProvenanceError         # Provenance log failure
```

Every error includes the specific location and cause:

```python
try:
    si.read("nonexistent.h5ad")
except FormatAdapterError as e:
    print(type(e).__name__, ":", e)
    # H5adAdapterError : Failed to read H5AD file 'nonexistent.h5ad': [Errno 2] No such file or directory
```

---

## CLI reference

```bash
# Detect format of a file or directory
scinterop detect data.h5ad
scinterop detect cellranger_output/

# Convert formats
scinterop convert input.h5ad output.rds --seurat
scinterop convert input.h5ad output_mtx/
scinterop convert input.rds output.h5ad --r-exe /opt/R/4.3/bin/Rscript

# Convert options:
#   --seurat        When writing RDS, create a Seurat object (not plain R list)
#   --r-exe PATH    Path to Rscript (overrides SCINTEROP_R_EXE)
#   --python-exe    Path to Python (overrides SCINTEROP_PYTHON_EXE)
#   --debug         Keep temporary files for debugging

# Run an external script
scinterop run analysis.R --executor r --exe /opt/R/4.3/bin/Rscript
scinterop run analysis.py --executor python --conda-env sc-env
scinterop run analysis.py --executor python --exe /opt/python/3.11/bin/python
```

---

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `SCINTEROP_R_EXE` | `Rscript` | Path to Rscript executable |
| `SCINTEROP_PYTHON_EXE` | `python` | Path to Python executable |
| `SCINTEROP_CONDA_EXE` | `conda` | Conda/micromamba executable (for conda-run mode) |
| `SCINTEROP_SCRATCH` | `/tmp/scinterop_$USER` | Base directory for temporary files |
| `SCINTEROP_DEBUG` | (unset) | Set to `1`/`true`/`yes` to keep temp files on success |

---

## Smoke test

Run this to verify the entire package works end-to-end without any external dependencies (no R, no anndata):

```python
"""
Smoke test: run with:
  PYTHONPATH=/path/to/scinterop python smoke_test.py
"""
import scinterop as si
import numpy as np
import pandas as pd
from scipy import sparse as sp
import tempfile
import os

print("=" * 50)
print("scinterop smoke test")
print("=" * 50)

# 1. Create a CanonicalObject
np.random.seed(42)
X = sp.csr_matrix(np.random.poisson(0.5, size=(20, 10)).astype(np.float64))
obs = pd.DataFrame({"cluster": ["A", "B"] * 10}, index=[f"cell_{i}" for i in range(20)])
var = pd.DataFrame(index=[f"gene_{i}" for i in range(10)])
obsm = {"X_pca": np.random.randn(20, 3)}
obj = si.CanonicalObject(X=X, obs=obs, var=var, obsm=obsm)
print(f"[OK] Created CanonicalObject: {obj.shape}")

# 2. Validate
result = si.validate(obj)
assert result.valid, f"Validation failed: {result.errors}"
si.assert_valid(obj)
print("[OK] Validation passed")

# 3. Detect format
for path in ["test.h5ad", "test.rds", "matrix.mtx", "matrix.mtx.gz"]:
    d = si.detect(path)
    assert d.fmt, f"Detection failed for {path}"
print("[OK] Format detection works on all extensions")

# 4. MTX round-trip
tmpdir = tempfile.mkdtemp()
si.mtx.write_mtx(obj, tmpdir)
obj2 = si.mtx.read_mtx(tmpdir)
assert obj2.shape == (20, 10), f"Shape mismatch: {obj2.shape}"
assert np.allclose(obj.X.toarray(), obj2.X.toarray()), "X mismatch"
print("[OK] MTX write + read round-trip")

# 5. Top-level read (MTX)
obj3 = si.read(tmpdir)
assert obj3.shape == (20, 10)
print("[OK] si.read() on MTX directory")

# 6. MTX → MTX convert (via si.read + si.mtx.write, no actual conversion needed)
tmpdir2 = tempfile.mkdtemp()
si.convert(tmpdir, tmpdir2)
obj4 = si.read(tmpdir2)
assert obj4.shape == (20, 10)
print("[OK] si.convert() MTX → MTX")

# 7. Provenance
log_dir = tempfile.mkdtemp()
from scinterop.provenance import write_conversion_record, read_provenance_log
write_conversion_record(
    input_path="/smoke/test.h5ad",
    input_format="h5ad",
    output_path="/smoke/test.rds",
    output_format="rds",
    success=True,
    runtime_s=0.5,
    log_dir=log_dir,
)
records = read_provenance_log(os.path.join(log_dir, "provenance.jsonl"))
assert len(records) == 1
assert records[0]["success"] is True
print("[OK] Provenance logging")

# 8. Cache manager
from scinterop.cache import ScratchManager
mgr = ScratchManager(base=tempfile.mkdtemp())
ctx = mgr.tempdir("smoke")
assert ctx.path.exists()
mgr.cleanup(ctx, keep=False)
assert not ctx.path.exists()
print("[OK] Cache manager (create + cleanup)")

# 9. Object convenience methods
assert obj.to_anndata is not None
assert obj.to_seurat is not None
assert obj.to_mtx is not None
print("[OK] Object .to_* methods exist")

# 10. Copy
obj_copy = obj.copy()
assert obj_copy.shape == obj.shape
print("[OK] Object copy")

print()
print("=" * 50)
print("All smoke tests passed!")
print("=" * 50)
```

**To run with H5AD support:**
```bash
pip install anndata
PYTHONPATH=/path/to/scinterop python smoke_test.py
# Also run the H5AD-specific tests:
python -m pytest /path/to/scinterop/scinterop/tests/test_h5ad.py -v
```

**To run the full test suite:**
```bash
cd /path/to/scinterop
pip install -e .[dev]    # installs pytest
python -m pytest
```

---

## Extending with a new format

Adding support for a new format (e.g., Loom, h5seurat, Zarr) requires:

1. **Create the adapter module** with two functions:

```python
# scinterop/loom.py
from scinterop.schema import CanonicalObject
from scinterop.errors import FormatAdapterError

def read_loom(path):
    """Read Loom file -> CanonicalObject"""
    # ... parsing logic ...
    return CanonicalObject(X=..., obs=..., var=...)

def write_loom(obj, path):
    """CanonicalObject -> Loom file"""
    # ... writing logic ...
    return str(path)
```

2. **Register in `detect.py`**: add the extension to `EXTENSION_MAP`:

```python
EXTENSION_MAP = {
    ".h5ad": "h5ad",
    ".rds": "rds",
    ".mtx": "mtx",
    ".loom": "loom",          # new
}
```

3. **Register in `__init__.py`**: add the format to `read()` and `convert()`:

```python
def read(path, **kwargs):
    ...
    elif fmt == "loom":
        return loom.read_loom(path, **kwargs)
    ...
```

That's it. The adapter owns all format-specific logic; the rest of the package is format-agnostic.

---

## Citation

If you use scinterop in your work, please cite it:

```bibtex
@misc{scinterop,
  author = {Muhammad Dawood},
  title  = {scinterop: Single-cell data interoperability between R and Python formats},
  year   = {2026},
  publisher = {GitHub},
  url    = {https://github.com/imuhdawood/scinterop}
}
```

See [CITATION.cff](CITATION.cff) for the machine-readable citation metadata.
