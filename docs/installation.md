# Installation

## From PyPI

```bash
pip install scinterop
```

## With optional dependencies

```bash
# H5AD support
pip install scinterop[anndata]

# Development (testing, docs)
pip install scinterop[dev]
```

## From source

```bash
git clone https://github.com/imuhdawood/scinterop.git
cd scinterop
pip install -e .[dev]
```

## R dependencies (required for RDS format)

scinterop bridges to R via `Rscript`. You need:

- R (>= 4.3)
- `anndata` R package
- `Seurat` R package (for Seurat objects)

```r
install.packages(c("anndata", "Seurat"))
```

## Conda

An `environment.yml` is provided for reproducible environments:

```bash
conda env create -f environment.yml
conda activate scinterop
```
