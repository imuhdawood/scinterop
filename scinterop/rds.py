"""RDS/QS format adapter — read/write Seurat R objects via an R bridge.

Supports both ``.rds`` and ``.qs`` (qs serialization) file formats.
Requires R with the ``anndata`` and (for Seurat objects) ``Seurat``
R packages installed. For ``.qs`` files, the ``qs`` R package is also
required. Conversion is done by writing an intermediate H5AD file and
invoking R as a subprocess.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from .cache import ScratchManager
from .errors import RdsAdapterError
from .h5ad import read_h5ad, write_h5ad
from .r_runner import run_r, resolve_r_exe
from .schema import CanonicalObject
from .validate import assert_valid

logger = logging.getLogger(__name__)


_scratch = ScratchManager()


READ_R_SCRIPT_TEMPLATE = r"""{header}
# Read RDS/QS and save as H5AD via anndata R package

input_path <- {input_rds!r}
output_h5ad <- {output_h5ad!r}

# Load input file (RDS or QS)
if (grepl("\\.qs$", input_path, ignore.case = TRUE)) {{
    if (!requireNamespace("qs", quietly = TRUE)) {{
        stop("R package 'qs' is required for .qs files. Install with: install.packages('qs')")
    }}
    library(qs)
    obj <- qread(input_path)
}} else {{
    obj <- readRDS(input_path)
}}

# Determine if Seurat
is_seurat <- inherits(obj, "Seurat")

# Load anndata
suppressPackageStartupMessages({{
    if (!requireNamespace("anndata", quietly = TRUE)) {{
        stop("R package 'anndata' is required. Install with: install.packages('anndata')")
    }}
    library(anndata)
}})

if (is_seurat) {{
    # Seurat object -> H5AD
    if (!requireNamespace("Seurat", quietly = TRUE)) {{
        stop("R package 'Seurat' is required to read Seurat RDS files")
    }}
    library(Seurat)

    # Seurat v5+: use layer; Seurat v4: use slot
    if (packageVersion("SeuratObject") >= "5.0.0") {{
        counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
    }} else {{
        counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
    }}
    mat <- t(counts)

    meta <- obj[[]]
    meta_names <- rownames(meta)
    if (is.null(meta_names)) {{
        meta_names <- colnames(obj)
    }}

    gene_names <- rownames(counts)
    if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(nrow(counts)))

    adata <- AnnData(
        X = mat,
        obs = data.frame(row.names = meta_names),
        var = data.frame(row.names = gene_names)
    )

    for (col_name in colnames(meta)) {{
        if (is.numeric(meta[[col_name]]) || is.character(meta[[col_name]]) || is.factor(meta[[col_name]])) {{
            adata$obs[[col_name]] <- meta[[col_name]]
        }}
    }}

    reductions <- Reductions(obj)
    for (red in reductions) {{
        embed <- Embeddings(obj, red)
        if (!is.null(rownames(embed))) {{
            adata$obsm[[paste0("X_", tolower(red))]] <- as.matrix(embed)
        }}
    }}

    adata$write_h5ad(output_h5ad, compression = "gzip")
    cat("Seurat -> H5AD conversion complete:", output_h5ad, "\n")

}} else {{
    # Assume named list or supported structure
    supported <- c("X", "obs", "var", "obsm", "layers", "uns", "raw")
    if (!is.list(obj)) {{
        stop("RDS must be a Seurat object or a named list with components: ",
             paste(supported, collapse = ", "))
    }}

    missing <- setdiff(supported[1:3], names(obj))
    if (length(missing) > 0) {{
        stop("RDS list is missing required components: ", paste(missing, collapse = ", "))
    }}

    X <- obj[["X"]]
    obs_df <- obj[["obs"]]
    var_df <- obj[["var"]]

    if (is.null(rownames(obs_df))) {{
        rownames(obs_df) <- paste0("cell_", seq_len(nrow(obs_df)))
    }}
    if (is.null(rownames(var_df))) {{
        rownames(var_df) <- paste0("gene_", seq_len(nrow(var_df)))
    }}

    adata <- AnnData(
        X = X,
        obs = obs_df,
        var = var_df
    )

    if (!is.null(obj[["obsm"]])) {{
        for (key in names(obj[["obsm"]])) {{
            adata$obsm[[key]] <- obj[["obsm"]][[key]]
        }}
    }}

    if (!is.null(obj[["layers"]])) {{
        for (key in names(obj[["layers"]])) {{
            adata$layers[[key]] <- obj[["layers"]][[key]]
        }}
    }}

    if (!is.null(obj[["uns"]])) {{
        adata$uns <- obj[["uns"]]
    }}

    adata$write_h5ad(output_h5ad, compression = "gzip")
    cat("R list -> H5AD conversion complete:", output_h5ad, "\n")
}}
"""

WRITE_LIST_R_SCRIPT_TEMPLATE = r"""{header}
# Convert H5AD to R named list and save as RDS/QS

input_h5ad <- {input_h5ad!r}
output_path <- {output_rds!r}

suppressPackageStartupMessages({{
    if (!requireNamespace("anndata", quietly = TRUE)) {{
        stop("R package 'anndata' is required. Install with: install.packages('anndata')")
    }}
    library(anndata)
}})

adata <- read_h5ad(input_h5ad)

X <- adata$X
obs <- adata$obs
var <- adata$var

result <- list(
    X = X,
    obs = obs,
    var = var
)

obsm_names <- names(adata$obsm)
if (length(obsm_names) > 0) {{
    result$obsm <- list()
    for (key in obsm_names) {{
        result$obsm[[key]] <- adata$obsm[[key]]
    }}
}}

layer_names <- names(adata$layers)
if (length(layer_names) > 0) {{
    result$layers <- list()
    for (key in layer_names) {{
        result$layers[[key]] <- adata$layers[[key]]
    }}
}}

result$uns <- adata$uns

if (grepl("\\.qs$", output_path, ignore.case = TRUE)) {{
    if (!requireNamespace("qs", quietly = TRUE)) {{
        stop("R package 'qs' is required for .qs files. Install with: install.packages('qs')")
    }}
    library(qs)
    qsave(result, output_path)
}} else {{
    saveRDS(result, output_path)
}}
cat("R list saved to:", output_path, "\n")
"""

WRITE_SEURAT_R_SCRIPT_TEMPLATE = r"""{header}
# Convert H5AD to Seurat object and save as RDS/QS

input_h5ad <- {input_h5ad!r}
output_path <- {output_rds!r}

suppressPackageStartupMessages({{
    if (!requireNamespace("anndata", quietly = TRUE)) {{
        stop("R package 'anndata' is required. Install with: install.packages('anndata')")
    }}
    if (!requireNamespace("Seurat", quietly = TRUE)) {{
        stop("R package 'Seurat' is required. Install with: install.packages('Seurat')")
    }}
    if (!requireNamespace("Matrix", quietly = TRUE)) {{
        stop("R package 'Matrix' is required. Install with: install.packages('Matrix')")
    }}
    library(anndata)
    library(Seurat)
    library(Matrix)
}})

to_dgC <- function(x) {{
    if (inherits(x, "dgCMatrix")) return(x)
    if (inherits(x, "dgRMatrix")) return(as(as(x, "TsparseMatrix"), "dgCMatrix"))
    if (inherits(x, "Matrix")) return(as(x, "dgCMatrix"))
    if (is.matrix(x)) return(as(Matrix(x, sparse = TRUE), "dgCMatrix"))
    stop("Unsupported matrix type: ", class(x))
}}

adata <- read_h5ad(input_h5ad)

# Expression matrix: cells x genes -> genes x cells (keep sparse)
counts <- t(adata$X)
counts <- to_dgC(counts)

# Metadata
meta <- adata$obs
if (ncol(meta) == 0 || all(sapply(meta, class) == "list")) {{
    # Empty or all-list metadata; create minimal df
    meta <- data.frame(row.names = rownames(meta))
}}

# Set dimnames
gene_names <- rownames(adata$var)
if (is.null(gene_names)) gene_names <- paste0("gene_", seq_len(ncol(counts)))
cell_names <- rownames(adata$obs)
if (is.null(cell_names)) cell_names <- paste0("cell_", seq_len(ncol(counts)))

rownames(counts) <- gene_names
colnames(counts) <- cell_names
rownames(meta) <- cell_names

# Create Seurat object
seurat_obj <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    assay = "RNA"
)

# Add reductions from obsm
obsm_names <- names(adata$obsm)
for (key in obsm_names) {{
    red_name <- tolower(sub("^X_", "", key))
    key_prefix <- paste0(toupper(red_name), "_")

    mat <- as.matrix(adata$obsm[[key]])
    rownames(mat) <- cell_names
    colnames(mat) <- paste0(key_prefix, seq_len(ncol(mat)))

    seurat_obj[[red_name]] <- CreateDimReducObject(
        embeddings = mat,
        key = key_prefix,
        assay = DefaultAssay(seurat_obj)
    )
}}

if (grepl("\\.qs$", output_path, ignore.case = TRUE)) {{
    if (!requireNamespace("qs", quietly = TRUE)) {{
        stop("R package 'qs' is required for .qs files. Install with: install.packages('qs')")
    }}
    library(qs)
    qsave(seurat_obj, output_path)
}} else {{
    saveRDS(seurat_obj, output_path)
}}
cat("Seurat object saved to:", output_path, "\n")
"""


def read_rds(path: str | Path, *, r_exe: str | None = None) -> CanonicalObject:
    """Read an RDS or QS file (Seurat or plain R list) into a CanonicalObject.

    Uses an intermediate H5AD file written by the R ``anndata`` package.

    Args:
        path: Path to an ``.rds`` or ``.qs`` file.
        r_exe: Path to Rscript executable.

    Returns:
        A validated CanonicalObject.

    Raises:
        RdsAdapterError: If the file doesn't exist, R fails, or the
            intermediate H5AD is not produced.
    """
    path = Path(path)
    if not path.exists():
        raise RdsAdapterError(f"File does not exist: {path}")

    resolved_r = resolve_r_exe(r_exe)
    ctx = _scratch.tempdir(prefix="rds_read")

    try:
        temp_h5ad = ctx.path / "temp_rds_output.h5ad"
        header = _r_header("RDS Reader", resolved_r)

        r_script = READ_R_SCRIPT_TEMPLATE.format(
            header=header,
            input_rds=str(path.resolve()),
            output_h5ad=str(temp_h5ad.resolve()),
        )

        script_path = _scratch.write_script(ctx, r_script, "read_rds")
        run_r(script_path, r_exe=r_exe)

        if not temp_h5ad.exists():
            raise RdsAdapterError(
                "R bridge completed but no H5AD output was produced"
            )

        obj = read_h5ad(temp_h5ad)
        logger.info(
            "Read RDS file: %s -> CanonicalObject (%d cells x %d genes)",
            path, obj.n_cells, obj.n_genes,
        )
        return obj

    except Exception as e:
        raise RdsAdapterError(
            f"Failed to read RDS file '{path}': {e}"
        ) from e
    finally:
        _scratch.cleanup(ctx)


OUTPUT_EXTENSIONS: set[str] = {".rds", ".qs"}


def write_rds(
    obj: CanonicalObject,
    path: str | Path,
    *,
    seurat: bool = False,
    r_exe: str | None = None,
) -> str:
    """Write a CanonicalObject to an RDS or QS file.

    Uses an intermediate H5AD file and the R ``anndata`` package.

    Args:
        obj: The object to write.
        path: Output ``.rds`` or ``.qs`` file path (extension added if missing).
        seurat: If True, create a Seurat object; otherwise a plain R list.
        r_exe: Path to Rscript executable.

    Returns:
        The resolved path the file was written to.

    Raises:
        RdsAdapterError: If writing the intermediate H5AD fails, the R
            script fails, or the output is not produced.
    """
    assert_valid(obj)
    path = Path(path)
    if path.suffix.lower() not in OUTPUT_EXTENSIONS:
        path = path.with_suffix(".rds")

    resolved_r = resolve_r_exe(r_exe)
    ctx = _scratch.tempdir(prefix="rds_write")

    try:
        temp_h5ad = ctx.path / "temp_input.h5ad"
        write_h5ad(obj, temp_h5ad)

        header = _r_header("RDS Writer", resolved_r)

        if seurat:
            r_script = WRITE_SEURAT_R_SCRIPT_TEMPLATE.format(
                header=header,
                input_h5ad=str(temp_h5ad.resolve()),
                output_rds=str(path.resolve()),
            )
            label = "Seurat"
        else:
            r_script = WRITE_LIST_R_SCRIPT_TEMPLATE.format(
                header=header,
                input_h5ad=str(temp_h5ad.resolve()),
                output_rds=str(path.resolve()),
            )
            label = "R list"

        script_path = _scratch.write_script(ctx, r_script, "write_rds")
        run_r(script_path, r_exe=r_exe)

        if not path.exists():
            raise RdsAdapterError(
                "R bridge completed but no RDS output was produced"
            )

        logger.info(
            "Wrote %s RDS file: %s (%d cells x %d genes)",
            label, path, obj.n_cells, obj.n_genes,
        )
        return str(path)

    except Exception as e:
        raise RdsAdapterError(
            f"Failed to write RDS file '{path}': {e}"
        ) from e
    finally:
        _scratch.cleanup(ctx)


def _r_header(task: str, r_exe: str) -> str:
    return f"""# Auto-generated by scinterop RDS adapter ({task})
# R executable: {r_exe}
# Generated for debugging purposes
options(warn = 1)
"""
