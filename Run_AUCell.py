# AUCell functions for scRNA-seq gene set enrichment
import pandas as pd
from pyscenic.aucell import aucell, derive_auc_threshold
from ctxcore.genesig import GeneSignature

# Configuration
config = {
    "adata_path": f"{dataPath}/adata.h5ad",
    "gmt_file": f"{dataPath}/cell_type_markers_set.gmt"
}

# Load AnnData
import scanpy as sc
adata = sc.read_h5ad(config["adata_path"])

# Load gene signatures from GMT
signatures = GeneSignature.from_gmt(
    config["gmt_file"],
    field_separator="\t",
    gene_separator="\t"
)
print(f"Loaded {len(signatures)} gene signatures.")

# Prepare sparse count matrix
sparse_mtx = pd.DataFrame.sparse.from_spmatrix(
    adata.layers["counts"],
    index=adata.obs_names,
    columns=adata.var["gene_symbols"]
)

# Derive AUC threshold (uses default percentile thresholds)
percentiles = derive_auc_threshold(sparse_mtx)

# Run AUCell to compute enrichment scores
aucs_mtx = aucell(
    sparse_mtx,
    signatures,
    auc_threshold=percentiles[1],  # Using the 95th percentile by default
    num_workers=8
)

# Merge AUCell scores into AnnData observations
adata.obs = adata.obs.merge(aucs_mtx, left_index=True, right_index=True)

# Save the updated AnnData object
adata.write_h5ad("Aucell_adata.h5ad")