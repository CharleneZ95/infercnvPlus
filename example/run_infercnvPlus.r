
library(infercnvPlus)

# Run example with built-in data
data(expr, genomic_pos, ref_obs, package = "infercnvPlus")

# Data tranforming: genes(rows) X cells(columns)
## For 10X counts data 
expr_tr <- umi_to_log2tpm(expr)

## For Smart-seq2 TPM values
expr_tr <- log2(expr + 1)

# Calucate cnv score
cnv_obj <- inferCNV(data = expr_tr,
                    gene_pos = genomic_pos,
                    cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                    reference_obs = ref_obs,
                    window_size = 101,
                    out_path = "output_dir", # dir is auto-created for storing outputs
                    noise_filter = NULL,
                    vis_bounds = "-1,1")

# Cluster cells and visualize
cnv_obj <- visualCNV(data = cnv_obj,
                     cutree_k = 2,
                     out_file = "plot_cnv.png")

# Extract cells from the specific subtrees
cnv_obj <- extractCells(data = cnv_obj,
                        subtrees = 2,
                        lab_to_rm = "ref")

# Get cell barcode
cells <- cnv_obj$target

