
library(infercnvPlus)

data(npc_expr, ref_obs, genomic_pos, package = 'infercnvPlus')

# Calucate cnv score
cnv_obj <- inferCNV(data = npc_expr,
                    gene_pos = genomic_pos,
                    cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                    reference_obs = ref_obs,
                    window_size = 101,
                    out_path = 'output_dir', # dir is auto-created for storing outputs
                    noise_filter = 0,
                    vis_bounds = '-1,1')

# Cluster cells and visualize
cnv_obj <- visualCNV(data = cnv_obj,
                     cutree_k = 2,
                     out_file = 'plot_cnv.png')

# Extract cells from the specific subtrees 
cnv_obj <- extractCells(data = cnv_obj, subtrees = c(2,3))

