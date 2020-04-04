
context("All functions")

load('../testdata/test_data.rda')


test_that("inferCNV", {
	expr_tr <- umi_to_log2tpm(expr)
    res1 <- inferCNV(data = expr_tr,
                     gene_pos = genomic_pos,
                     cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                     reference_obs = ref_obs,
                     window_size = 101,
                     out_path = "output_dir", # dir is auto-created for storing outputs
                     noise_filter = NULL,
                     vis_bounds = "-1,1")
    expect_that(res1$cnv_score_raw, equals(cnv_obj$cnv_score_raw))
})

test_that("denoiseVis", {
    res2 <- denoiseVis(data = cnv_obj, noise_threshold = 0.1)
    expect_that(res2$cnv_score_vis_filter0.1, equals(cnv_obj$cnv_score_vis_filter0.1))
})

test_that("visualCNV", {
    res3 <- visualCNV(data = cnv_obj)
    expect_that(res3$cluster$hclust, equals(cnv_obj$cluster$hclust))
})


test_that("pruneCutree ", {
    res4 <- pruneCutree(data = cnv_obj,
                        k = 2,
                        out_path = 'output_dir')
    expect_that(res4$cluster$cutree, equals(cnv_obj$cluster$cutree))
})


test_that("extractCells", {
    res5 <- extractCells(data = cnv_obj, subtrees = 2, lab_to_rm = 'ref')
    expect_that(res5$target, equals(cnv_obj$target))
})

