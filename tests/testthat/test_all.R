
context("All functions")

load('../testdata/test_data.rda')

test_that("inferCNV", {
	res1 <- inferCNV(data = npc_expr,
                      gene_pos = genomic_pos,
                      cutoff = 0.1,
                      reference_obs = ref_obs,
                      window_size=101,
                      contig_tail=NULL)
	expect_that(res1$cnv_score, equals(infercnv_obj$cnv_score))
})

test_that("denoiseVis", {
    res2 <- denoiseVis(data=infercnv_obj, noise_threshold = 0.2)
    expect_that(res2$cnv_score_vis_filter0.2, equals(infercnv_obj$cnv_score_vis_filter0.2))
})

test_that("visualCNV", {
    res3 <- visualCNV(data = infercnv_obj, cutree_k = 2)
    
    expect_that(res3$cluster$hclust, equals(infercnv_obj$cluster$hclust))
})


test_that("pruneCutree ", {
    res4 <- pruneCutree (data = infercnv_obj,
                         k = 3,
                         out_path = 'output_dir')
    
    expect_that(res4$cluster$cutree, equals(infercnv_obj$cluster$cutree))
})


test_that("extractCells", {
    res5 <- extractCells(data = infercnv_obj, subtrees = c(2,3))
    expect_that(res5$target, equals(infercnv_obj$target))
})
