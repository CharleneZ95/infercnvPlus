#!/usr/bin/env Rscript

#' @import logging
#' @importFrom stats filter
#' @importFrom magrittr "%>%"
#' @importFrom dendextend set
#' @importFrom dendextend color_branches
#' @importFrom dendextend cutree prune
#' @importFrom grid grid.text gpar unit
#' @importFrom ComplexHeatmap draw
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom ComplexHeatmap decorate_annotation
#' @importFrom Seurat GetAssayData
#' 
#' @title Infer CNV changes given a matrix of RNASeq counts. Output a pdf and matrix of final values.
#'
#' @param data expression matrix (genes X cells), assumed to be log2(TPM+1).
#' @param gene_pos ordering of the genes (data's rows)
#'                       according to their genomic location
#'                       To include all genes use 0.
#' @param cutoff cut-off for the average expression of genes to be
#'                   used for CNV inference.
#' @param reference_obs column names of the subset of cells (data's columns)
#'                          that should be used as references.
#'                          If not given, the average of all cells will
#'                          be the reference.
#' @param window_size length of the window for the moving average
#'                          (smoothing). Should be an odd integer.
#' @param out_path the path to what to save the pdf as. The raw data is
#'                     also written to this path but with the extension .txt .
#' @param contig_tail length of the tail removed from the ends of contigs.
#' @param noise_filter the minimum difference a value can be from the
#'                            average reference in order for it not to be
#'                            removed as noise.
#' @param vis_bounds Used as upper and lower bounds for values in the visualization.
#'                      Should be given in the form of '-1,1' (lower bound, upper bound).
#'
#' @return
#' Returns an 'infercnv' object including:
#'     1. CNV matrix before visualization.
#'     2. CNV matrix after denoise and outlier removal for visualization.
#'     3. Chromosome order.
#'     4. Names of cells in reference groups.
#' 
#' @export
#' 
inferCNV <- function(data = NULL,
                     gene_pos = NULL,
                     cutoff = NULL,
                     reference_obs,
                     window_size = 101,
                     out_path = 'output_dir',
                     contig_tail = NULL,
                     noise_filter = NULL,
                     vis_bounds = "-1,1") {
    
    logging::loginfo("inferCNV: Start.")
    
    # Check arguments ----------------
    logging::loginfo(":Check_arguments.")
    
    # check data
    if (is.null(data)) {
        stop("data: Expression matrix (genes X cells), assumed to be log2(TPM+1).")
    }
    
    # genomic pisition file
    if (is.null(gene_pos)) {
        stop("gene_pos: A genomic position file is needed.")
    }
    
    
    # Require the cut off to be above 0
    if (is.null(cutoff) || cutoff < 0) {
        stop("cutoff: Please enter a value", "greater or equal to zero for the cut off.")
    }
    
    # Warn that an average of the cells is used in the absence of normal /
    # reference cells
    if (is.null(reference_obs)) {
        logging::logwarn(paste0("reference_obs: No reference cells were given, the average of the cells ", 
            "will be used."))
        reference_obs <- colnames(data)
    }
    
    # Make sure the given reference cells are in the matrix.
    if (length(setdiff(reference_obs, colnames(data))) > 0) {
        missing_reference_cell <- setdiff(reference_obs, colnames(data))
        error_message <- paste0("Please make sure that all the reference cell ", 
            "names match a cell in your data matrix. ", "Attention to:", paste(missing_reference_cell, 
                collapse = ","), '.')
        stop(error_message)
    }
    
    # Require the name of output path
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    
    # check contig_tail
    if (is.null(contig_tail)) {
        logging::logwarn(paste("contig_tail: contig_tail wasn't", "assigned, (window_size - 1) / 2", 
            "will be used."))
        contig_tail = (window_size - 1)/2
    }
    
    # check noise filter
    if (is.null(noise_filter) || noise_filter == 0) {
        noise_filter <- NULL
    }
    
    # Run infercnv
    ret_list <- calcCNV(data = data, 
                        gene_pos = gene_pos, 
                        cutoff = cutoff, 
                        reference_obs = reference_obs, 
                        window_size = window_size, 
                        out_path = out_path, 
                        contig_tail = contig_tail, 
                        noise_filter = noise_filter, 
                        vis_bounds = vis_bounds)
    
    return(ret_list)
}


#' @import RColorBrewer
#' 
#' 
#' @title Clustering and plotting cells based on cnv score matrix using ComplexHeatmap package.
#'
#' Args:
#' @param data an 'infercnv' object as produced by inferCNV.
#' @param assay assay which assigned to use for plotting.
#' @param ref_lab label for reference cells.
#' @param obs_lab label for observations.
#' @param dist_method distance measure used in clustering cells, possible
#'                      values are 'correlation' for Pearson correlation 
#'                      and all the distances supported by stats::dist.
#' @param clustering_method clustering method used. Accepts the same values as stats::hclust.
#' @param cutree_k an integer scalar or vector with the desired number of groups.
#' @param colors vector of colors used in heatmap.
#' @param border whether draw border. The value can be logical or a string of color.
#' @param out_file filename to save plot.
#' @param out_path output directory.
#' @param plot_dend whether to plot dendrogram or not. The value should be logical.
#' 
#' @return 
#' Returns an 'infercnv' object including:
#'     1. Clustering result.
#'     2. Dendlist of 'cutree'.
#' 
#' @export
#' 
visualCNV <- function(data,
                      assay = NULL,
                      ref_lab = "ref",
                      obs_lab = "obs",
                      distance_method = 'euclidean',
                      clustering_method = "ward.D2",
                      cutree_k = NULL,
                      colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11),
                      border = TRUE,
                      plot_dend = FALSE,
                      out_file = 'plot_cnv.png',
                      out_path = 'output_dir') {
    
    # Check arguments --------------------
    logging::loginfo("visualCNV: Start.")
    
    # Require an 'infercnv' object
    if (class(data) != "infercnv") {
        stop("Require an infercnv object.")
    }
    
    # check assay
    if (is.null(assay)) {
        noise_filter <- data$parameter$noise_filter
        assay <- ifelse(is.null(noise_filter), 'cnv_score_vis',
            paste0("cnv_score_vis_filter", noise_filter))
        logging::logwarn(paste0("assay: Since assay was not assigned, \"", assay, 
            "\" will be used."))
    }
    
    # Agglomeration method to be used.
    METHODS <- c("ward.D", "single", "complete", "average", "mcquitty", "median", 
        "centroid", "ward.D2")
    if (clustering_method %!in% METHODS) {
        error_message <- paste("clustering_method: This should be (an unambiguous", 
            "abbreviation of) one of \"ward.D\",", "\"ward.D2\", \"single\", \"complete\",", 
            "\"average\" (= UPGMA), \"mcquitty\" (= WPGMA),", "\"median\" (= WPGMC) or \"centroid\" (= UPGMC).")
        stop(error_message)
    }
    
    # Whether to cutree or not
    if (!is.null(cutree_k) & !is.numeric(cutree_k)) {
        stop("cutree_k: Numeric scala is required.")
    }
    
    # Output path
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    
    ret_list <- plotCNV(data = data, 
                        assay = assay, 
                        ref_lab = ref_lab, 
                        obs_lab = obs_lab, 
                        dist_method = distance_method, 
                        clustering_method = clustering_method, 
                        cutree_k = cutree_k, 
                        colors = colors, 
                        border = border,
                        plot_dend = plot_dend, 
                        out_file = out_file, 
                        out_path = out_path)
    
    return(ret_list)
}


#' Remove values that are too close to the average and are considered noise.
#'
#' Args:
#' @param smooth_matrix a matrix of values, smoothed, and with average
#'                          reference removed. Row = Genes, Col = Cells.
#' @param noise_threshold the minimum difference a value can be from the
#'                            average reference in order for it not to be
#'                            removed as noise.
#' @param vis_bounds used as upper and lower bounds for values in the visualization.
#'                      Should be given in the form of '-1,1' (lower bound, upper bound).
#'
#' @return 
#' Returns an 'infercnv' object including: denoised matrix.
#'
#' @export
#' 
denoiseVis <- function(data = NULL, 
                       noise_threshold = NULL, 
                       vis_bounds = '-1,1') {
    
    if (is.null(data) || class(data) != "infercnv") {
        stop("data: Data should be an \"infercnv\" object.")
    }
    
    bounds <- strsplit(vis_bounds, ',')[[1]]
    lower_bound <- as.numeric(bounds[1])
    upper_bound <- as.numeric(bounds[2])
    if (!is.numeric(lower_bound) || !is.numeric(upper_bound)) {
        stop("Should be given in the form of \'-1,1\'.")
    }
   
    logging::loginfo("Remove noise.")
    ret_list <- data
    data <- data$cnv_score_raw
    
    # Noise filter
    if (!is.null(noise_threshold)) {
        data[abs(data) < noise_threshold] <- 0
    }
    
    # Handle outliers
    logging::loginfo("Remove outlier and normalize data.")
    ## Hard threshold given bounds
    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound
    
    # return data
    label <- ifelse(is.null(noise_threshold), "cnv_score_vis",
        paste0("cnv_score_vis_filter", noise_threshold))
    ret_list[[label]] <- data
    ret_list[["parameter"]][["noise_filter"]] <- noise_threshold
    ret_list[["parameter"]][["vis_bounds"]] <- c(lower_bound, upper_bound)
    return(ret_list)
}


#' @title Prune cutree to dendlist
#'
#' Args:
#' @param data an 'infercnv' object as produced by inferCNV.
#' @param k an integer scalar or vector with the desired number of groups.
#' @param out_path output directory.
#' @param plot_dend whether to plot dendrogram or not. The value should be logical.
#'
#' @return Return an 'infercnv' object with dendlist.
#'
#' @export
#'
pruneCutree <- function(data, 
                        k, 
                        out_path = "output_dir",
                        plot_dend = FALSE) {
    # data processing
    dend <- data[["cluster"]][["hclust"]] %>% 
        as.dendrogram() %>% 
        set("labels_to_character") %>% 
        color_branches(k = k, groupLabels = TRUE)

    # Plot dendrogram with groupLabels for each subtrees
    if (plot_dend) {
        png(file.path(out_path, "plot_dendrogram.png"), width = 7, height = 10, unit = "in", 
            res = 300)
        plot(dend, main = "Original dendrogram", horiz = TRUE)
        dev.off()
    }

    # Prune cutree to dendlist
    clusters <- dendextend::cutree(dend, k, order_clusters_as_data = FALSE)
    split <- data.frame(clusters, paste0("C", as.numeric(clusters)))
    dendl <- vector("list", k)
    for (i in 1:k) {
        leves_to_prune <- labels(dend)[clusters != i]
        dendl[[i]] <- prune(dend, leves_to_prune)
    }
    
    class(dendl) <- "dendlist"
    data[["cluster"]][["cutree"]] <- dendl
    data[["parameter"]][["cutree_k"]] <- k
    data[["cluster"]][["dend"]] <- dend
    return(data)
}


#' Subtract cells from specific branches
#'
#' Args:
#' @param data an 'infercnv' object as produced by inferCNV.
#' @param subtrees a numeric scalar or vector specify those subtrees to extract. 
#' @param lab_to_rm a character scalar with the label of cells to exclude.
#'
#' @return Return an 'infercnv' object with target cells.
#'
#' @export
#'
extractCells <- function(data, 
                         subtrees = NULL, 
                         lab_to_rm = "ref") {
    
    # data should be a infercnv object
    if (class(data) != "infercnv") {
        stop("data: should be an \"infercnv\" object.")
    }
    
    # subtrees should be numeric
    if (!is.numeric(subtrees)) {
        stop("subtrees: should be numeric.")
    }
    
    ret_list <- extract_helper(data, subtrees, lab_to_rm)
    return(ret_list)
}


#' Transform UMI to log2(TPM+1)
#' 
#' Args:
#' @param data UMI counts matrix (genes X cells).
#' 
#' @return Return a log2-transformed TPM matrix. 
#' 
#' @export 
#' 
umi_to_log2tpm <- function(data) {
    data <- log2( 1e6 * ( sweep( as.matrix(data), 2, colSums(data), '/' ) ) + 1)
    return(data)
}



#' Import a seurat object and convert it to matrix for inferCNV.
#' 
#' Args:
#' @param data a Seurat object.
#' 
#' @return Return a TPM log2-transformed matrix.
#' 
#' @export 
#' 
importSrat <- function(obj,
                       slot = 'counts',
                       assay = 'RNA') {
    rawdata <- GetAssayData(object=obj, slot = slot, assay = assay)
    data <- umi_to_log2tpm( as.matrix(rawdata) )
    return(data)
}
