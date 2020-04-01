#!/usr/bin/env Rscript

logging::basicConfig(level = "INFO")

# opposite to %in%
"%!in%" <- function(x, y) !(x %in% y)

# Return the indices of the rows that average above the cut off
#
# Args:
# data: Data to measure the average row and evaluate
#                 against the cutoff. Row = Genes, Col = Cells.
# cutoff: Threshold to be above to be kept.
#
# Returns:
# Returns a vector of row indicies to keep (are above the cutoff).
#
above_cutoff <- function(data, cutoff) {
    
    average_gene <- log2(rowMeans(((2^data) - 1), na.rm = TRUE) + 1)
    logging::loginfo(":Averages (counts).")
    # Find averages above a certain threshold
    genes <- rownames(data)[which(average_gene > cutoff)]
    if (length(genes) > 0) {
        return(genes)
    } else {
        return(NULL)
    }
}


# Center data and threshold (both negative and postive values)
#
# Args:
# center_data: Matrix to center. Row = Genes, Col = Cells.
# threshold: Values will be required to be with -/+1 *
#                      threshold after centering.
# Returns:
# Centered and thresholded matrix.
#
center_with_threshold <- function(center_data) {
    
    # Center data (automatically ignores zeros)
    center_data <- center_data - rowMeans(center_data, na.rm = TRUE)
    # Cap values between threshold and -threshold and recenter
    center_data[center_data > 3] <- 3
    center_data[center_data < -3] <- -3
    center_data <- center_data - rowMeans(center_data, na.rm = TRUE)
    return(center_data)
}


# Smooth a matrix by column using a simple moving average.
# Tails of the averages use a window length that is truncated to
# available data.
#
# Args:
# data: Data matrix to smooth. Row = Genes, Col = Cells.
# window_size: Length of window to use for the moving average.
#        Should be a positive, odd integer.
#
# Returns:
# Matrix with columns smoothed with a simple moving average.
#
smooth_window <- function(data, window_size) {
    
    if (window_size < 2) 
        return(data)
    if (window_size > nrow(data)) 
        return(data)
    
    tail_length <- (window_size - 1)/2
    num_genes <- nrow(data)
    data_sm <- apply(data, 2, smooth_window_helper, window_size = window_size)
    # Fix ends
    data_end <- apply(data, 2, smooth_ends_helper, obs_tails = tail_length)
    
    for (row_end in 1:tail_length) {
        end_bound <- (num_genes - row_end) + 1
        data_sm[row_end, ] <- data_end[row_end, ]
        data_sm[end_bound, ] <- data_end[end_bound, ]
    }
    
    # Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)
    return(data_sm)
}

# Helper function for smoothing the ends of a moving average.
#
# Args:
# obs_data: Data to smooth
# obs_tails: Length of the tail to smooth.
#
# Returns:
# Data smoothed.
#
smooth_ends_helper <- function(obs_data, obs_tails) {
    end_data <- rep(NA, length(obs_data))
    obs_count <- length(obs_data)
    for (tail_end in 1:obs_tails) {
        bounds <- tail_end - 1
        end_tail <- obs_count - bounds
        end_data[tail_end] <- mean(obs_data[(tail_end - bounds):(tail_end + bounds)], 
            na.rm = TRUE)
        end_data[end_tail] <- mean(obs_data[(end_tail - bounds):(end_tail + bounds)], 
            na.rm = TRUE)
    }
    return(end_data)
}

# Smooth vector of values over the given window length.
#
# Args:
# obs_data: Vector of data to smooth with a moving average.
# window_size: Length of the window for smoothing.
#        Must be and odd, positive, integer.
#
# Returns:
# Vector of values smoothed with a moving average.
#
smooth_window_helper <- function(obs_data, window_size) {
    return(stats::filter(obs_data, rep(1/window_size, window_size), sides = 2))
}

# Center data after smoothing. Center with in cells using median.
#
# Args:
# data_smoothed: Matrix to center.
#                          Row = Genes, Col = cells.
#
# Returns:
# Matrix that is median centered.
#             Row = Genes, Col = cells.
#
center_smoothed <- function(data_smoothed) {
    
    logging::loginfo(":Center_smoothed.")
    # Center within columns
    row_median <- apply(data_smoothed, 2, median)
    return(t(apply(data_smoothed, 1, "-", row_median)))
}


# Remove the average of the genes of the reference observations from all
# observations' expression. Normalization by column.
#
# Args:
# average_data: Matrix containing the data to remove average
#               from (this includes the reference observations).
#               Row = Genes, Col = Cells.
# ref_observations: Indices of reference observations.
#                   Only these are used in the average.
# ref_groups: A list of vectors of indices refering to the
#             different groups of the reference indices.
#
# Returns:
# Expression with the average gene expression in the reference
#          observations removed.
#
average_over_ref <- function(average_data, ref_observations) {
    # r = genes, c = cells Max and min mean gene expression within reference groups.
    ref_average <- rowMeans(average_data[, ref_observations, drop = FALSE])
    # Remove the average gene expression of the reference groups from the For each
    # gene.
    average_data <- sweep(average_data, 1, ref_average, "-")
    return(average_data)
}


# Remove the tails of values of a specific chromosome.
# The smooth_matrix values are expected to be in genomic order.
# If the tail is too large and no contig will be left 1/3 of the
# contig is left.
#
# Args:
# smooth_matrix: Smoothed values in genomic order.
#                          Row = Genes, Col = Cells.
# chr: Indices of the chr in which the tails are to be removed.
# tail_length: Length of the tail to remove on both ends of the
#                        chr indices.
# Returns:
# Indices to remove.
#
remove_tails <- function(smooth_matrix, chr, tail_length) {
    
    chr_length <- length(chr)
    if ((tail_length < 3) || (chr_length < 3)) {
        return(c())
    }
    if (chr_length < (tail_length * 2)) {
        tail_length <- floor(chr_length/3)
    }
    remove_indices <- -1 * chr[1:tail_length]
    remove_indices <- c(remove_indices, -1 * (chr[((chr_length + 1) - tail_length):chr_length]))
    return(remove_indices)
}


# Order the data and subset the data to data in the genomic position file.
#
# Args:
# data: Data (expression) matrix where the row names should be in
#                 the row names of the genomic_position file.
# genomic_position: Data frame read in from the genomic position file
#
# Returns:
# A matrix of expression in the order of the
# genomic_position file. NULL is returned if the genes in both
# data parameters do not match.
# 
order_reduce <- function(data, genomic_position) {
    
    ret_results <- list(expr = NULL, order = NULL, chr_order = NULL)
    if (is.null(data) || is.null(genomic_position)) {
        return(ret_results)
    }
    
    # Drop pos_gen entries that are position 0
    remove_by_position <- -1 * which(genomic_position[2] + genomic_position[3] == 
        0)
    if (length(remove_by_position)) {
        genomic_position <- genomic_position[remove_by_position, , drop = FALSE]
    }
    
    # Reduce to genes in pos file
    keep_genes <- row.names(data)[which(row.names(data) %in% row.names(genomic_position))]
    
    # Keep genes found in position file
    if (length(keep_genes)) {
        ret_results$expr <- data[keep_genes, , drop = FALSE]
        ret_results$order <- genomic_position[keep_genes, , drop = FALSE]
    } else {
        logging::loginfo(paste("inferCNV:order_reduce: The position file ", "and the expression file row (gene) names do not match."))
        return(list(expr = NULL, order = NULL, chr_order = NULL))
    }
    
    # Set the chr to factor so the order can be arbitrarily set and sorted.
    chr_levels <- unique(genomic_position[["chr"]])
    ret_results$order[["chr"]] <- factor(ret_results$order[["chr"]], levels = chr_levels)
    
    # Sort genomic position file and expression file to genomic position file Order
    # genes by genomic region
    order_names <- row.names(ret_results$order)[with(ret_results$order, order(chr, 
        start, stop))]
    ret_results$expr <- ret_results$expr[order_names, , drop = FALSE]
    
    # This is the contig order, will be used in visualization.  Get the contig order
    # in the same order as the genes.
    ret_results$order <- ret_results$order[order_names, , drop = FALSE]
    ret_results$chr_order <- ret_results$order[1]
    
    # Remove any gene without position information Genes may be sorted correctly by
    # not have position information Here they are removed.
    logging::loginfo(paste0(":Order and subset genes, ", "new dimensions (genes,cells) = ", 
        paste(dim(data), collapse = ","), "."))
    return(ret_results)
}


#' @title Calculate CNV scores given a matrix of RNASeq counts. Output a matrix of final values.
#'
#' @param data: expression matrix (genes X samples),
#'                 assumed to be log2(TPM+1) .
#' @param gene_order: ordering of the genes (data's rows)
#'                       according to their genomic location
#'                       To include all genes use 0.
#' @param cutoff: cut-off for the average expression of genes to be
#'                   used for CNV inference.
#' @param reference_obs: Column names of the subset of samples (data's columns)
#'                          that should be used as references.
#'                          If not given, the average of all samples will
#'                          be the reference.
#' @param window_size: length of the window for the moving average
#'                          (smoothing). Should be an odd integer.
#' @param out_path: the path to what to save the pdf as. The raw data is
#'                     also written to this path but with the extension .txt .
#' @param contig_tail: length of the tail removed from the ends of contigs.
#' @param noise_filter: the minimum difference a value can be from the
#'                            average reference in order for it not to be
#'                            removed as noise.
#' @param vis_bounds: Used as upper and lower bounds for values in the visualization.
#'                      Should be given in the form of '-1,1' (lower bound, upper bound).
#'
#' @return
#' Returns an 'infercnv' object including:
#'     CNV matrix before visualization.
#'     CNV matrix after denoise and outlier removal for visualization.
#'     Chromosome order
#'     Names of cells in reference groups.
#' 
calcCNV <- function(data,
                    gene_pos,
                    cutoff,
                    reference_obs,
                    window_size,
                    out_path,
                    contig_tail,
                    noise_filter,
                    vis_bounds) {
   
    ret_list <- list()
    
    # Data
    data <- as.matrix(data)
    logging::loginfo(paste(":Original matrix dimensions (genes,cells) =", paste(dim(data), 
        collapse = ",")))
    
    # Order and reduce the expression to the genomic file.
    colnames(gene_pos) <- c("chr", "start", "stop")
    order_ret <- order_reduce(data = data, genomic_position = gene_pos)
    data <- order_ret$expr
    gene_order <- order_ret$order
    if (is.null(data)) {
        error_message <- paste("None of the genes in the expression data", "matched the genes in the reference genomic", 
            "position file. Analysis Stopped.")
        stop(error_message)
    }
    
    # Reduce by cutoff
    keep_genes <- above_cutoff(data, cutoff)
    if (!is.null(keep_genes)) {
        data <- data[keep_genes, , drop = FALSE]
        gene_order <- gene_order[keep_genes, , drop = FALSE]
        logging::loginfo(paste0(":Reduce by cutoff, ", "new dimensions (genes,cells) = ", 
            paste(dim(data), collapse = ","), "."))
        logging::logdebug(":Keeping genes.")
    } else {
        logging::loginfo(":Reduce by cutoff.")
        logging::logwarn(paste("::No genes left to keep.", "Stoping."))
        stop(998)
    }
    
    # Reduce contig info
    chr_order <- gene_order[1]
    gene_order <- NULL
    
    # Center data (automatically ignores zeros)
    data <- center_with_threshold(data)
    logging::loginfo(paste0(":Outlier removal, ", " Min=", round(min(data), 3), " Max=", 
        round(max(data), 3), "."))
    
    # Smooth the data with gene windows
    data_smoothed <- smooth_window(data, window_size)
    data <- NULL
    logging::loginfo(":Smoothed data.")
    
    # Center cells/observations after smoothing. This helps reduce the effect of
    # complexity.
    data_smoothed <- center_smoothed(data_smoothed)
    
    # Remove average reference
    data_smoothed <- average_over_ref(average_data = data_smoothed, ref_observations = reference_obs)
    
    logging::loginfo(paste0(":Remove average, ", " Min=", round(min(data_smoothed), 3), 
        " Max=", round(max(data_smoothed), 3), "."))
    
    # Remove Ends
    remove_indices <- c()
    logging::loginfo(":Remove tails.")
    for (chr in unlist(unique(chr_order))) {
        logging::loginfo(paste0("::Remove tail contig ", chr, "."))
        remove_chr <- remove_tails(data_smoothed, which(chr_order == chr), contig_tail)
        remove_indices <- c(remove_indices, remove_chr)
    }
    if (length(remove_indices) > 0) {
        chr_order <- chr_order[remove_indices, ]
        data_smoothed <- data_smoothed[remove_indices, , drop = FALSE]
    }

    
    # Save objects
    ret_list[["cnv_score_raw"]] <- as.matrix(data_smoothed)
    ret_list[["chr_order"]] <- paste(as.vector(as.matrix(chr_order)))
    ret_list[["reference_obs"]] <- reference_obs
    ret_list[["parameter"]] <- list(gene_cutoff = cutoff, window_size = window_size, 
        contig_tail = contig_tail)
    class(ret_list) <- "infercnv"

    # Denoise and outliers removal for visualization
    ret_list <- denoiseVis(data = ret_list, 
                           noise_threshold = noise_filter, 
                           vis_bounds = vis_bounds)

    assay <- ifelse(is.null(noise_filter), "cnv_score_vis",
        paste0("cnv_score_vis_filter", noise_filter))

    finalMat <- ret_list[[assay]]
    
    # Log output
    logging::loginfo(paste0(":Final dimensions (genes,cells) = ", paste(dim(finalMat), 
        collapse = ","), " Min=", round(min(finalMat), 3), " Max=", round(max(finalMat), 
        3), "."))
    
    return(ret_list)
}


#' Clustering and plotting cells based on cnv score matrix using ComplexHeatmap package.
#'
#' Args:
#' @param data: an 'infercnv' object as produced by inferCNV.
#' @param assay: assay which assigned to use for plotting.
#' @param ref_lab: label for reference cells.
#' @param obs_lab: label for observations.
#' @param dist_method: distance measure used in clustering cells, possible
#'                      values are 'correlation' for Pearson correlation 
#'                      and all the distances supported by stats::dist.
#' @param clustering_method: clustering method used. Accepts the same values as stats::hclust.
#' @param cutree_k: an integer scalar or vector with the desired number of groups
#' @param colors: vector of colors used in heatmap.
#' @param border: whether draw border. The value can be logical or a string of color.
#' @param out_file: filename to save plot
#' @param out_path: output directory
#'
#' @return An 'infercnv' object 
#' 
plotCNV <- function(data,
                    assay,
                    ref_lab,
                    obs_lab,
                    dist_method,
                    clustering_method,
                    cutree_k,
                    colors,
                    border,
                    plot_dend,
                    out_file, 
                    out_path) {
    
    # data processing -----------------
    logging::loginfo(":Data processing")
    ret_list <- data
    ret_list[["cluster"]] <- list()
    
    plot_data <- as.matrix(t(data[[assay]]))
    hc <- hclust(dist(plot_data, method = dist_method), method = clustering_method)
    ret_list[["cluster"]][["hclust"]] <- hc
    
    dend <- as.dendrogram(hc)
    if (is.numeric(cutree_k)) {
        ret_list <- pruneCutree(ret_list, k = cutree_k, out_path = out_path, plot_dend=plot_dend)
        dend <- ret_list[["cluster"]][["dend"]]
        ret_list[["cluster"]][["dend"]] <- NULL
    }
    
    # plot ----------------------
    logging::loginfo(":Plot cnv")

    # set annotation column annotation: chromosome
    chr_order <- data$chr_order
    chr <- ifelse(chr_order %in% c(seq(1, 22, 2), "X"), "odd", "even")
    col_anno <- HeatmapAnnotation(chr = chr, 
        col = list(chr = c(odd = "gray60", even = "gray80")), 
        show_legend = FALSE, show_annotation_name = TRUE)

    ## row annotation: origin
    ref_obs <- as.character(data$reference_obs)
    origin <- rep(obs_lab, nrow(plot_data))
    origin[which(rownames(plot_data) %in% ref_obs)] <- ref_lab
    row_anno <- rowAnnotation(origin=origin, col = list(origin = c('obs'='#FFFF00', 
        'ref'='#696969')), annotation_legend_param = list(title_gp=gpar(fontsize=12, fontface='bold'), 
        labels_gp=gpar(fontsize=11), legend_height = unit(4, "cm"), grid_width=unit(0.5, "cm")), 
        show_legend=TRUE, show_annotation_name=FALSE)
   
    # set legend parameters
    vis_bounds <- data[["parameter"]][["vis_bounds"]]
    limits <- c(min(vis_bounds), 0, max(vis_bounds))
    
    lg_param <- list(title='cnv',
                     at=limits,
                     labels=limits,
                     title_gp=gpar(fontsize=12, fontface='bold'),
                     labels_gp=gpar(fontsize=11),
                     legend_height = unit(3, "cm"),
                     grid_width=unit(0.5, "cm"),
                     border='black')
    # plot
    png(file.path(out_path, out_file), width = 9.5, height = 6, units = "in", res = 300)
    
    ht <- Heatmap(plot_data,
                  name = "cnv",
                  col = colors,
                  cluster_columns = FALSE,
                  cluster_rows = dend,
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  row_split = cutree_k,
                  heatmap_legend_param = lg_param,
                  top_annotation = col_anno,
                  left_annotation = row_anno,
                  row_dend_width = unit(0.15, "npc"),
                  row_gap = unit(0, "mm"),
                  border=TRUE)
    draw(ht)
    
    # add chromosome labels
    chr_lab <- c()
    x_axis <- c()
    for (c in unique(chr_order)) {
        perc <- length(which(chr_order == c))/length(chr_order)
        chr_lab <- c(chr_lab, ifelse(perc >= 0.02, c, "."))
        x_axis <- c(x_axis, median(which(chr_order == c))/length(chr_order))
    }

    decorate_annotation("chr", {
        grid.text(chr_lab, x = unit(x_axis, "npc"), gp=gpar(fontsize=10.5))
    })

    dev.off()
    return(ret_list)
}


# Subtract cells from specific branches 
#
# Args: 
# data: an 'infercnv' object as produced by inferCNV.
# subtrees: specify the numbers of subtrees to extract. 
# lab_to_rm: a character scalar with the label of cells to exclude.
#
# Returns: 
# An 'infercnv' object with target cells
#
extract_helper <- function(data, 
                           subtrees = NULL, 
                           lab_to_rm = "ref") {
    # extract leaves
    dendlist <- data[["cluster"]][["cutree"]]
    dendl <- dendlist[as.numeric(subtrees)]
    leaves <- unlist(sapply(dendl, labels))
    cells <- grep(lab_to_rm, leaves, value = TRUE, invert = TRUE)
    data[["target"]] <- cells
    data[["parameter"]][["subtrees"]] <- subtrees
    return(data)
    
}
