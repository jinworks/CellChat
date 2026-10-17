#' Normalize data using a scaling factor
#'
#' @param data.raw input raw data
#' @param scale.factor the scaling factor used for each cell
#' @param do.log whether to do log transformation with pseudocount 1
#' @param do.sparse whether to use sparse format
#' @export
#'
normalizeData <- function(data.raw, scale.factor = 10000, do.log = TRUE, do.sparse = TRUE) {
  # Scale counts within a sample
  library.size <- Matrix::colSums(data.raw)
  #scale.factor <- median(library.size)
  expr <- Matrix::t(Matrix::t(data.raw) / library.size) * scale.factor
  if (do.log) {
    data.norm <-log1p(expr)
  }
  if (do.sparse) {
    data.input <- as(data.norm, "dgCMatrix")
  }
  return(data.norm)
}


#' Scale the data
#'
#' @param data.use input data
#' @param do.center whether center the values
#' @export
#'
scaleData <- function(data.use, do.center = T) {
  data.use <- Matrix::t(scale(Matrix::t(data.use), center = do.center, scale = TRUE))
  return(data.use)
}


#' Scale a data matrix
#'
#' @param x data matrix
#' @param scale the method to scale the data
#' @param na.rm whether remove na
#' @importFrom Matrix rowMeans colMeans rowSums colSums
#' @return
#' @export
#'
#' @examples
scaleMat <- function(x, scale, na.rm=TRUE){

  av <- c("none", "row", "column", 'r1', 'c1')
  i <- pmatch(scale, av)
  if(is.na(i) )
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  scale <- av[i]

  switch(scale, none = x
         , row = {
           x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), '-',check.margin = FALSE)
           sx <- apply(x, 1L, sd, na.rm = na.rm)
           sweep(x, 1L, sx, "/", check.margin = FALSE)
         }
         , column = {
           x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), '-',check.margin = FALSE)
           sx <- apply(x, 2L, sd, na.rm = na.rm)
           sweep(x, 2L, sx, "/", check.margin = FALSE)
         }
         , r1 = sweep(x, 1L, rowSums(x, na.rm = na.rm), '/', check.margin = FALSE)
         , c1 = sweep(x, 2L, colSums(x, na.rm = na.rm), '/', check.margin = FALSE)
  )
}

#' Downsampling single cell data using geometric sketching algorithm
#'
#' USERs need to install the python package `pip install geosketch` (https://github.com/brianhie/geosketch)
#'
#' @param object A data matrix (should have row names; samples in rows, features in columns) or a Seurat object.
#'
#' When object is a PCA or UMAP space, please set `do.PCA = FALSE`
#'
#' When object is a data matrix (cells in rows and genes in columns), it is better to use the highly variable genes. PCA will be done on this input data matrix.
#' @param percent the percent of data to sketch
#' @param idents A vector of identity classes to keep for sketching
#' @param do.PCA whether doing PCA on the input data
#' @param dimPC the number of components to use
#' @importFrom reticulate import
#' @return A vector of cell names to use for downsampling
#' @export
#'
sketchData <- function(object, percent, idents = NULL, do.PCA = TRUE, dimPC = 30) {
  # pip install geosketch
  geosketch <- reticulate::import('geosketch')
  if (is(object,"Seurat")) {
    sketch.size <- as.integer(percent*ncol(object))
    if (!is.null(idents)) {
      object <- subset(object, idents = idents)
    }
    object <- object %>% #Seurat::NormalizeData(verbose = FALSE) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      RunPCA(pc.genes = object@var.genes, npcs = dimPC, verbose = FALSE)

    X.pcs <- object@reductions$pca@cell.embeddings
    cells.all <- Cells(object)

  } else {
    # Get top PCs
    if (do.PCA) {
      X.pcs <- runPCA(object, dimPC = dimPC)
    } else {
      X.pcs <- object
    }

    # Sketch percent of data.
    sketch.size <- as.integer(percent*nrow(X))
    cells.all <- rownames(object)
  }
  sketch.index <- geosketch$gs(X.pcs, sketch.size)
  sketch.index <- unlist(sketch.index) + 1
  sketch.cells <- cells.all[sketch.index]
  return(sketch.cells)
}


#' Add the cell information into meta slot
#'
#' @param object CellChat object
#' @param meta cell information to be added
#' @param meta.name the name of column to be assigned
#'
#' @return
#' @export
#'
#' @examples
addMeta <- function(object, meta, meta.name = NULL) {
  if (is.null(x = meta.name) && is.atomic(x = meta)) {
    stop("'meta.name' must be provided for atomic meta types (eg. vectors)")
  }
  if (inherits(x = meta, what = c("matrix", "Matrix"))) {
    meta <- as.data.frame(x = meta)
  }

  if (is.null(x = meta.name)) {
    meta.name <- names(meta)
  } else {
    names(meta) <- meta.name
  }
  object@meta <- meta
  return(object)
}


#' Set the default identity of cells
#' @param object CellChat object
#' @param ident.use the name of the variable in object.meta;
#' @param levels set the levels of factor
#' @param display.warning whether display the warning message
#' @return
#' @export
#'
#' @examples
setIdent <- function(object, ident.use = NULL, levels = NULL, display.warning = TRUE){
  if (!is.null(ident.use)) {
    object@idents <- as.factor(object@meta[[ident.use]])
  }

  if (!is.null(levels)) {
    object@idents <- factor(object@idents, levels = levels)
  }
  if ("0" %in% as.character(object@idents)) {
    stop("Cell labels cannot contain `0`! ")
  }
  if (length(object@net) > 0) {
    if (all(dimnames(object@net$prob)[[1]] %in% levels(object@idents) )) {
      message("Reorder cell groups! ")
      cat("The cell group order before reordering is ", dimnames(object@net$prob)[[1]],'\n')
      # idx <- match(dimnames(object@net$prob)[[1]], levels(object@idents))
      idx <- match(levels(object@idents), dimnames(object@net$prob)[[1]])
      object@net$prob <- object@net$prob[idx, , ]
      object@net$prob <- object@net$prob[, idx, ]
      object@net$pval <- object@net$pval[idx, , ]
      object@net$pval <- object@net$pval[, idx, ]
      cat("The cell group order after reordering is ", dimnames(object@net$prob)[[1]],'\n')
    } else {
      message("Rename cell groups but do not change the order! ")
      cat("The cell group order before renaming is ", dimnames(object@net$prob)[[1]],'\n')
      dimnames(object@net$prob) <- list(levels(object@idents), levels(object@idents), dimnames(object@net$prob)[[3]])
      dimnames(object@net$pval) <- dimnames(object@net$prob)
      cat("The cell group order after renaming is ", dimnames(object@net$prob)[[1]],'\n')
    }
    if (display.warning) {
      warning("All the calculations after `computeCommunProb` should be re-run!!
    These include but not limited to `computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`.")
    }


  }
  return(object)
}


#' Add a reduced space of the data into CellChat object
#'
#' @param object CellChat object from a single dataset
#' @param dr A data frame (rows are cells with rownames) consisting of a low-dimensional space for visualization
#' @param dr.name A char name of the reduction method for the input `dr`
#' @param seu.obj A Seurat object with the reduced space of the data
#' @param dr.use A char name of the reduction method to use when taking `seu.obj` as input. By default, all reduced space in `seu.obj` will be added in `object@dr`
#' @param force.add Whether to force to add a new reduced space when a reduced space exists in `object@dr`
#' @return
#' @export
#' @examples
#' \dontrun{
#' cellChat <- addReduction(object = cellchat, dr = cell.embeddings, dr.name = "umap")
#'
#' cellChat <- addReduction(object = cellchat, seu.obj = seu.obj)
#' }
addReduction <- function(object, dr = NULL, dr.name = NULL, seu.obj = NULL, dr.use = NULL, force.add = FALSE) {
  if (length(names(object@dr)) > 0) {
    if (!force.add) {
      stop(paste0("The `object@dr` contains the following reduced space: ", toString(names(object@dr)), ". Please set `force.add = TRUE` if intending to add a new reduced space. \n"))
    }
  }
  if (!is.null(dr)) {
    if (is.null(dr.name)) {
      stop("When inputing `dr`, please also provide the `dr.name`! \n")
    }
    dr <- as.data.frame(dr)
    if (all(colnames(object@data.signaling) %in% rownames(dr))) {
      object@dr[[dr.name]] <- dr[colnames(object@data.signaling), ]
    } else {
      stop("Some cell barcodes in the CellChat object are not the rownames of the input `dr`. Please check the input `dr` and make sure it contains all cells in the CellChat analysis. \n")
    }
  } else if(!is.null(seu.obj)) {
    if (!is(seu.obj,"Seurat")) {
      stop("The input `seu.obj` can be only the Seurat object. \n")
    }
    reductions <- names(seu.obj@reductions)
    if (length(reductions) == 0) {
      stop("The input `seu.obj` does not contain any low-dimensional space. Please generate a low-dimensional space for visualization. \n")
    }
    if (!is.null(dr.use)) {
      reductions <- intersect(reductions, dr.use)
    }
    if (length(reductions) == 0) {
      stop("The input `dr.use` is not in the reduced space in `seu.obj`. \n")
    }
    for (i in 1:length(reductions)) {
      dr.name <- reductions[i]
      dr = seu.obj@reductions[[dr.name]]@cell.embeddings
      if (all(colnames(object@data.signaling) %in% rownames(dr))) {
        cat(paste0(dr.name, " is now added in `object@dr` as a low-dimensional space. \n"))
        object@dr[[dr.name]] <- dr[colnames(object@data.signaling), ]
      } else {
        stop("Some cell barcodes in the CellChat object are not in the input `seu.obj`. Please check the input `seu.obj` and make sure it contains all cells in the CellChat analysis. \n")
      }
    }
  } else {
    stop("Please input either `dr` or `seu.obj`! \n")
  }
  return(object)
}


#' Update and re-order the cell group names after running `computeCommunProb`
#'
#' @param object CellChat object
#' @param old.cluster.name A vector defining old cell group labels in `object@idents`; Default = NULL, which will use `levels(object@idents)`
#' @param new.cluster.name A vector defining new cell group labels to rename
#' @param new.order reset order of cell group labels
#' @param new.cluster.metaname assign a name of the new labels, which will be the column name of new labels in `object@meta`
#' @return An updated CellChat object
#' @export
#'
updateClusterLabels <- function(object, old.cluster.name = NULL, new.cluster.name = NULL, new.order = NULL, new.cluster.metaname = "new.labels") {
  if (is.null(old.cluster.name)) {
    old.cluster.name <- levels(object@idents)
  }
  if (new.cluster.metaname %in% colnames(object@meta)) {
    stop("Please define another `new.cluster.metaname` as it exists in `colnames(object@meta)`!")
  }
  if (!is.null(new.cluster.name)) {
    labels.new <- plyr::mapvalues(object@idents, from = old.cluster.name, to = new.cluster.name)
    object@meta[[new.cluster.metaname]] <- labels.new
    object <- setIdent(object, ident.use = new.cluster.metaname, display.warning = FALSE)
  } else {
    new.cluster.metaname <- NULL
    cat("Only reorder cell groups but do not rename cell groups!")
  }

  if (!is.null(new.order)) {
    object <- setIdent(object, ident.use = new.cluster.metaname, levels = new.order, display.warning = FALSE)
  }
  message("We now re-run computeCommunProbPathway`,`aggregateNet`, and `netAnalysis_computeCentrality`...")
  object <- computeCommunProbPathway(object)
  ## calculate the aggregated network by counting the number of links or summarizing the communication probability
  object <- aggregateNet(object)
  # network importance analysis
  object <-netAnalysis_computeCentrality(object, slot.name = "netP")
  return(object)
}





#' Subset the expression data of signaling genes for saving computation cost
#'
#' @param object CellChat object
#' @param features default = NULL: subset the expression data of signaling genes in CellChatDB.use
#'
#' @return An updated CellChat object by assigning a subset of the data into the slot `data.signaling`
#' @export
#'
subsetData <- function(object, features = NULL) {
  interaction_input <- object@DB$interaction
  if (object@options$datatype != "RNA") {
    if ("annotation" %in% colnames(interaction_input) == FALSE) {
      warning("A column named `annotation` is required in `object@DB$interaction` when running CellChat on spatial transcriptomics! The `annotation` column is now automatically added and all L-R pairs are assigned as `Secreted Signaling`, which means that these L-R pairs are assumed to mediate diffusion-based cellular communication.")
      interaction_input$annotation <- "Secreted Signaling"
    }
  }
  if ("annotation" %in% colnames(interaction_input) == TRUE) {
    if (length(unique(interaction_input$annotation)) > 1) {
      interaction_input$annotation <- factor(interaction_input$annotation, levels = c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"))
      interaction_input <- interaction_input[order(interaction_input$annotation), , drop = FALSE]
      interaction_input$annotation <- as.character(interaction_input$annotation)
    }
    object@DB$interaction <- interaction_input
  }

  if (is.null(features)) {
    DB <- object@DB
    gene.use_input <- extractGene(DB)
    gene.use <- intersect(gene.use_input, rownames(object@data))
  } else {
    gene.use <- intersect(features, rownames(object@data))
  }
  object@data.signaling <- object@data[rownames(object@data) %in% gene.use, ]
  return(object)
}



#' Identify over-expressed signaling genes associated with each cell group
#'
#' USERS can use customized gene set as over-expressed signaling genes by setting `object@var.features[[features.name]] <- features.sig`
#' The adjusted p value (see `p.adjust.method`) can be obtained via `object@var.features[[paste0(features.name, ".info")]]`. Note that by default `features.name = "features"`
#'
#' @param object CellChat object
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the slot 'data.signaling' is used
#' @param group.by cell group information; default is `object@idents`; otherwise it should be one of the column names of the meta slot
#' @param idents.use a subset of cell groups used for analysis
#' @param invert whether to invert the idents.use
#' @param group.dataset dataset origin information in a merged CellChat object; set it as one of the column names of meta slot when identifying the highly enriched genes in one dataset for each cell group
#' @param pos.dataset the dataset name used for identifying highly enriched genes in this dataset for each cell group
#' @param group.DE.combined Whether to perform differential expression between conditions by ignoring cell group information. By default, group.DE.combined = FALSE, which will perform differential expression analysis between two biological conditions for each cell group;
#' When group.DE.combined = TRUE, it will perform DE analysis by combining all cell groups together.
#'
#' @param features.name a char name used for storing the over-expressed signaling genes in `object@var.features[[features.name]]`
#' @param only.pos Only return positive markers
#' @param features features used for identifying Over Expressed genes. default use all features
#' @param return.object whether to return the object; otherwise return a data frame consisting of over-expressed signaling genes associated with each cell group
#' @param neg.dataset the dataset name(s) compared to `pos.dataset`. Default = NULL: all the other datasets are pooled together.
#' When analyzing three or more conditions, setting e.g. pos.dataset = "TrtA" and neg.dataset = "Ctrl" performs the comparison TrtA vs. Ctrl only, ignoring the cells from the other conditions.
#' @param thresh.pc Threshold of the fraction of cells expressed in one cluster, i.e., thresh.pc = 0.1
#' @param thresh.fc Threshold of Log Fold Change, i.e., thresh.fc = 0.1. Increase it (e.g., thresh.fc = 0.25) to remove the genes that are only slightly different between groups or conditions
#' @param thresh.p Threshold of p-values, i.e., thresh.p = 0.05. NB: this threshold is applied to the raw (non-adjusted) p-values; see `thresh.p.adj`
#' @param thresh.p.adj Threshold of the p-values adjusted for multiple testing (column `pvalues.adj`), i.e., thresh.p.adj = 0.05. Default = NULL (no filtering on adjusted p-values, as in the previous versions).
#' Using thresh.p.adj is recommended: with thousands of genes tested in each cell group, filtering on raw p-values only keeps many false positives.
#' @param p.adjust.method the method used to adjust the p-values for multiple testing within each cell group (see `stats::p.adjust`); Default = "BH" (false discovery rate).
#' NB: when comparing conditions, cells are considered as independent observations (no biological replicates), which makes the p-values over-optimistic; effect size thresholds (`thresh.fc`, `thresh.pc`) are therefore important.
#' @param group.replicate the column name of `object@meta` giving the biological replicates (samples), e.g., the sample column of the Seurat metadata such as "orig.ident".
#' When set together with `group.dataset`, the differential expression between conditions is performed at the level of the replicates (pseudobulk) instead of the cells:
#' for each cell group, the expression of each replicate is the log of the average normalized expression of its cells (see `min.cells.replicate`), and the replicates are compared between conditions.
#' This takes into account the biological variability between replicates and avoids the over-optimistic p-values of the cell-level tests. At least two replicates per condition are required for each cell group.
#' logFC is the difference of the mean pseudobulk expression between conditions (natural log), and pct.1/pct.2 are the fractions of expressing cells in `pos.dataset` and in the other condition(s).
#' @param replicate.method the test used with `group.replicate`: "limma" (moderated t-test from the limma package; recommended for few replicates) or "t.test" (Welch t-test, or paired t-test when `paired.replicate = TRUE`)
#' @param paired.replicate whether the replicate labels identify the same individuals/donors in the different conditions (paired design, e.g. before/after treatment). Default = FALSE: replicates are considered as different samples in each condition, even if they share the same label (e.g., "rep1" in each condition).
#' @param min.cells.replicate the minimum number of cells of a replicate in a cell group for this replicate to be used in the pseudobulk analysis of this cell group. Default = 10
#' @param min.pct.test with `group.replicate`: the genes expressed in less than this fraction of the cells of a cell group (all the compared cells pooled, whatever their condition) are not tested in this cell group
#' and are not counted in the multiple-testing adjustment, e.g., min.pct.test = 0.05. This independent filtering (Bourgon et al., PNAS 2010) does not use the condition labels, hence it keeps the control of the false discovery rate
#' while increasing the power: the adjustment no longer pays for the genes that are hardly expressed in a cell group, which cannot be detected as differentially expressed and can hardly be involved in an inferred communication.
#' Default = NULL: all the genes are counted in the adjustment (conservative). NB: `thresh.pc` uses the fraction of expressing cells within each condition and the genes it removes are still counted in the adjustment.
#' @param covariates column names of `object@meta` used as covariates in the pseudobulk analysis (requires `group.replicate` and replicate.method = "limma"), e.g., covariates = "sex".
#' The covariates should be constant within each replicate. The model is then `~ covariates + condition`: the conditions are compared within each level of the covariates (e.g., males with males and females with females),
#' which removes the variability due to the covariates, and logFC is the effect of the condition adjusted for the covariates.
#' @param do.DE Whether to perform differential expression analysis. By default do.DE = TRUE; When do.DE = FALSE, selecting over-expressed genes that are expressed in more than `min.cells` cells.
#' @param do.fast If do.fast = TRUE, then perform a ultra-fast Wilcoxon test using presto package; otherwise using stats package. These two methods produce different logFC values, and the presto::wilcoxauc method gives smaller values.
#' @param min.cells the minmum number of expressed cells required for the genes that are considered for cell-cell communication analysis
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom stats sd wilcox.test p.adjust
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, two new elements named 'features.name' and paste0(features.name, ".info") will be added into the list `object@var.features`
#' `object@var.features[[features.name]]` is a vector consisting of the identified over-expressed signaling genes;
#' `object@var.features[[paste0(features.name, ".info")]]` is a data frame returned from the differential expression analysis
#' @export
#'
identifyOverExpressedGenes <- function(object, data.use = NULL, group.by = NULL, idents.use = NULL, invert = FALSE,
                                       group.dataset = NULL, pos.dataset = NULL, neg.dataset = NULL, group.DE.combined = FALSE,
                                       features.name = "features",  only.pos = TRUE, features = NULL, return.object = TRUE,
                                       thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05, thresh.p.adj = NULL, p.adjust.method = "BH",
                                       group.replicate = NULL, replicate.method = c("limma", "t.test"), paired.replicate = FALSE, min.cells.replicate = 10, min.pct.test = NULL, covariates = NULL,
                                       do.DE = TRUE, do.fast = TRUE, min.cells = 10) {
  p.adjust.method <- match.arg(p.adjust.method, stats::p.adjust.methods)
  if (!is.null(min.pct.test) && is.null(group.replicate)) {
    warning("`min.pct.test` is only used with `group.replicate` (pseudobulk analysis); it is ignored.", call. = FALSE)
  }
  if (!is.list(object@var.features)) {
    stop("Please update your CellChat object via `updateCellChat()`")
  }
  if (is.null(data.use)) {
    X <- object@data.signaling
    if (nrow(X) < 3) {stop("Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks OK.")}
  } else {
    X <- data.use
  }

  if (is.null(features)) {
    features.use <- row.names(X)
  } else {
    features.use <- intersect(features, row.names(X))
  }
  data.use <- X[features.use, , drop = FALSE]

  if (do.DE) {
    # select genes based on differential expression
    data.use <- as.matrix(data.use)

    if (is.null(group.by)) {
      labels <- object@idents
      if (!is.factor(labels)) {
        message("Use the joint cell labels from the merged CellChat object")
        labels <- object@idents$joint
      }
    } else {
      labels <- object@meta[[group.by]]
    }
    if (!is.factor(labels)) {
      labels <- factor(labels)
    }
    level.use <- levels(labels)[levels(labels) %in% unique(labels)]
    if (!is.null(idents.use)) {
      if (invert) {
        level.use <- level.use[!(level.use %in% idents.use)]
      } else {
        level.use <- level.use[level.use %in% idents.use]
      }
    }
    numCluster <- length(level.use)

    if (!is.null(group.dataset)) {
      labels.dataset <- as.character(object@meta[[group.dataset]])
      if (!is.null(group.replicate)) {
        if (!(group.replicate %in% colnames(object@meta))) {
          stop("`group.replicate` should be one of the column names of `object@meta`, e.g., the sample column of the Seurat metadata ('orig.ident', 'sample', ...)")
        }
        labels.replicate <- as.character(object@meta[[group.replicate]])
        if (!is.null(covariates)) {
          if (!all(covariates %in% colnames(object@meta))) {
            stop("`covariates` should be column names of `object@meta`: ", toString(setdiff(covariates, colnames(object@meta))), " not found.")
          }
          covariates.data <- object@meta[, covariates, drop = FALSE]
        } else {
          covariates.data <- NULL
        }
      } else if (!is.null(covariates)) {
        stop("`covariates` can only be used with `group.replicate` (pseudobulk analysis)")
      }
      if (!(pos.dataset %in% unique(labels.dataset))) {
        cat("Please set pos.dataset to be one of the following dataset names: ", unique(as.character(labels.dataset)))
        stop()
      }
      if (!is.null(neg.dataset)) {
        if (!all(neg.dataset %in% unique(labels.dataset)) || pos.dataset %in% neg.dataset) {
          stop("Please set neg.dataset to be one or several of the following dataset names (different from pos.dataset): ", toString(setdiff(unique(labels.dataset), pos.dataset)))
        }
        # only keep the cells from the two compared conditions
        cells.keep <- which(labels.dataset %in% c(pos.dataset, neg.dataset))
        data.use <- data.use[, cells.keep, drop = FALSE]
        labels <- droplevels(labels[cells.keep])
        labels.dataset <- labels.dataset[cells.keep]
        if (!is.null(group.replicate)) {
          labels.replicate <- labels.replicate[cells.keep]
          if (!is.null(covariates.data)) {
            covariates.data <- covariates.data[cells.keep, , drop = FALSE]
          }
        }
        level.use <- level.use[level.use %in% unique(labels)]
        numCluster <- length(level.use)
      }
      labels.dataset[labels.dataset != pos.dataset] <- toString(setdiff(unique(labels.dataset), pos.dataset))
      labels.dataset <- factor(labels.dataset, levels = c(pos.dataset, setdiff(unique(labels.dataset), pos.dataset)))
      if (nlevels(labels.dataset) < 2) {
        stop("No cells from other datasets than `pos.dataset` are available for the comparison!")
      }
    }

    if (!is.null(group.replicate)) {
      if (is.null(group.dataset)) {
        stop("`group.replicate` is used to compare conditions: please also set `group.dataset` and `pos.dataset`")
      }
      replicate.method <- match.arg(replicate.method)
      clusters.use <- labels
      level.use.pb <- level.use
      if (group.DE.combined) {
        # ignore the cell group information
        clusters.use <- factor(rep("all", length(labels)))
        level.use.pb <- "all"
      }
      markers.all <- pseudobulkDE(data.use, clusters = clusters.use, level.use = level.use.pb,
                                  conditions = labels.dataset, pos.dataset = pos.dataset, replicates = labels.replicate,
                                  method = replicate.method, paired = paired.replicate, min.cells = min.cells.replicate, covariates = covariates.data,
                                  thresh.pc = thresh.pc, min.pct.test = min.pct.test, p.adjust.method = p.adjust.method, n.features = length(features.use))
      if (group.DE.combined && nrow(markers.all) > 0) {
        markers.all <- do.call(rbind, lapply(level.use, function(x) {markers.all$clusters <- x; markers.all}))
      }
      if (nrow(markers.all) > 0) {
        markers.all <- markers.all[markers.all$pvalues < thresh.p & abs(markers.all$logFC) >= thresh.fc, , drop = FALSE]
        if (!is.null(thresh.p.adj)) {
          markers.all <- markers.all[markers.all$pvalues.adj < thresh.p.adj, , drop = FALSE]
        }
        if (only.pos) {
          markers.all <- markers.all[markers.all$logFC > 0, , drop = FALSE]
        }
      }
      markers.all$datasets <- factor(ifelse(markers.all$logFC > 0, pos.dataset, levels(labels.dataset)[2]), levels = levels(labels.dataset))
      markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), , drop = FALSE]
      rownames(markers.all) <- NULL
      object@var.features[[features.name]] <- unique(markers.all$features)
      object@var.features[[paste0(features.name, ".info")]] <- markers.all

    } else if (do.fast) {
      presto.check <- rlang::is_installed(c("presto"))
      if (!presto.check) {
        stop(
          "For a faster implementation of the Wilcoxon Test, please install the presto package",
          "\n--------------------------------------------",
          "\n devtools::install_github('immunogenomics/presto')",
          "\n--------------------------------------------",
          "\n Otherwise, plase set `do.fast = FALSE` for running the standard Wilcoxon Test!\n"
        )
      }
      rename.presto <- function(x, group.name) {
        colnames(x) <- plyr::mapvalues(colnames(x),from = c("group","feature","pval","logFC","pct_in","pct_out","padj"), to = c(group.name,"features","pvalues","logFC","pct.1", "pct.2","pvalues.adj"), warn_missing = FALSE)
        x
      }
      if (is.null(group.dataset)) {
        genes.de <- presto::wilcoxauc(data.use, labels, groups_use = level.use)
        genes.de <- rename.presto(genes.de, "clusters")
        # adjust the p-values within each cell group
        genes.de$pvalues.adj <- stats::ave(genes.de$pvalues, genes.de$clusters, FUN = function(x) stats::p.adjust(x, method = p.adjust.method))
      } else if ((!is.null(group.dataset)) & (group.DE.combined == FALSE)) {
        genes.de <- data.frame()
        for (i in 1:numCluster) {
          idx <- which(labels == level.use[i])
          labels.i <- droplevels(labels.dataset[idx])
          if (nlevels(labels.i) < 2) {
            next
          }
          genes.de.i <- presto::wilcoxauc(data.use[ ,idx, drop = FALSE], labels.i)
          # presto returns the results of each dataset vs. the others; keep `pos.dataset` vs. the other dataset(s) only,
          # otherwise each gene is reported twice and the genes up-regulated in the other dataset get a positive logFC
          genes.de.i <- genes.de.i[genes.de.i$group == pos.dataset, , drop = FALSE]
          genes.de.i$clusters <- level.use[i]
          genes.de.i$padj <- stats::p.adjust(genes.de.i$pval, method = p.adjust.method)
          genes.de <- rbind(genes.de, genes.de.i)
        }
        genes.de <- rename.presto(genes.de, "datasets")
      } else if ((!is.null(group.dataset)) & (group.DE.combined == TRUE)) {
        genes.de.c <- presto::wilcoxauc(data.use, labels.dataset)
        genes.de.c <- genes.de.c[genes.de.c$group == pos.dataset, , drop = FALSE]
        genes.de.c$padj <- stats::p.adjust(genes.de.c$pval, method = p.adjust.method)
        genes.de <- data.frame()
        for (i in 1:numCluster) {
          genes.de.c$clusters <- level.use[i]
          genes.de <- rbind(genes.de, genes.de.c)
        }
        genes.de <- rename.presto(genes.de, "datasets")
      }
      if (nrow(genes.de) == 0) {
        stop("No genes can be tested. Please check the cell groups and the datasets used for the comparison!")
      }
      if (!is.null(group.dataset)) {
        # direction of the change: `pos.dataset` when the gene is up-regulated in `pos.dataset`, otherwise the other dataset(s)
        genes.de$datasets <- ifelse(genes.de$logFC > 0, pos.dataset, levels(labels.dataset)[2])
      }
      genes.de$logFC_abs <- abs(genes.de$logFC)
      pct.max <- apply(genes.de[, c("pct.1", "pct.2"), drop = FALSE], MARGIN = 1, FUN = max)
      genes.de$pct.max <- pct.max
      markers.all <- dplyr::filter(genes.de, pvalues < thresh.p, logFC_abs >= thresh.fc, pct.max > thresh.pc*100)
      if (!is.null(thresh.p.adj)) {
        markers.all <- dplyr::filter(markers.all, pvalues.adj < thresh.p.adj)
      }
      if (!is.null(group.dataset)) {
        markers.all$datasets <- factor(markers.all$datasets, levels = levels(labels.dataset))
        markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), ]
      } else {
        markers.all <- markers.all[order(markers.all$pvalues), ]
      }
      markers.all <- dplyr::select(markers.all, -dplyr::any_of(c("logFC_abs","statistic","pct.max","group")))

      if (only.pos & nrow(markers.all) > 0) {
        markers.all <- subset(markers.all, subset = logFC > 0)
      }
      markers.all$features <- as.character(markers.all$features)

      features.sig <- markers.all$features
      object@var.features[[features.name]] <- unique(features.sig)
      features.name <- paste0(features.name, ".info")
      object@var.features[[features.name]] <- markers.all

    } else {

      my.sapply <- ifelse(
        test = future::nbrOfWorkers() == 1,
        yes = pbapply::pbsapply,
        no = future.apply::future_sapply
      )

      mean.fxn <- function(x) {
        return(log(x = mean(x = expm1(x = x)) + 1))
      }
      labels <- as.character(labels)
      genes.de <- vector("list", length = numCluster)
      for (i in 1:numCluster) {
        features <- features.use
        if (is.null(group.dataset)) {
          cell.use1 <- which(labels == level.use[i])
          cell.use2 <- base::setdiff(1:length(labels), cell.use1)
        } else if ((!is.null(group.dataset)) & (group.DE.combined == FALSE)) {
          cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
          cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
        } else if ((!is.null(group.dataset)) & (group.DE.combined == TRUE)) {
          cell.use1 <- which(labels.dataset == pos.dataset)
          cell.use2 <- which(labels.dataset != pos.dataset)
        }
        if (length(cell.use1) == 0 || length(cell.use2) == 0) {
          next
        }

        # feature selection (based on percentages)
        thresh.min <- 0
        pct.1 <- round(
          x = rowSums(data.use[features, cell.use1, drop = FALSE] > thresh.min) /
            length(x = cell.use1),
          digits = 3
        )
        pct.2 <- round(
          x = rowSums(data.use[features, cell.use2, drop = FALSE] > thresh.min) /
            length(x = cell.use2),
          digits = 3
        )
        data.alpha <- cbind(pct.1, pct.2)
        colnames(x = data.alpha) <- c("pct.1", "pct.2")
        alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
        names(x = alpha.min) <- rownames(x = data.alpha)
        features <- names(x = which(x = alpha.min > thresh.pc))
        if (length(x = features) == 0) {
          #stop("No features pass thresh.pc threshold")
          next
        }

        # feature selection (based on average difference)
        data.1 <- apply(X = data.use[features, cell.use1, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
        data.2 <- apply(X = data.use[features, cell.use2, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
        FC <- (data.1 - data.2)
        if (only.pos) {
          features.diff <- names(which(FC > thresh.fc))
        } else {
          features.diff <- names(which(abs(FC) > thresh.fc))
        }

        features <- intersect(x = features, y = features.diff)
        if (length(x = features) == 0) {
          #  stop("No features pass thresh.fc threshold")
          next
        }

        data1 <- data.use[features, cell.use1, drop = FALSE]
        data2 <- data.use[features, cell.use2, drop = FALSE]

        pvalues <- unlist(
          x = my.sapply(
            X = 1:nrow(x = data1),
            FUN = function(x) {
              # return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
              return(wilcox.test(data1[x, ], data2[x, ])$p.value)
            }
          )
        )

        # the genes removed by the pct/logFC pre-filtering are counted as tests (conservative adjustment)
        pval.adj = stats::p.adjust(
          p = pvalues,
          method = p.adjust.method,
          n = length(features.use)
        )
        genes.de[[i]] <- data.frame(clusters = level.use[i], features = as.character(rownames(data1)), pvalues = pvalues, logFC = FC[features], data.alpha[features,, drop = F],pvalues.adj = pval.adj, stringsAsFactors = FALSE)
      }

      markers.all <- data.frame()
      for (i in 1:numCluster) {
        gde <- genes.de[[i]]
        if (!is.null(gde)) {
          gde <- gde[order(gde$pvalues, -gde$logFC), ]
          gde <- subset(gde, subset = pvalues < thresh.p)
          if (!is.null(thresh.p.adj)) {
            gde <- subset(gde, subset = pvalues.adj < thresh.p.adj)
          }
          if (nrow(gde) > 0) {
            markers.all <- rbind(markers.all, gde)
          }
        }
      }
      if (only.pos & nrow(markers.all) > 0) {
        markers.all <- subset(markers.all, subset = logFC > 0)
      }
      if (!is.null(group.dataset) & nrow(markers.all) > 0) {
        markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
        markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), pos.dataset)
        markers.all$datasets <- factor(markers.all$datasets, levels = levels(labels.dataset))
        markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), ]
      }
      markers.all$features <- as.character(markers.all$features)

      features.sig <- markers.all$features
      object@var.features[[features.name]] <- features.sig
      features.name <- paste0(features.name, ".info")
      object@var.features[[features.name]] <- markers.all
    }

  } else {
    # select genes if they are exprssed in at least `min.cells` cells
    markers.all <- data.frame(features = as.character(rownames(data.use)), nCells = rowSums(data.use > 0))
    markers.all <- dplyr::filter(markers.all, nCells >= min.cells)
    features.sig <- markers.all$features
    object@var.features[[features.name]] <- unique(features.sig)
    features.name <- paste0(features.name, ".info")
    object@var.features[[features.name]] <- markers.all
  }


  if (do.DE && nrow(markers.all) == 0) {
    warning("No genes pass the thresholds of the differential expression analysis (thresh.p = ", thresh.p,
            ifelse(is.null(thresh.p.adj), "", paste0(", thresh.p.adj = ", thresh.p.adj)), ", thresh.fc = ", thresh.fc, ", thresh.pc = ", thresh.pc,
            ifelse(is.null(group.replicate), "", paste0(", min.cells.replicate = ", min.cells.replicate)), "). Please consider relaxing them.", call. = FALSE)
  }
  if (return.object) {
    return(object)
  } else {
    return(markers.all)
  }
}


#' Pseudobulk differential expression between conditions using biological replicates
#'
#' Internal function used by `identifyOverExpressedGenes` when `group.replicate` is set.
#'
#' @param data.use log-normalized expression matrix (genes x cells)
#' @param clusters cell group of each cell
#' @param level.use the cell groups to analyze
#' @param conditions a factor giving the condition of each cell; the first level is `pos.dataset` and the second level the other condition(s)
#' @param pos.dataset the condition of interest
#' @param replicates the replicate (sample) of each cell
#' @param method "limma" or "t.test"
#' @param paired whether the replicate labels identify the same donors in the different conditions
#' @param min.cells the minimum number of cells of a replicate in a cell group
#' @param covariates a data frame of covariates for each cell (constant within each replicate), or NULL
#' @param thresh.pc genes expressed in less than this fraction of cells in both conditions are not tested
#' @param min.pct.test genes expressed in less than this fraction of the cells of a cell group (conditions pooled) are not tested and not counted in the adjustment (independent filtering); NULL to count all the genes
#' @param p.adjust.method the method used to adjust the p-values within each cell group
#' @param n.features the number of features used for the adjustment (genes not tested are counted as tests)
#' @return a data frame with the columns clusters, features, pvalues, logFC, pct.1, pct.2, pvalues.adj, n.replicates.1, n.replicates.2
#' @keywords internal
pseudobulkDE <- function(data.use, clusters, level.use, conditions, pos.dataset, replicates, method = c("limma", "t.test"),
                         paired = FALSE, min.cells = 10, covariates = NULL, thresh.pc = 0, min.pct.test = NULL, p.adjust.method = "BH", n.features = nrow(data.use)) {
  method <- match.arg(method)
  if (!is.null(min.pct.test) && (!is.numeric(min.pct.test) || length(min.pct.test) != 1 || min.pct.test < 0 || min.pct.test >= 1)) {
    stop("`min.pct.test` should be a fraction of cells in [0, 1), e.g., 0.05")
  }
  if (!is.null(covariates) && method != "limma") {
    stop("`covariates` require replicate.method = 'limma'")
  }
  if (method == "limma" && !requireNamespace("limma", quietly = TRUE)) {
    message("The limma package is not installed: the Welch t-test is used instead. Install limma with `BiocManager::install('limma')`.")
    method <- "t.test"
  }
  is.pos <- as.character(conditions) == pos.dataset
  # the pseudobulk expression is computed from log-normalized data
  if (any(is.na(data.use))) {
    warning("The expression data contain NA values; they are set to 0 for the pseudobulk analysis.", call. = FALSE)
    data.use[is.na(data.use)] <- 0
  }
  if (any(data.use < 0)) {
    stop("The expression data contain negative values: `group.replicate` requires log-normalized data (e.g., the 'data' layer of Seurat, not the scaled data).")
  }
  if (max(data.use) > 50) {
    warning("The maximum expression value is ", round(max(data.use)), ": the data do not seem to be log-normalized, which is required by the pseudobulk analysis.", call. = FALSE)
  }
  # a sample is a replicate in a given condition (replicates may share labels between conditions)
  samples <- paste(ifelse(is.pos, "pos", "neg"), replicates, sep = "__")
  # covariates of each sample (they should be constant within a sample)
  if (!is.null(covariates)) {
    covariates[] <- lapply(covariates, function(x) if (is.numeric(x)) x else as.character(x))
    n.values <- sapply(covariates, function(x) tapply(x, samples, function(y) length(unique(y))))
    n.values <- matrix(n.values, ncol = ncol(covariates), dimnames = list(NULL, colnames(covariates)))
    if (any(n.values > 1)) {
      stop("The covariate(s) ", toString(colnames(covariates)[colSums(n.values > 1) > 0]), " are not constant within the replicates; they cannot be used as sample-level covariates.")
    }
    covariates.sample <- covariates[!duplicated(samples), , drop = FALSE]
    rownames(covariates.sample) <- samples[!duplicated(samples)]
  }
  reps.per.condition <- tapply(as.character(replicates), ifelse(is.pos, pos.dataset, "the other condition(s)"), function(x) sort(unique(x)))
  if (any(lengths(reps.per.condition) < 2)) {
    stop("At least two biological replicates are required in each condition, but `group.replicate` gives: ",
         paste0(names(reps.per.condition), ": ", sapply(reps.per.condition, toString), collapse = "; "), ". ",
         "NB: when the CellChat object is created without sample information, the column `samples` contains 'sample1' for all cells. ",
         "Please use the metadata column giving the samples (e.g., group.replicate = 'orig.ident'), check it with `table(object@meta[[group.dataset]], object@meta[[group.replicate]])`, ",
         "or remove `group.replicate` to perform the cell-level test.", call. = FALSE)
  }
  if (!paired) {
    reps.pos <- unique(replicates[is.pos]); reps.neg <- unique(replicates[!is.pos])
    if (length(intersect(reps.pos, reps.neg)) > 0) {
      message("Some replicate labels are found in both conditions (", toString(utils::head(intersect(reps.pos, reps.neg), 3)),
              "); they are treated as different samples. Set `paired.replicate = TRUE` if they are the same donors.")
    }
  }
  res <- list()
  skipped <- c()
  for (cl in level.use) {
    idx <- which(clusters == cl)
    if (length(idx) == 0) next
    n.cells <- table(samples[idx])
    samples.use <- names(n.cells)[n.cells >= min.cells]
    sample.pos <- grepl("^pos__", samples.use)
    if (sum(sample.pos) < 2 || sum(!sample.pos) < 2) {
      skipped <- c(skipped, paste0(cl, " (", sum(sample.pos), " vs. ", sum(!sample.pos), ")"))
      next
    }
    idx <- idx[samples[idx] %in% samples.use]
    # pseudobulk expression: log of the average normalized expression of the cells of each replicate
    pb <- sapply(samples.use, function(x) {
      log1p(Matrix::rowMeans(expm1(data.use[, idx[samples[idx] == x], drop = FALSE])))
    })
    pb <- matrix(pb, nrow = nrow(data.use), dimnames = list(rownames(data.use), samples.use))
    donors <- sub("^(pos|neg)__", "", samples.use)
    if (paired) {
      donors.both <- intersect(donors[sample.pos], donors[!sample.pos])
      if (length(donors.both) < 2) {
        skipped <- c(skipped, paste0(cl, " (", length(donors.both), " paired donor(s))"))
        next
      }
      keep <- donors %in% donors.both
      pb <- pb[, keep, drop = FALSE]; sample.pos <- sample.pos[keep]; donors <- donors[keep]; samples.use <- samples.use[keep]
      idx <- idx[samples[idx] %in% samples.use]
    }
    # fraction of expressing cells in each condition
    pct.1 <- Matrix::rowMeans(data.use[, idx[is.pos[idx]], drop = FALSE] > 0)
    pct.2 <- Matrix::rowMeans(data.use[, idx[!is.pos[idx]], drop = FALSE] > 0)
    testable <- pmax(pct.1, pct.2) > thresh.pc & apply(pb, 1, stats::var) > 0 & rowSums(!is.finite(pb)) == 0
    n.tests <- max(n.features, sum(testable))
    if (!is.null(min.pct.test)) {
      # independent filtering: the fraction of expressing cells of the cell group, all the compared cells pooled, does not depend on the condition labels;
      # the genes below `min.pct.test` are neither tested nor counted in the adjustment
      pct.pooled <- Matrix::rowMeans(data.use[, idx, drop = FALSE] > 0)
      testable <- testable & pct.pooled >= min.pct.test
      n.tests <- sum(pct.pooled >= min.pct.test)
    }
    features <- rownames(pb)[which(testable)]
    if (length(features) == 0) next
    pb <- pb[features, , drop = FALSE]
    if (paired) {
      pos.mat <- pb[, sample.pos, drop = FALSE][, order(donors[sample.pos]), drop = FALSE]
      neg.mat <- pb[, !sample.pos, drop = FALSE][, order(donors[!sample.pos]), drop = FALSE]
      logFC <- rowMeans(pos.mat - neg.mat)
    } else {
      logFC <- rowMeans(pb[, sample.pos, drop = FALSE]) - rowMeans(pb[, !sample.pos, drop = FALSE])
    }
    if (method == "limma") {
      design.data <- data.frame(condition = factor(ifelse(sample.pos, "pos", "neg"), levels = c("neg", "pos")), row.names = samples.use)
      terms <- c()
      if (paired) {
        design.data$donor <- factor(donors)
        terms <- c(terms, "donor")
      }
      if (!is.null(covariates)) {
        cov.use <- covariates.sample[samples.use, , drop = FALSE]
        # covariates without variation in this cell group, or nested within the donors, are not informative
        keep.cov <- sapply(cov.use, function(x) length(unique(x)) > 1)
        if (paired) {
          keep.cov <- keep.cov & sapply(cov.use, function(x) any(tapply(x, donors, function(y) length(unique(y))) > 1))
        }
        cov.use <- cov.use[, keep.cov, drop = FALSE]
        if (ncol(cov.use) > 0) {
          design.data <- cbind(design.data, cov.use)
          terms <- c(terms, colnames(cov.use))
        }
      }
      terms <- if (length(terms) > 0) paste0("`", terms, "`") else character(0)
      design <- stats::model.matrix(stats::as.formula(paste("~", paste(c(terms, "condition"), collapse = " + "))), data = design.data)
      if (qr(design)$rank < ncol(design)) {
        stop("Cell group ", cl, ": the condition is confounded with the covariates (e.g., all the samples of one condition have the same sex); the effect of the condition cannot be estimated.")
      }
      if (nrow(design) - ncol(design) < 1) {
        skipped <- c(skipped, paste0(cl, " (no residual degree of freedom)"))
        next
      }
      fit <- limma::eBayes(limma::lmFit(pb, design), trend = TRUE)
      pvalues <- fit$p.value[, "conditionpos"]
      # effect of the condition adjusted for the covariates/donors
      logFC <- fit$coefficients[, "conditionpos"]
    } else {
      pvalues <- sapply(features, function(g) {
        if (paired) {
          d <- pos.mat[g, ] - neg.mat[g, ]
          if (stats::sd(d) == 0) return(ifelse(mean(d) == 0, 1, 0))
          stats::t.test(pos.mat[g, ], neg.mat[g, ], paired = TRUE)$p.value
        } else {
          x <- pb[g, sample.pos]; y <- pb[g, !sample.pos]
          if (stats::sd(x) == 0 && stats::sd(y) == 0) return(ifelse(mean(x) == mean(y), 1, 0))
          stats::t.test(x, y)$p.value
        }
      })
    }
    pvalues[is.na(pvalues)] <- 1
    res[[cl]] <- data.frame(clusters = cl, features = features, pvalues = as.numeric(pvalues), logFC = as.numeric(logFC),
                            pct.1 = as.numeric(pct.1[features]), pct.2 = as.numeric(pct.2[features]),
                            pvalues.adj = stats::p.adjust(as.numeric(pvalues), method = p.adjust.method, n = max(n.tests, length(features))),
                            n.replicates.1 = sum(sample.pos), n.replicates.2 = sum(!sample.pos), stringsAsFactors = FALSE)
  }
  if (length(skipped) > 0) {
    message(length(skipped), " of ", length(level.use), " cell group(s) skipped (less than two replicates with at least ", min.cells,
            " cells in each condition; number of usable replicates in brackets): ", toString(skipped))
  }
  if (length(res) == 0) {
    stop("No cell group has at least two replicates with at least `min.cells.replicate` = ", min.cells, " cells in each condition: ",
         toString(skipped), ". Please decrease `min.cells.replicate` or check the replicates with `table(object@meta[[group.replicate]], object@idents$joint)`.", call. = FALSE)
  }
  do.call(rbind, res)
}


#' Identify over-expressed ligands and (complex) receptors associated with each cell group
#'
#' This function identifies the over-expressed ligands and (complex) receptors based on the identified signaling genes from 'identifyOverExpressedGenes'.
#'
#' @param object CellChat object
#' @param features.name a char name used for storing the over-expressed ligands and receptors in `object@var.features[[paste0(features.name, ".LR")]]`
#' @param features a vector of features to use. default use all over-expressed genes in `object@var.features[[features.name]]`
#' @param return.object whether returning a CellChat object. If FALSE, it will return a data frame containing over-expressed ligands and (complex) receptors associated with each cell group
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom dplyr select
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, a new element named paste0(features.name, ".LR") will be added into the list `object@var.features`
#' @export
#'
identifyOverExpressedLigandReceptor <- function(object, features.name = "features", features = NULL, return.object = TRUE) {

  features.name.LR <- paste0(features.name, ".LR")
  features.name <- paste0(features.name, ".info")
  DB <- object@DB
  interaction_input <- DB$interaction
  complex_input <- DB$complex
  pairLR <- select(interaction_input, ligand, receptor)
  LR.use <- unique(c(pairLR$ligand, pairLR$receptor))
  if (is.null(features)) {
    if (is.list(object@var.features)) {
      markers.all <- object@var.features[[features.name]] # use the updated CellChat object 12/2020
    } else {
      stop("Please update your CellChat object via `updateCellChat()`")
    }

  } else {
    features.use <- features
    rm(features)
    markers.all <- subset(markers.all, subset = features %in% features.use)
  }

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  complexSubunits <- complex_input[, grepl("subunit" , colnames(complex_input))]

  markers.all.new <- data.frame()
  for (i in 1:nrow(markers.all)) {
    if (markers.all$features[i] %in% LR.use) {
      markers.all.new <- rbind(markers.all.new, markers.all[i, , drop = FALSE])
    } else {
      index.sig <- unlist(
        x = my.sapply(
          X = 1:nrow(complexSubunits),
          FUN = function(x) {
            complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
            complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
            if (markers.all$features[i] %in% complexsubunitsV) {
              return(x)
            }
          }
        )
      )
      complexSubunits.sig <- rownames(complexSubunits[index.sig,])
      markers.all.complex <- data.frame()
      for (j in 1:length(complexSubunits.sig)) {
        markers.all.complex <- rbind(markers.all.complex, markers.all[i, , drop = FALSE])
      }
      markers.all.complex$features <- complexSubunits.sig
      markers.all.new <- rbind(markers.all.new, markers.all.complex)
    }
  }

  object@var.features[[features.name.LR]] <- markers.all.new

  if (return.object) {
    return(object)
  } else {
    return(markers.all.new)
  }
}



#' Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#'
#' @param object CellChat object
#' @param features.name a char name used for assess the results in `object@var.features[[features.name]]`
#' @param features a vector of features to use. default use all over-expressed genes in `object@var.features[[features.name]]`
#' @param variable.both variable.both = TRUE will require that both ligand and receptor from one pair are over-expressed;
#'
#' variable.both = FALSE will only require that either ligand or receptor from one pair is over-expressed, leading to more over-expressed ligand-receptor interactions (pairs) for further analysis.
#' @param return.object whether returning a CellChat object. If FALSE, it will return a data frame containing the over-expressed ligand-receptor pairs
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom dplyr select
#'
#' @return A CellChat object or a data frame. If returning a CellChat object, a new element named 'LRsig' will be added into the list `object@LR`
#' @export
#'
identifyOverExpressedInteractions <- function(object, features.name = "features", variable.both = TRUE, features = NULL, return.object = TRUE) {
  gene.use <- row.names(object@data.signaling)
  DB <- object@DB
  if (is.null(features)) {
    if (is.list(object@var.features)) {
      features.sig <- object@var.features[[features.name]] # use the updated CellChat object 12/2020
    } else {
      stop("Please update your CellChat object via `updateCellChat()`")
    }

  } else {
    features.sig <- features
  }

  interaction_input <- DB$interaction
  complex_input <- DB$complex
  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )
  complexSubunits <- complex_input[, grepl("subunit" , colnames(complex_input))]
  index.sig <- unlist(
    x = my.sapply(
      X = 1:nrow(complexSubunits),
      FUN = function(x) {
        complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
        complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
        if (length(intersect(complexsubunitsV, features.sig)) > 0 & all(complexsubunitsV %in% gene.use)) {
          return(x)
        }
      }
    )
  )
  complexSubunits.sig <- complexSubunits[index.sig,]

  index.use <- unlist(
    x = my.sapply(
      X = 1:nrow(complexSubunits),
      FUN = function(x) {
        complexsubunitsV <- unlist(complexSubunits[x,], use.names = F)
        complexsubunitsV <- complexsubunitsV[complexsubunitsV != ""]
        if (all(complexsubunitsV %in% gene.use)) {
          return(x)
        }
      }
    )
  )
  complexSubunits.use <- complexSubunits[index.use,]

  pairLR <- select(interaction_input, ligand, receptor)

  if (variable.both) {
    index.sig <- unlist(
      x = my.sapply(
        X = 1:nrow(pairLR),
        FUN = function(x) {
          if (all(unlist(pairLR[x,], use.names = F) %in% c(features.sig, rownames(complexSubunits.sig)))) {
            return(x)
          }
        }
      )
    )
  } else {
    index.sig <- unlist(
      x = my.sapply(
        X = 1:nrow(pairLR),
        FUN = function(x) {
          # if (all(unlist(pairLR[x,], use.names = F) %in% c(features.sig, rownames(complexSubunits.sig)))) {
          if (all(unlist(pairLR[x,], use.names = F) %in% c(gene.use, rownames(complexSubunits.use))) & (length(intersect(unlist(pairLR[x,], use.names = F), c(features.sig, rownames(complexSubunits.sig)))) > 0)) {
            return(x)
          }
        }
      )
    )
  }

  pairLRsig <- interaction_input[index.sig, ]
  object@LR$LRsig <- pairLRsig
  cat("The number of highly variable ligand-receptor pairs used for signaling inference is", nrow(pairLRsig), '\n')
  if (nrow(pairLRsig) == 0) {
    warning("No ligand-receptor pairs are selected: ", length(features.sig), " over-expressed signaling genes are available in `object@var.features[['", features.name, "']]`. ",
            "Please relax the thresholds of `identifyOverExpressedGenes` (the default thresholds are recommended before `computeCommunProb`).", call. = FALSE)
  }
  if (return.object) {
    return(object)
  } else {
    return(pairLRsig)
  }
}


#' Smooth the gene expression data
#'
#' A diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network.
#'
#' This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors
#'
#' @param object  CellChat object
#' @param method When method = "netSmooth", smoothing a gene’s expression values based on its neighbors defined in a high-confidence experimentally validated protein-protein network.
#' @param adj adjacency matrix of protein-protein interaction network to use
#' @param alpha numeric in [0,1] alpha = 0: no smoothing; a larger value alpha results in increasing levels of smoothing.
#' @param normalizeAdjMatrix    how to normalize the adjacency matrix
#'                              possible values are 'rows' (in-degree)
#'                              and 'columns' (out-degree)
#' @return a smoothed gene expression matrix
#' @export
#'
# This function is adapted from https://github.com/BIMSBbioinfo/netSmooth
smoothData <- function(object, method = c("netSmooth"), adj = NULL, alpha=0.5, normalizeAdjMatrix=c('rows','columns')){
  if ("data.smooth" %in% methods::slotNames(object) == FALSE) {
    stop("`object@data.smooth` is missing. Please update the CellChat object via `updateCellChat`! \n")
  }
  data <- as.matrix(object@data.signaling)
  normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
  if (method == "netSmooth") {
    if (is.null(adj)) stop("Please provide the `adj`. \n")
    stopifnot(is(adj, 'matrix') | is(adj, 'sparseMatrix'))
    stopifnot((is.numeric(alpha) & (alpha > 0 & alpha < 1)))
    if(sum(Matrix::rowSums(adj)==0)>0) stop("PPI cannot have zero rows/columns")
    if(sum(Matrix::colSums(adj)==0)>0) stop("PPI cannot have zero rows/columns")
  }
  if(is.numeric(alpha)) {
    if(alpha<0 | alpha > 1) {
      stop('alpha must be between 0 and 1')
    }
    data.projected <- projectAndRecombine(data, adj, alpha,normalizeAdjMatrix=normalizeAdjMatrix)
  } else stop("unsupported alpha value: ", class(alpha))
  object@data.smooth <- data.projected
  return(object)
}

#' Perform network projecting on network when the network genes and the
#' experiment genes aren't exactly the same.
#'
#' The gene network might be defined only on a subset of genes that are
#' measured in any experiment. Further, an experiment might not measure all
#' genes that are present in the network. This function projects the experiment
#' data onto the gene space defined by the network prior to projecting. Then,
#' it projects the projected data back into the original dimansions.
#'
#' @param gene_expression  gene expession data to be projected
#'                         [N_genes x M_samples]
#' @param adj_matrix  adjacenty matrix of network to perform projecting over.
#'                    Will be column-normalized.
#'                    Rownames and colnames should be genes.
#' @param alpha  network projecting parameter (1 - restart probability in random
#'                walk model.
#' @param projecting.function  must be a function that takes in data, adjacency
#'                            matrix, and alpha. Will be used to perform the
#'                            actual projecting.
#' @param normalizeAdjMatrix    which dimension (rows or columns) should the
#'                              adjacency matrix be normalized by. rows
#'                              corresponds to in-degree, columns to
#'                              out-degree.
#' @return  matrix with network-projected gene expression data. Genes that are
#'          not present in projecting network will retain original values.
#' @keywords internal
#'
projectAndRecombine <- function(gene_expression, adj_matrix, alpha,
                                projecting.function=randomWalkBySolve,
                                normalizeAdjMatrix=c('rows','columns')) {
  normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
  gene_expression_in_A_space <- projectOnNetwork(gene_expression,rownames(adj_matrix))
  gene_expression_in_A_space_project <- projecting.function(gene_expression_in_A_space, adj_matrix, alpha, normalizeAdjMatrix)
  gene_expression_project <- projectFromNetworkRecombine(gene_expression, gene_expression_in_A_space_project)
  return(gene_expression_project)
}


#' Project the gene expression matrix onto a lower space
#' of the genes defined in the projecting network
#' @param gene_expression    gene expression matrix
#' @param new_features       the genes in the network, on which to project
#'                           the gene expression matrix
#' @param missing.value      value to assign to genes that are in network,
#'                           but missing from gene expression matrix
#' @return the gene expression matrix projected onto the gene space defined by new_features
#' @keywords internal
projectOnNetwork <- function(gene_expression, new_features, missing.value=0) {
 # data_in_new_space = matrix(rep(0, length(new_features)*dim(gene_expression)[2]),nrow=length(new_features))
  data_in_new_space = matrix(0, ncol=dim(gene_expression)[2], nrow=length(new_features))
  rownames(data_in_new_space) <- new_features
  colnames(data_in_new_space) <- colnames(gene_expression)
  genes_in_both <- intersect(rownames(data_in_new_space),rownames(gene_expression))
  data_in_new_space[genes_in_both,] <- gene_expression[genes_in_both,]
  genes_only_in_network <- setdiff(new_features, rownames(gene_expression))
  data_in_new_space[genes_only_in_network,] <- missing.value
  return(data_in_new_space)
}

#' project data on graph by solving the linear equation (I - alpha*A) * E_sm = E * (1-alpha)

#' @param E      initial data matrix [NxM]
#' @param A      adjacency matrix of graph to network project on will be column-normalized.
#' @param alpha  projecting coefficient (1 - restart probability of random walk)
#' @return network-projected gene expression
#' @keywords internal
randomWalkBySolve <- function(E, A, alpha, normalizeAjdMatrix=c('rows','columns')) {
  normalizeAjdMatrix <- match.arg(normalizeAjdMatrix)
  if (normalizeAjdMatrix=='rows') {
    Anorm <- l1NormalizeRows(A)
  } else if (normalizeAjdMatrix=='columns') {
    Anorm <- l1NormalizeColumns(A)
  }
  eye <- diag(dim(A)[1])
  AA <- eye - alpha*Anorm
  BB <- (1-alpha) * E
  return(solve(AA, BB))
}

#' Column-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' column sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeColumns(A)
#' @return column-normalized sparse matrix object
#' @keywords internal
l1NormalizeColumns <- function(A) {
  return(Matrix::t(Matrix::t(A)/Matrix::colSums(A)))
}

#' Row-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' row sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeRows(A)
#' @return row-normalized sparse matrix object
#' @keywords internal
l1NormalizeRows <- function(A) {
  return(A/Matrix::rowSums(A))
}

#' Combine gene expression from projected space (that of the network) with the
#' expression of genes that were not projected (not present in network)
#' @keywords internal
#' @param original_expression    the non-projected expression
#' @param projected_expression    the projected gene expression, in the space
#'                               of the genes defined by the network
#' @return a matrix in the dimensions of original_expression, where values that
#'         are present in projected_expression are copied from there.
projectFromNetworkRecombine <- function(original_expression, projected_expression) {
  data_in_original_space <- original_expression
  genes_in_both <- intersect(rownames(original_expression),rownames(projected_expression))
  data_in_original_space[genes_in_both,] <- as.matrix(projected_expression[genes_in_both,])
  return(data_in_original_space)
}


#' Dimension reduction using PCA
#'
#' @param data.use input data (samples in rows, features in columns)
#' @param do.fast whether do fast PCA
#' @param dimPC the number of components to keep
#' @param seed.use set a seed
#' @param weight.by.var whether use weighted pc.scores
#' @importFrom stats prcomp
#' @importFrom irlba irlba
#' @return
#' @export
#'
#' @examples
runPCA <- function(data.use, do.fast = T, dimPC = 50, seed.use = 42, weight.by.var = T) {
  set.seed(seed = seed.use)
  if (do.fast) {
    dimPC <- min(dimPC, ncol(data.use) - 1)
    pca.res <- irlba::irlba(data.use, nv = dimPC)
    sdev <- pca.res$d/sqrt(max(1, nrow(data.use) - 1))
    if (weight.by.var){
      pc.scores <- pca.res$u %*% diag(pca.res$d)
    } else {
      pc.scores <- pca.res$u
    }
  } else {
    dimPC <- min(dimPC, ncol(data.use) - 1)
    pca.res <- stats::prcomp(x = data.use, rank. = dimPC)
    sdev <- pca.res$sdev
    if (weight.by.var) {
      pc.scores <- pca.res$x %*% diag(pca.res$sdev[1:dimPC]^2)
    } else {
      pc.scores <- pca.res$x
    }
  }
  rownames(pc.scores) <- rownames(data.use)
  colnames(pc.scores) <- paste0('PC', 1:ncol(pc.scores))
  return(pc.scores)
}


#' Run UMAP
#' @param data.use input data matrix
#' @param n_neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n_components The dimension of the space to embed into.
#' @param metric This determines the choice of metric used to measure distance in the input space.
#' @param n_epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning_rate The initial learning rate for the embedding optimization.
#' @param min_dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set_op_mix_ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local_connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion_strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative_sample_rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @param metric_kwds,angular_rp_forest,verbose other parameters used in UMAP
#' @import reticulate
#' @export
#'
runUMAP <- function(
  data.use,
  n_neighbors = 30L,
  n_components = 2L,
  metric = "correlation",
  n_epochs = NULL,
  learning_rate = 1.0,
  min_dist = 0.3,
  spread = 1.0,
  set_op_mix_ratio = 1.0,
  local_connectivity = 1L,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric_kwds = NULL,
  angular_rp_forest = FALSE,
  verbose = FALSE){
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).")
  }
  set.seed(seed.use)
  reticulate::py_set_seed(seed.use)
  umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(n_neighbors),
    n_components = as.integer(n_components),
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    min_dist = min_dist,
    spread = spread,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    metric_kwds = metric_kwds,
    angular_rp_forest = angular_rp_forest,
    verbose = verbose
  )
  Rumap <- umap$fit_transform
  umap_output <- Rumap(t(data.use))
  colnames(umap_output) <- paste0('UMAP', 1:ncol(umap_output))
  rownames(umap_output) <- colnames(data.use)
  return(umap_output)
}

.error_if_no_Seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat installation required for working with Seurat objects")
  }
}


#' Color interpolation
#'
#' This function is modified from https://rdrr.io/cran/circlize/src/R/utils.R
#' Colors are linearly interpolated according to break values and corresponding colors through CIE Lab color space (`colorspace::LAB`) by default.
#' Values exceeding breaks will be assigned with corresponding maximum or minimum colors.
#'
#' @param breaks A vector indicating numeric breaks
#' @param colors A vector of colors which correspond to values in ``breaks``
#' @param transparency A single value in ``[0, 1]``. 0 refers to no transparency and 1 refers to full transparency
#' @param space color space in which colors are interpolated. Value should be one of "RGB", "HSV", "HLS", "LAB", "XYZ", "sRGB", "LUV", see `colorspace::color-class` for detail.
#' @importFrom colorspace coords RGB HSV HLS LAB XYZ sRGB LUV hex
#' @importFrom grDevices col2rgb
#' @return It returns a function which accepts a vector of numeric values and returns interpolated colors.
#' @export
#' @examples
#' \dontrun{
#' col_fun = colorRamp3(c(-1, 0, 1), c("green", "white", "red"))
#' col_fun(c(-2, -1, -0.5, 0, 0.5, 1, 2))
#' }
colorRamp3 = function(breaks, colors, transparency = 0, space = "LAB") {

  if(length(breaks) != length(colors)) {
    stop("Length of `breaks` should be equal to `colors`.\n")
  }

  colors = colors[order(breaks)]
  breaks = sort(breaks)

  l = duplicated(breaks)
  breaks = breaks[!l]
  colors = colors[!l]

  if(length(breaks) == 1) {
    stop("You should have at least two distinct break values.")
  }


  if(! space %in% c("RGB", "HSV", "HLS", "LAB", "XYZ", "sRGB", "LUV")) {
    stop("`space` should be in 'RGB', 'HSV', 'HLS', 'LAB', 'XYZ', 'sRGB', 'LUV'")
  }

  colors = t(grDevices::col2rgb(colors)/255)

  attr = list(breaks = breaks, colors = colors, transparency = transparency, space = space)

  if(space == "LUV") {
    i = which(apply(colors, 1, function(x) all(x == 0)))
    colors[i, ] = 1e-5
  }

  transparency = 1-ifelse(transparency > 1, 1, ifelse(transparency < 0, 0, transparency))[1]
  transparency_str = sprintf("%X", round(transparency*255))
  if(nchar(transparency_str) == 1) transparency_str = paste0("0", transparency_str)

  fun = function(x = NULL, return_rgb = FALSE, max_value = 1) {
    if(is.null(x)) {
      stop("Please specify `x`\n")
    }

    att = attributes(x)
    if(is.data.frame(x)) x = as.matrix(x)

    l_na = is.na(x)
    if(all(l_na)) {
      return(rep(NA, length(l_na)))
    }

    x2 = x[!l_na]

    x2 = ifelse(x2 < breaks[1], breaks[1],
                ifelse(x2 > breaks[length(breaks)], breaks[length(breaks)],
                       x2
                ))
    ibin = .bincode(x2, breaks, right = TRUE, include.lowest = TRUE)
    res_col = character(length(x2))
    for(i in unique(ibin)) {
      l = ibin == i
      res_col[l] = .get_color(x2[l], breaks[i], breaks[i+1], colors[i, ], colors[i+1, ], space = space)
    }
    res_col = paste(res_col, transparency_str[1], sep = "")

    if(return_rgb) {
      res_col = t(grDevices::col2rgb(as.vector(res_col), alpha = TRUE)/255)
      return(res_col)
    } else {
      res_col2 = character(length(x))
      res_col2[l_na] = NA
      res_col2[!l_na] = res_col

      attributes(res_col2) = att
      return(res_col2)
    }
  }

  attributes(fun) = attr
  return(fun)
}

.restrict_in = function(x, lower, upper) {
  x[x > upper] = upper
  x[x < lower] = lower
  x
}

# x: vector
# break1 single value
# break2 single value
# rgb1 vector with 3 elements
# rgb2 vector with 3 elements
.get_color = function(x, break1, break2, col1, col2, space) {

  col1 = colorspace::coords(as(colorspace::sRGB(col1[1], col1[2], col1[3]), space))
  col2 = colorspace::coords(as(colorspace::sRGB(col2[1], col2[2], col2[3]), space))

  res_col = matrix(ncol = 3, nrow = length(x))
  for(j in 1:3) {
    xx = (x - break2)*(col2[j] - col1[j]) / (break2 - break1) + col2[j]
    res_col[, j] = xx
  }

  res_col = get(space)(res_col)
  res_col = colorspace::coords(as(res_col, "sRGB"))
  res_col[, 1] = .restrict_in(res_col[,1], 0, 1)
  res_col[, 2] = .restrict_in(res_col[,2], 0, 1)
  res_col[, 3] = .restrict_in(res_col[,3], 0, 1)
  colorspace::hex(colorspace::sRGB(res_col))
}

#' Update the cell-cell communication array from a customized cell-cell-communication scores between different cell groups
#'
#' Users may also check the `updateCellChatDB` function for integrating other resources or utilizing a custom database
#'
#' @param object CellChat object
#' @param net a data frame with at least five columns named as `source`,`target`,`ligand`,`receptor` and `score`, which defines the customized cell-cell-communication scores between different cell groups.
#' a p-value column named `pval`, and additional columns named `interaction_name` and `interaction_name_2` can be also provided.
#' @return a CellChat object with updated slot `net` and slot `DB` if db is not NULL.
#' @export

updateCCC_score <- function(object, net) {
  df.net <- net
  if (all(c("source","target","ligand","receptor","score") %in% colnames(df.net)) == FALSE) {
    stop("The input `net` must contain at least five columns named as source,target,ligand,receptor,score")
  }
  if (all(c("interaction_name") %in% colnames(df.net)) == FALSE) {
    df.net$interaction_name <- paste0(toupper(df.net$ligand), "_", toupper(df.net$receptor))
  }
  if (all(c("interaction_name_2") %in% colnames(df.net)) == FALSE) {
    df.net$interaction_name_2 <- paste0(df.net$ligand, " - ", df.net$receptor)
  }
  if (all(c("pval") %in% colnames(df.net)) == FALSE) {
    df.net$pval <- rep(0, nrow(df.net))
  }
  df.net$prob <- df.net$score

  LR <- unique(df.net$interaction_name)
  cell.levels <- levels(object@idents)
  numCluster <- length(cell.levels)
  mat.prob.all <- array(0, dim = c(numCluster,numCluster,length(LR)))
  mat.pval.all <- mat.prob.all
  for (i in 1:length(LR)) {
    df.i <- df.net[df.net$interaction_name == LR[i], , drop = FALSE]
    mat.prob <- matrix(0, nrow = numCluster, ncol = numCluster)
    mat.pval <- mat.prob
    for (j in 1:nrow(df.i)) {
      idx.s <- which(df.i$source[j] == cell.levels)
      idx.t <- which(df.i$target[j] == cell.levels)
      mat.prob[idx.s, idx.t] <- df.i$prob[j]
      mat.pval[idx.s, idx.t] <- df.i$pval[j]
    }
    mat.prob.all[,,i] <- mat.prob
    mat.pval.all[,,i] <- mat.pval
  }

  dimnames(mat.prob.all) <- list(cell.levels, cell.levels, LR)
  dimnames(mat.pval.all) <- dimnames(mat.prob.all)
  net <- list("prob" = mat.prob.all, "pval" = mat.pval.all)
  object@net <- net

  return(object)
}

#' Preprocessing multi-omics data and preparing the L-R database
#'
#' @param data.list a list consisting of multi-omics data (e.g., RNA & ADT)
#' @param db one of the CellChatDB databases: CellChatDB.human, CellChatDB.mouse, CellChatDB.zebrafish
#' @param cutoff the cutoff value of low protein expression
#' @param do.sparse whether to use sparse format
#' @export
#'
preProcMultiomics <- function(data.list, db, cutoff = 0.5, do.sparse = TRUE) {
  # normalize the data
  data.input.rna <- data.list[[1]]
  data.input.adt <- data.list[[2]]
  data.input.rna = data.input.rna/max(data.input.rna)
  data.input.adt = data.input.adt/max(data.input.adt)
  data.input.adt.temp = data.input.adt
  X = data.input.adt
  for (i in 1:nrow(X)) {
    data.input.adt.temp[i,] = (X[i,]-min(X[i,]))/(max(X[i,])-min(X[i,]))
  }
  data.input.adt[data.input.adt.temp < cutoff] <- 0
  if (do.sparse) {
    data.input = rbind(data.input.rna, as(data.input.adt, "dgCMatrix"))
  } else {
    data.input = rbind(as.matrix(data.input.rna), as.matrix(data.input.adt))
  }
  
  # create a new L-R database
  proteins <- rownames(data.input.adt)
  geneInfo.subset <- db$geneInfo[db$geneInfo$AntibodyName %in% proteins, ]
  proteins.nonmapping <- setdiff(proteins, geneInfo.subset$AntibodyName)
  if (length(proteins.nonmapping) > 0) {
    warning(cat("The following antibodies are not found in `CellChatDB$geneInfo$AntibodyName`: ", toString(proteins.nonmapping), "! Please manually add them via the function `updateCellChatDB`. \n"))
  }
  out <- extractLRfromGenes(geneSet = geneInfo.subset$Symbol, db)
  LR.use <- out$LR.use
  idx <- match(LR.use$ligand, geneInfo.subset$Symbol)
  LR.use$ligand[!is.na(idx)] <- geneInfo.subset$AntibodyName[idx[!is.na(idx)]]
  idx <- match(LR.use$receptor, geneInfo.subset$Symbol)
  LR.use$receptor[!is.na(idx)] <- geneInfo.subset$AntibodyName[idx[!is.na(idx)]]
  
  db.use <- db
  db.use$interaction <- LR.use
  db.use$geneInfo <- dplyr::add_row(db.use$geneInfo, Symbol = geneInfo.subset$AntibodyName)
  
  return(list(data.input = data.input, db.use = db.use))
}


                    


#' Check that exactly two datasets are given for a pairwise comparison
#'
#' Functions showing differences between two datasets (e.g., `netVisual_diffInteraction`) only use the first two values of `comparison`.
#' With three or more conditions, these functions should be run for each pair of datasets, e.g., comparison = c(1,2), c(1,3) and c(2,3).
#'
#' @param comparison a numerical vector giving the datasets for comparison
#' @param n.datasets the number of datasets in the merged object
#' @return invisible NULL; raises an error if `comparison` is not valid
#' @keywords internal
checkPairwiseComparison <- function(comparison, n.datasets) {
  if (length(comparison) != 2) {
    pairs <- utils::combn(seq_len(max(n.datasets, 2)), 2)
    pairs <- paste0("c(", pairs[1, ], ",", pairs[2, ], ")", collapse = ", ")
    stop("This function compares exactly two datasets, but `comparison` contains ", length(comparison), " values. ",
         "For more than two conditions, please run it for each pair of datasets, e.g., comparison = ", pairs, ". ",
         "Use `rankNet`, `netVisual_bubble` or `compareInteractions` to compare all the conditions at once.", call. = FALSE)
  }
  if (!is.numeric(comparison) || any(comparison < 1) || any(comparison > n.datasets) || comparison[1] == comparison[2]) {
    stop("`comparison` should contain the indices of two different datasets among 1:", n.datasets, call. = FALSE)
  }
  invisible(NULL)
}
