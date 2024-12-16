
#' Compute the communication probability/strength between any interacting cell groups
#'
#' To further speed up on large-scale datasets, USER can downsample the data using the function 'subset' from Seurat package (e.g., pbmc.small <- subset(pbmc, downsample = 500)), or using the function `sketchData` from CellChat, in particular for the large cell clusters;
#'
#'
#' @param object CellChat object
#' @param type Methods for computing the average gene expression per cell group. By default = "triMean", producing fewer but stronger interactions;
#' When setting `type = "truncatedMean"`, a value should be assigned to 'trim',  producing more interactions.
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed
#' @param LR.use A subset of ligand-receptor interactions used in inferring communication network
#' @param raw.use Whether use the raw data (i.e., `object@data.signaling`) or the smoothed data (i.e., `object@data.smooth`).
#' Set raw.use = FALSE to use the projected data when analyzing single-cell data with shallow sequencing depth because the projected data could help to reduce the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors.
#' @param population.size Whether consider the proportion of cells in each group across all sequenced cells.
#' Set population.size = FALSE if analyzing sorting-enriched single cells, to remove the potential artifact of population size.
#' Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, with the reason that abundant cell populations tend to send collectively stronger signals than the rare cell populations.
#'
#' Parameters for spatial data analysis:
#' @param distance.use Whether to use distance constraints to compute communication probability. Setting `distance.use = TRUE` indicates that the cell-cell communication probability is inversely proportional to the computed distance.
#' Setting `distance.use = FALSE` will only filter out interactions between spatially distant regions, but not add distance constraints.
#' @param interaction.range The maximum interaction/diffusion length of ligands (Unit: microns). This hard threshold is used to filter out the connections between spatially distant regions
#' @param scale.distance A scale or normalization factor for the spatial distances when setting `distance.use = TRUE`. For example, scale.distance equals 1, 0.1, 0.01, 0.001, 0.11, or 0.011. We choose this values such that the minimum value of the scaled distances is in [1,2]. This value is not necessary when setting `distance.use = FALSE`.
#'
#' When comparing communication across different CellChat objects, the same scale factor should be used. For a single CellChat analysis, different scale factors will not affect the ranking of the signaling based on their interaction strength.
#'
#' @param k.min The minimum number of interacting cell pairs required for defining spatially proximal cell groups.
#' @param contact.dependent Whether using the `contact-dependent` manner for inference signaling, that is determining interacting cell pairs by requiring cells to be in direct membrane-membrane contact. By default `contact.dependent = TRUE` when inferring contact-dependent and juxtacrine signaling (that is "Cell-Cell Contact" signaling classified in CellChatDB$interaction$annotation).
#' If only focusing on `Secreted Signaling`, the `contact-dependent` manner will be not used except for setting `contact.dependent.forced = TRUE`.
#' @param contact.range The interaction range (Unit: microns) to restrict the contact-dependent signaling when `contact.dependent = TRUE`.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine signaling can only happens when the two cells are contact to each other.
#'
#' Typically, `contact.range = 10`, which is a typical human cell size. However, for low-resolution spatial data such as 10X visium, it should be the cell center-to-center distance (i.e., `contact.range = 100` for visium data).  The function `computeCellDistance` can compute the center-to-center distance.
#'
#' @param contact.knn.k Number of neighbors to restrict the contact-dependent signaling within the neatest neighbors when `contact.dependent = TRUE`. By default, CellChat uses `contact.range` to restrict the contact-dependent signaling; however, users can also provide a value of `contact.knn.k`, in order to determine interacting cell pairs based on the k-nearest neighbors (knn).
#' For 10X visium, contact.knn.k = 6. For other spatial technologies, this value may be hard to determine because the sequenced cells/spots are usually not regularly arranged.
#' @param do.symmetric Whether converting the adjacent matrix into symmetric one when determining spatially proximal cell groups. Default is TRUE, indicating that if adj(i,j) or adj(j,i) is zero, then both are zeros.
#'
#' @param contact.dependent.forced Whether forcing to use the `contact-dependent` manner for inference signaling for all L-R pairs including secreted signaling. Users can set `contact.dependent.forced = TRUE` if also preferring interactions within a contact manner for `Secreted Signaling`.
#'
#' @param nboot Threshold of p-values
#' @param seed.use Set a random seed. By default, set the seed to 1.
#' @param Kh Parameter in Hill function
#' @param n Parameter in Hill function
#'
#'
#' @importFrom future.apply future_sapply
#' @importFrom progressr progressor
#' @importFrom stats aggregate
#' @importFrom Matrix crossprod
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return A CellChat object with updated slot 'net':
#'
#' object@net$prob is the inferred communication probability (strength) array, where the first, second and third dimensions represent a source, target and ligand-receptor pair, respectively.
#'
#' USER can access all the inferred cell-cell communications using the function 'subsetCommunication(object)', which returns a data frame.
#'
#' object@net$pval is the corresponding p-values of each interaction
#'
#' @export
#'
computeCommunProb <- function(object, type = c("triMean", "truncatedMean","thresholdedMean", "median"), trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.01, k.min = 10, 
			      contact.dependent = TRUE, contact.range = NULL, contact.knn.k = NULL, contact.dependent.forced = FALSE, do.symmetric = TRUE,
                              nboot = 100, seed.use = 1L, Kh = 0.5, n = 1) {
  type <- match.arg(type)
  cat(type, "is used for calculating the average gene expression per cell group.", "\n")
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    thresholdedMean = function(x) thresholdedMean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  if (raw.use) {
    data <- as.matrix(object@data.signaling)
  } else {
    data <- as.matrix(object@data.smooth)
  }
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    if (length(unique(LR.use$annotation)) > 1) {
      LR.use$annotation <- factor(LR.use$annotation, levels = c("Secreted Signaling","ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"))
      LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
      LR.use$annotation <- as.character(LR.use$annotation)
    }
    pairLR.use <- LR.use
  }
  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor

  ptm = Sys.time()

  pairLRsig <- pairLR.use
  group <- object@idents
  geneL <- as.character(pairLRsig$ligand)
  geneR <- as.character(pairLRsig$receptor)
  nLR <- nrow(pairLRsig)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!
         You may need to drop unused levels using 'droplevels' function. e.g.,
         `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`")
  }

  data.use <- data/max(data)
  nC <- ncol(data.use)

  # compute the average expression per group
  data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)
  data.use.avg <- t(data.use.avg[,-1])
  colnames(data.use.avg) <- levels(group)
  # compute the expression of ligand or receptor
  dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
  dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
  # take account into the effect of co-activation and co-inhibition receptors
  dataRavg.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "A")
  dataRavg.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avg, pairLRsig, type = "I")
  dataRavg <- dataRavg * dataRavg.co.A.receptor/dataRavg.co.I.receptor

  dataLavg2 <- t(replicate(nrow(dataLavg), as.numeric(table(group))/nC))
  dataRavg2 <- dataLavg2

  # compute the expression of agonist and antagonist
  index.agonist <- which(!is.na(pairLRsig$agonist) & pairLRsig$agonist != "")
  index.antagonist <- which(!is.na(pairLRsig$antagonist) & pairLRsig$antagonist != "")
  # quantify the communication probability

  # compute the spatial constraint
  if (object@options$datatype != "RNA") {
    data.spatial <- object@images$coordinates
    if ("spatial.factors" %in% names(object@images)) {
      ratio <- object@images$spatial.factors$ratio
      tol <- object@images$spatial.factors$tol
    } else {
      stop("`object@images$spatial.factors` is missing. Please update the object via `updateCellChat`! \n")
    }

    meta.t = data.frame(group = group, samples = object@meta$samples, row.names = rownames(object@meta))
    res <- computeRegionDistance(coordinates = data.spatial, meta = meta.t, interaction.range = interaction.range, ratio = ratio, tol = tol, k.min = k.min, contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k)
    d.spatial <- res$d.spatial # NaN if no nearby cell pairs
    adj.contact <- res$adj.contact # zeros if no nearby cell pairs
    if (distance.use) {
      print(paste0('>>> Run CellChat on spatial transcriptomics data using distances as constraints of the computed communication probability <<< [', Sys.time(),']'))
      d.spatial <- d.spatial * scale.distance
      diag(d.spatial) <- NaN
      d.min <- min(d.spatial, na.rm = TRUE)
      if (d.min < 1) {
        cat("The suggested minimum value of scaled distances is in [1,2], and the calculated value here is ", d.min,"\n")
        stop("Please increase the value of `scale.distance` and use a value that is slighly smaller than ", format(1/d.min, digits = 2) ,"\n")
      }
      P.spatial <- 1/d.spatial
      P.spatial[is.na(d.spatial)] <- 0
      diag(P.spatial) <- max(P.spatial) # if this value is 1, the self-connections will have more larger weight.
      d.spatial <- d.spatial/scale.distance # This is only for saving the data
    } else {
      print(paste0('>>> Run CellChat on spatial transcriptomics data without distance values as constraints of the computed communication probability <<< [', Sys.time(),']'))
      P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
      P.spatial[is.na(d.spatial)] <- 0 # diagonal is 1
    }

  } else {
    print(paste0('>>> Run CellChat on sc/snRNA-seq data <<< [', Sys.time(),']'))
    d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
    P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
    adj.contact <- matrix(1, nrow = numCluster, ncol = numCluster)
    contact.dependent = FALSE; contact.dependent.forced = FALSE; contact.range = NULL; contact.knn.k = NULL;
    distance.use = NULL; interaction.range = NULL; ratio = NULL; tol = NULL; k.min = NULL;
  }

  if (object@options$datatype == "RNA") {
    nLR1 <- nLR
  } else {
    if (contact.dependent.forced == TRUE) {
      cat("Force to run CellChat in a `contact-dependent` manner for all L-R pairs including secreted signaling.\n")
      P.spatial <- P.spatial * adj.contact
      nLR1 <- nLR
    } else { # contact.dependent.forced == F
      if (contact.dependent == TRUE && length(unique(pairLRsig$annotation)) > 0) {
        if (all(unique(pairLRsig$annotation) %in% c("Cell-Cell Contact"))) {
          cat("All the input L-R pairs are `Cell-Cell Contact` signaling. Run CellChat in a contact-dependent manner. \n")
          P.spatial <- P.spatial * adj.contact
          nLR1 <- nLR
        } else if (all(unique(pairLRsig$annotation) %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
          cat("Molecules of the input L-R pairs are diffusible. Run CellChat in a diffusion manner based on the `interaction.range`.\n")
          nLR1 <- nLR
        } else {
          cat("The input L-R pairs have both secreted signaling and contact-dependent signaling. Run CellChat in a contact-dependent manner for `Cell-Cell Contact` signaling, and in a diffusion manner based on the `interaction.range` for other L-R pairs. \n")
          nLR1 <- max(which(pairLRsig$annotation %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
        }
      } else { # contact.dependent == F or there is no `annotation` column in the database
        cat("Run CellChat in a diffusion manner based on the `interaction.range` for all L-R pairs. Setting `contact.dependent = TRUE` if preferring a contact-dependent manner for `Cell-Cell Contact` signaling. \n")
        nLR1 <- nLR
      }
    }
  }

  Prob <- array(0, dim = c(numCluster,numCluster,nLR))
  Pval <- array(0, dim = c(numCluster,numCluster,nLR))

  set.seed(seed.use)
  permutation <- replicate(nboot, sample.int(nC, size = nC))
  p <- progressr::progressor(nboot)
  data.use.avg.boot <- future.apply::future_sapply(
    X = 1:nboot,
    FUN = function(nE) {
      p()
      groupboot <- group[permutation[, nE]]
      data.use.avgB <- aggregate(t(data.use), list(groupboot), FUN = FunMean)
      data.use.avgB <- t(data.use.avgB[,-1])
      return(data.use.avgB)
    },
    future.seed = TRUE,
    simplify = FALSE
  )
  pb <- txtProgressBar(min = 0, max = nLR, style = 3, file = stderr())

  for (i in 1:nLR) {
    # ligand/receptor
    dataLR <- Matrix::crossprod(matrix(dataLavg[i,], nrow = 1), matrix(dataRavg[i,], nrow = 1))
    P1 <- dataLR^n/(Kh^n + dataLR^n)
    P1_Pspatial <- P1*P.spatial
    if (sum(P1_Pspatial) == 0) {
      Pnull = P1_Pspatial
      Prob[ , , i] <- Pnull
      p = 1
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    } else {
      if (i > nLR1) {
        P.spatial <- P.spatial * adj.contact
      }
      # agonist and antagonist
      if (is.element(i, index.agonist)) {
        data.agonist <- computeExpr_agonist(data.use = data.use.avg, pairLRsig, cofactor_input, index.agonist = i, Kh = Kh,  n = n)
        P2 <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
      } else {
        P2 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      if (is.element(i, index.antagonist)) {
        data.antagonist <- computeExpr_antagonist(data.use = data.use.avg, pairLRsig, cofactor_input,  index.antagonist = i, Kh = Kh,  n = n)
        P3 <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
      } else {
        P3 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }
      # number of cells
      if (population.size) {
        P4 <- Matrix::crossprod(matrix(dataLavg2[i,], nrow = 1), matrix(dataRavg2[i,], nrow = 1))
      } else {
        P4 <- matrix(1, nrow = numCluster, ncol = numCluster)
      }

      # Pnull = P1*P2*P3*P4
      Pnull = P1*P2*P3*P4*P.spatial
      Prob[ , , i] <- Pnull

      Pnull <- as.vector(Pnull)

      p <- progressr::progressor(nboot)
      Pboot <- future.apply::future_sapply(
        X = 1:nboot,
        FUN = function(nE) {
          p()
          data.use.avgB <- data.use.avg.boot[[nE]]
          dataLavgB <- computeExpr_LR(geneL[i], data.use.avgB, complex_input)
          dataRavgB <- computeExpr_LR(geneR[i], data.use.avgB, complex_input)
          # take account into the effect of co-activation and co-inhibition receptors
          dataRavgB.co.A.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avgB, pairLRsig[i, , drop = FALSE], type = "A")
          dataRavgB.co.I.receptor <- computeExpr_coreceptor(cofactor_input, data.use.avgB, pairLRsig[i, , drop = FALSE], type = "I")
          dataRavgB <- dataRavgB * dataRavgB.co.A.receptor/dataRavgB.co.I.receptor
          dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
          P1.boot <- dataLRB^n/(Kh^n + dataLRB^n)
          # agonist and antagonist
          if (is.element(i, index.agonist)) {
            data.agonist <- computeExpr_agonist(data.use = data.use.avgB, pairLRsig, cofactor_input, index.agonist = i, Kh = Kh,  n = n)
            P2.boot <- Matrix::crossprod(matrix(data.agonist, nrow = 1))
          } else {
            P2.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
          }
          if (is.element(i, index.antagonist)) {
            data.antagonist <- computeExpr_antagonist(data.use = data.use.avgB, pairLRsig, cofactor_input, index.antagonist = i, Kh = Kh,  n= n)
            P3.boot <- Matrix::crossprod(matrix(data.antagonist, nrow = 1))
          } else {
            P3.boot <- matrix(1, nrow = numCluster, ncol = numCluster)
          }

          if (population.size) {
            groupboot <- group[permutation[, nE]]
            dataLavg2B <- as.numeric(table(groupboot))/nC
            dataLavg2B <- matrix(dataLavg2B, nrow = 1)
            dataRavg2B <- dataLavg2B
            P4.boot = Matrix::crossprod(dataLavg2B, dataRavg2B)
          } else {
            P4.boot = matrix(1, nrow = numCluster, ncol = numCluster)
          }

          #  Pboot = P1.boot*P2.boot*P3.boot*P4.boot
          Pboot = P1.boot*P2.boot*P3.boot*P4.boot*P.spatial
          return(as.vector(Pboot))
        }
      )
      Pboot <- matrix(unlist(Pboot), nrow=length(Pnull), ncol = nboot, byrow = FALSE)
      nReject <- rowSums(Pboot - Pnull > 0)
      p = nReject/nboot
      Pval[, , i] <- matrix(p, nrow = numCluster, ncol = numCluster, byrow = FALSE)
    }
    setTxtProgressBar(pb = pb, value = i)
  }
  close(con = pb)
  Pval[Prob == 0] <- 1
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLRsig))
  dimnames(Pval) <- dimnames(Prob)
  net <- list("prob" = Prob, "pval" = Pval)
  execution.time = Sys.time() - ptm
  object@options$run.time <- as.numeric(execution.time, units = "secs")

  object@options$parameter <- list(type.mean = type, trim = trim, raw.use = raw.use, population.size = population.size,  nboot = nboot, seed.use = seed.use, Kh = Kh, n = n,
                                   distance.use = distance.use, interaction.range = interaction.range, ratio = ratio, tol = tol, k.min = k.min,
                                   contact.dependent = contact.dependent, contact.range = contact.range, contact.knn.k = contact.knn.k, contact.dependent.forced = contact.dependent.forced
                                   )
  if (object@options$datatype != "RNA") {
    object@images$distance <- d.spatial
  }
  object@net <- net
  print(paste0('>>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [', Sys.time(),']'))
  return(object)
}


#' Compute the communication probability on signaling pathway level by summarizing all related ligands/receptors
#'
#' @param object CellChat object
#' @param net A list from object@net; If net = NULL, net = object@net
#' @param pairLR.use A dataframe giving the ligand-receptor interactions; If pairLR.use = NULL, pairLR.use = object@LR$LRsig
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return A CellChat object with updated slot 'netP':
#'
#' object@netP$prob is the communication probability array on signaling pathway level; USER can convert this array to a data frame using the function 'reshape2::melt()',
#'
#' e.g., `df.netP <- reshape2::melt(object@netP$prob, value.name = "prob"); colnames(df.netP)[1:3] <- c("source","target","pathway_name")` or access all significant interactions using the function \code{\link{subsetCommunication}}
#'
#' object@netP$pathways list all the signaling pathways with significant communications.
#'
#' From version >= 1.1.0, pathways are ordered based on the total communication probabilities. NB: pathways with small total communication probabilities might be also very important since they might be specifically activated between only few cell types.
#'
#' @export
#'
computeCommunProbPathway <- function(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  prob <- net$prob
  prob[net$pval > thresh] <- 0

  LR <- dimnames(prob)[[3]]
  LR.sig <- LR[apply(prob, 3, sum) != 0]

  pathways <- unique(pairLR.use$pathway_name)
  group <- factor(pairLR.use$pathway_name, levels = pathways)
  prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
  pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]
  prob.pathways.sig <- prob.pathways[,,pathways.sig, drop = FALSE]
  idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing=TRUE, index.return = TRUE)$ix
  pathways.sig <- pathways.sig[idx]
  prob.pathways.sig <- prob.pathways.sig[, , idx]

  if (is.null(object)) {
    netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
    return(netP)
  } else {
    object@net$LRs <- LR.sig
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    return(object)
  }
}


#' Calculate the aggregated network by counting the number of links or summarizing the communication probability
#'
#' @param object CellChat object
#' @param sources.use,targets.use,signaling,pairLR.use Please check the description in function \code{\link{subsetCommunication}}
#' @param remove.isolate whether removing the isolate cell groups without any interactions when applying \code{\link{subsetCommunication}}
#' @param thresh threshold of the p-value for determining significant interaction
#' @param return.object whether return an updated CellChat object
#' @importFrom  dplyr group_by summarize groups
#' @importFrom stringr str_split
#'
#' @return Return an updated CellChat object:
#'
#' `object@net$count` is a matrix: rows and columns are sources and targets respectively, and elements are the number of interactions between any two cell groups. USER can convert a matrix to a data frame using the function `reshape2::melt()`
#'
#' `object@net$weight` is also a matrix containing the interaction weights between any two cell groups
#'
#' `object@net$sum` is deprecated. Use `object@net$weight`
#'
#' @export
#'
aggregateNet <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE) {
  net <- object@net
  if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
    prob <- net$prob
    pval <- net$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net$count <- apply(prob > 0, c(1,2), sum)
    net$weight <- apply(prob, c(1,2), sum)
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  } else {
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, "_", simplify = T)
    df.net2$source <- as.character(a[, 1])
    df.net2$target <- as.character(a[, 2])
    cells.level <- levels(object@idents)
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }

    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$count <- count
    net$weight <- prob
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  }
  if (return.object) {
    object@net <- net
    return(object)
  } else {
    return(net)
  }

}


#' Compute averaged expression values for each cell group
#'
#' @param object CellChat object
#' @param features a char vector giving the used features. default use all features
#' @param group.by cell group information; default is `object@idents` when input is a single object and `object@idents$joint` when input is a merged object; otherwise it should be one of the column names of the meta slot
#' @param type methods for computing the average gene expression per cell group.
#'
#' By default = "triMean", defined as a weighted average of the distribution's median and its two quartiles (https://en.wikipedia.org/wiki/Trimean);
#'
#' When setting `type = "truncatedMean"`, a value should be assigned to 'trim'. See the function `base::mean`.
#'
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.
#' @param slot.name the data in the slot.name to use
#' @param data.use a customed data matrix. Default: data.use = NULL and the expression matrix in the 'slot.name' is used
#'
#' @return Returns a matrix with genes as rows, cell groups as columns.

#' @export
#'
computeAveExpr <- function(object, features = NULL, group.by = NULL, type = c("triMean", "truncatedMean", "median"), trim = NULL,
                           slot.name = c("data.signaling", "data"), data.use = NULL) {
  type <- match.arg(type)
  slot.name <- match.arg(slot.name)
  FunMean <- switch(type,
                    triMean = triMean,
                    truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                    median = function(x) median(x, na.rm = TRUE))

  if (is.null(data.use)) {
    data.use <- slot(object, slot.name)
  }
  if (is.null(features)) {
    features.use <- row.names(data.use)
  } else {
    features.use <- intersect(features, row.names(data.use))
  }
  data.use <- data.use[features.use, , drop = FALSE]
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
  # compute the average expression per group
  data.use.avg <- aggregate(t(data.use), list(labels), FUN = FunMean)
  data.use.avg <- t(data.use.avg[,-1])
  rownames(data.use.avg) <- features.use
  colnames(data.use.avg) <- levels(labels)
  return(data.use.avg)
}



#' Compute the expression of complex in individual cells using geometric mean
#' @param complex_input the complex_input from CellChatDB
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex the names of complex
#' @return
#' @importFrom dplyr select starts_with
#' @importFrom future.apply future_sapply
#' @importFrom progressr progressor
#' @export
computeExpr_complex <- function(complex_input, data.use, complex) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))
  nrun <- nrow(Rsubunits)
  p <- progressr::progressor(nrun)
  data.complex <- future.apply::future_sapply(
    X = 1:nrun,
    FUN = function(x) {
      p()
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      return(geometricMean(data.use[RsubunitsV, , drop = FALSE]))
    }
  )
  data.complex <- t(data.complex)
  return(data.complex)
}

# Compute the average expression of complex per cell group using geometric mean
# @param complex_input the complex_input from CellChatDB
# @param data.use data matrix (rows are genes and columns are cells)
# @param complex the names of complex
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom dplyr select starts_with
# @importFrom future.apply future_sapply
# @importFrom progressr progressor
# @export
.computeExprGroup_complex <- function(complex_input, data.use, complex, group, FunMean) {
  Rsubunits <- complex_input[complex,] %>% dplyr::select(starts_with("subunit"))
  nrun <- nrow(Rsubunits)
  p <- progressr::progressor(nrun)
  data.complex <- future.apply::future_sapply(
    X = 1:nrun,
    FUN = function(x) {
      p()
      RsubunitsV <- unlist(Rsubunits[x,], use.names = F)
      RsubunitsV <- RsubunitsV[RsubunitsV != ""]
      RsubunitsV <- intersect(RsubunitsV, rownames(data.use))
      if (length(RsubunitsV) > 1) {
        data.avg <- aggregate(t(data.use[RsubunitsV, ,drop = FALSE]), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else if (length(RsubunitsV) == 1) {
        data.avg <- aggregate(matrix(data.use[RsubunitsV,], ncol = 1), list(group), FUN = FunMean)
        data.avg <- t(data.avg[,-1])
      } else {
        data.avg = matrix(0, nrow = 1, ncol = length(unique(group)))
      }
      return(geometricMean(data.avg))
    }
  )
  data.complex <- t(data.complex)
  return(data.complex)
}

#' Compute the expression of ligands or receptors using geometric mean
#' @param geneLR a char vector giving a set of ligands or receptors
#' @param data.use data matrix (row are genes and columns are cells or cell groups)
#' @param complex_input the complex_input from CellChatDB
#' @param group a factor defining the cell groups; If NULL, compute the expression of ligands or receptors in individual cells; otherwise, compute the average expression of ligands or receptors per cell group
#' @param FunMean the function for computing average expression per cell group
#' @return
#' @export
computeExpr_LR <- function(geneLR, data.use, complex_input){
  nLR <- length(geneLR)
  numCluster <- ncol(data.use)
  index.singleL <- which(geneLR %in% rownames(data.use))
  dataL1avg <- data.use[geneLR[index.singleL],]
  dataLavg <- matrix(nrow = nLR, ncol = numCluster)
  dataLavg[index.singleL,] <- dataL1avg
  index.complexL <- setdiff(1:nLR, index.singleL)
  if (length(index.complexL) > 0) {
    complex <- geneLR[index.complexL]
    data.complex <- computeExpr_complex(complex_input, data.use, complex)
    dataLavg[index.complexL,] <- data.complex
  }
  return(dataLavg)
}


#' Modeling the effect of coreceptor on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig a data frame giving ligand-receptor interactions
#' @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
#' @return
#' @importFrom future.apply future_sapply
#' @importFrom progressr progressor
#' @export
computeExpr_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I")) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {
    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    nrun <- nrow(coreceptor.ind)
    p <- progressr::progressor(nrun)
    data.coreceptor.ind <- future.apply::future_sapply(
      X = 1:nrun,
      FUN = function(x) {
        p()
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) == 1) {
          return(1 + data.use[coreceptor.indV, ])
        } else if (length(coreceptor.indV) > 1) {
          return(apply(1 + data.use[coreceptor.indV, ], 2, prod))
        } else {
          return(matrix(1, nrow = 1, ncol = ncol(data.use)))
        }
      }
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = ncol(data.use))
  }
  return(data.coreceptor)
}

# Modeling the effect of coreceptor on the ligand-receptor interaction
#
# @param data.use data matrix
# @param cofactor_input the cofactor_input from CellChatDB
# @param pairLRsig a data frame giving ligand-receptor interactions
# @param type when type == "A", computing expression of co-activation receptor; when type == "I", computing expression of co-inhibition receptor.
# @param group a factor defining the cell groups
# @param FunMean the function for computing mean expression per group
# @return
# @importFrom future.apply future_sapply
# @importFrom progressr progressor
# @export
.computeExprGroup_coreceptor <- function(cofactor_input, data.use, pairLRsig, type = c("A", "I"), group, FunMean) {
  type <- match.arg(type)
  if (type == "A") {
    coreceptor.all = pairLRsig$co_A_receptor
  } else if (type == "I"){
    coreceptor.all = pairLRsig$co_I_receptor
  }
  index.coreceptor <- which(!is.na(coreceptor.all) & coreceptor.all != "")
  if (length(index.coreceptor) > 0) {
    coreceptor <- coreceptor.all[index.coreceptor]
    coreceptor.ind <- cofactor_input[coreceptor, grepl("cofactor" , colnames(cofactor_input) )]
    nrun <- nrow(coreceptor.ind)
    p <- progressr::progressor(nrun)
    data.coreceptor.ind <- future.apply::future_sapply(
      X = 1:nrun,
      FUN = function(x) {
        p()
        coreceptor.indV <- unlist(coreceptor.ind[x,], use.names = F)
        coreceptor.indV <- coreceptor.indV[coreceptor.indV != ""]
        coreceptor.indV <- intersect(coreceptor.indV, rownames(data.use))
        if (length(coreceptor.indV) > 1) {
          data.avg <- aggregate(t(data.use[coreceptor.indV,]), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(apply(1 + data.avg, 2, prod))
        } else if (length(coreceptor.indV) == 1) {
          data.avg <- aggregate(matrix(data.use[coreceptor.indV,], ncol = 1), list(group), FUN = FunMean)
          data.avg <- t(data.avg[,-1])
          return(1 + data.avg)
        } else {
          return(matrix(1, nrow = 1, ncol = length(unique(group))))
        }
      }
    )
    data.coreceptor.ind <- t(data.coreceptor.ind)
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
    data.coreceptor[index.coreceptor,] <- data.coreceptor.ind
  } else {
    data.coreceptor <- matrix(1, nrow = length(coreceptor.all), ncol = length(unique(group)))
  }

  return(data.coreceptor)
}

#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
#' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_agonist <- function(data.use, pairLRsig, cofactor_input, group, index.agonist, Kh, FunMean, n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
#' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExprGroup_antagonist <- function(data.use, pairLRsig, cofactor_input, group, index.antagonist, Kh, FunMean, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    data.avg <- t(data.avg[,-1])
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
  }
  return(data.antagonist)
}


#' Modeling the effect of agonist on the ligand-receptor interaction
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.agonist the index of agonist in the database
#' @param Kh a parameter in Hill function
#' @param FunMean the function for computing mean expression per group
#' @param n Hill coefficient
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_agonist <- function(data.use, pairLRsig, cofactor_input, index.agonist, Kh,  n) {
  agonist <- pairLRsig$agonist[index.agonist]
  agonist.ind <- cofactor_input[agonist, grepl("cofactor" , colnames(cofactor_input))]
  agonist.indV <- unlist(agonist.ind, use.names = F)
  agonist.indV <- agonist.indV[agonist.indV != ""]
  agonist.indV <- intersect(agonist.indV, rownames(data.use))
  if (length(agonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[agonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- 1 + data.avg^n/(Kh^n + data.avg^n)
  } else if (length(agonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[agonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[agonist.indV,, drop = FALSE]
    data.agonist <- apply(1 + data.avg^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.agonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.agonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.agonist)
}

#' Modeling the effect of antagonist on the ligand-receptor interaction
#'
#' @param data.use data matrix
#' @param cofactor_input the cofactor_input from CellChatDB
#' @param pairLRsig the L-R interactions
#' @param group a factor defining the cell groups
#' @param index.antagonist the index of antagonist in the database
#' @param Kh a parameter in Hill function
#' @param n Hill coefficient
#' @param FunMean the function for computing mean expression per group
#' @return
#' @export
#' @importFrom stats aggregate
computeExpr_antagonist <- function(data.use, pairLRsig, cofactor_input, index.antagonist, Kh, n) {
  antagonist <- pairLRsig$antagonist[index.antagonist]
  antagonist.ind <- cofactor_input[antagonist, grepl( "cofactor" , colnames(cofactor_input) )]
  antagonist.indV <- unlist(antagonist.ind, use.names = F)
  antagonist.indV <- antagonist.indV[antagonist.indV != ""]
  antagonist.indV <- intersect(antagonist.indV, rownames(data.use))
  if (length(antagonist.indV) == 1) {
    # data.avg <- aggregate(matrix(data.use[antagonist.indV,], ncol = 1), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- Kh^n/(Kh^n + data.avg^n)
  } else if (length(antagonist.indV) > 1) {
    # data.avg <- aggregate(t(data.use[antagonist.indV,]), list(group), FUN = FunMean)
    # data.avg <- t(data.avg[,-1])
    data.avg <- data.use[antagonist.indV,, drop = FALSE]
    data.antagonist <- apply(Kh^n/(Kh^n + data.avg^n), 2, prod)
  } else {
    # data.antagonist = matrix(1, nrow = 1, ncol = length(unique(group)))
    data.antagonist = matrix(1, nrow = 1, ncol = ncol(data.use))
  }
  return(data.antagonist)
}


#' Compute the geometric mean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return
#' @export
geometricMean <- function(x,na.rm=TRUE){
  if (is.null(nrow(x))) {
    exp(mean(log(x),na.rm=na.rm))
  } else {
    exp(apply(log(x),2,mean,na.rm=na.rm))
  }
}


#' Compute the Tukey's trimean
#' @param x a numeric vector
#' @param na.rm whether remove na
#' @return
#' @importFrom stats quantile
#' @export
triMean <- function(x, na.rm = TRUE) {
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
}

#' Compute the average expression per cell group when the percent of expressing cells per cell group larger than a threshold
#' @param x a numeric vector
#' @param trim the percent of expressing cells per cell group to be considered as zero
#' @param na.rm whether remove na
#' @return
#' @importFrom Matrix nnzero
#' @export
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x)/length(x)
  if (percent < trim) {
    return(0)
  } else {
    return(mean(x, na.rm = na.rm))
  }
}

#' Filter cell-cell communication if there are only few number of cells in certain cell groups or inconsistent cell-cell communication across samples
#'
#' @param object CellChat object
#' @param min.cells The minmum number of cells required in each cell group for cell-cell communication
#' @param min.samples The minmum number of samples required for consistent cell-cell communication across samples (that is an interaction present in at least `min.samples` samples) when mutiple samples/replicates/batches are merged as an input for CellChat analysis.
#' @param rare.keep Whether to keep the interactions associated with the rare populations when min.samples >= 2. When a rare population is identified in the merged samples (say 15 cells in this rare population from two samples), it is likely to filter out the interactions associated with this rare population when setting min.samples >= 2. Setting `rare.keep = TRUE` to retain the identified interactions associated with this rare population.
#' @param nonFilter.keep Whether to keep the non-filtered cell-cell communication in the CellChat object. This is useful for avoiding re-running `computeCommunProb` if you want to adjust the parameters when running `filterCommunication`.
#' @return CellChat object with an updated slot net
#' @export
#'
filterCommunication <- function(object, min.cells = 10, min.samples = NULL, rare.keep = FALSE, nonFilter.keep = FALSE) {
  net <- object@net
  if (nonFilter.keep == TRUE) {
    cat("The non-filtered cell-cell communication is stored in `object@net$prob.nonFilter` and `object@net$pval.nonFilter`. \n")
    object@net$prob.nonFilter <- net$prob
    object@net$pval.nonFilter <- net$pval
  }
  num.interaction0 <- sum(net$prob > 0)
  cell.excludes <- which(as.numeric(table(object@idents)) <= min.cells)
  if (length(cell.excludes) > 0) {
    cat("The cell-cell communication related with the following cell groups are excluded due to the few number of cells: ", toString(levels(object@idents)[cell.excludes]), "!",'\t')
    net$prob[cell.excludes,,] <- 0
    net$prob[,cell.excludes,] <- 0
    num.interaction1 <- sum(net$prob > 0)
    pct.dicrease <- scales::percent((num.interaction0-num.interaction1)/num.interaction0, accuracy = .1)
    cat(paste0(pct.dicrease, " interactions are removed!",'\n'))
  } else {
    num.interaction1 <- num.interaction0
  }

  sample.info <- object@meta$samples
  sample.id <- levels(sample.info)
  if (is.null(min.samples)) {
    min.samples <- 1
  } else if (min.samples > length(sample.id)) {
    stop(paste0("There are only ", length(sample.id), " samples in the data. Please change the value of `min.samples`! "))
  }
  if (length(sample.id) >= 2 & min.samples >= 2) {
    if (object@options$parameter$raw.use == TRUE) {
      data <- as.matrix(object@data.signaling)
    } else {
      if ("data.smooth" %in% methods::slotNames(object) == FALSE) {
        stop("`object@data.smooth` is missing. Please update the CellChat object via `updateCellChat`! \n")
      }
      data <- as.matrix(object@data.smooth)
    }
    data.use <- data/max(data)
    group <- object@idents
    type <- object@options$parameter$type.mean
    trim <- object@options$parameter$trim
    FunMean <- switch(type,
                      triMean = triMean,
                      truncatedMean = function(x) mean(x, trim = trim, na.rm = TRUE),
                      thresholdedMean = function(x) thresholdedMean(x, trim = trim, na.rm = TRUE),
                      median = function(x) median(x, na.rm = TRUE))

    LR <- dimnames(net$prob)[[3]]
    idx.nonzero <- which(apply(net$prob, 3, sum) != 0)
    LR.nonzero <- LR[idx.nonzero] # only examine the L-R pairs with nonzero communication probabilities.

    interaction_input <- object@DB$interaction
    complex_input <- object@DB$complex
    geneIfo <- object@DB$geneInfo
    idx <- match(LR.nonzero, interaction_input$interaction_name)
    geneL <- as.character(interaction_input$ligand[idx])
    geneR <- as.character(interaction_input$receptor[idx])

    geneLR <- c(unique(geneL), unique(geneR))
    geneLR <- extractGeneSubset(geneLR, complex_input, geneIfo)
    data.use <- data.use[rownames(data.use) %in% geneLR, ]

    score.LR <- array(0, dim = c(nlevels(group),nlevels(group),length(LR.nonzero), length(sample.id)))
    LR.nonzero.all <- c()
    cell.excludes.sample <- c()
    for (i in 1:length(sample.id)) {
      cell.use <- which(sample.info == sample.id[i])
      group.use <- group[cell.use]
      group.use <- droplevels(group.use)
      # get the rare populations with few cells in each sample
      cell.excludes.sample.i <- which(as.numeric(table(object@idents[cell.use])) <= min.cells)
      cell.excludes.sample <- c(cell.excludes.sample, cell.excludes.sample.i)
      # compute average expression per cell group
      data.use.i <- data.use[, cell.use]
      data.use.avg <- aggregate(t(data.use.i), list(group.use), FUN = FunMean)
      data.use.avg <- t(data.use.avg[,-1])
      group.exist <- which(levels(group) %in% unique(group.use))
      if (length(group.exist) < nlevels(group)) {
        data.use.avg.temp <- matrix(0, nrow = nrow(data.use), ncol = nlevels(group))
        data.use.avg.temp[ , group.exist] <- data.use.avg
        rownames(data.use.avg.temp) <- rownames(data.use.avg)
        data.use.avg <- data.use.avg.temp
      }
      colnames(data.use.avg) <- levels(group)
      # compute the average expression of ligand or receptor in each cell group
      dataLavg <- computeExpr_LR(geneL, data.use.avg, complex_input)
      dataRavg <- computeExpr_LR(geneR, data.use.avg, complex_input)
      # compute the interaction scores for each ligand-receptor pair based on their expression
      for (jj in 1:length(LR.nonzero)) { # It is not good to use parallel here because it will change the order of LR
        score.LR[,,jj,i] <- Matrix::crossprod(matrix(dataLavg[jj, ], nrow = 1), matrix(dataRavg[jj, ], nrow = 1))
      }
      if (length(cell.excludes.sample.i) > 0) {
        cat(paste0("The number of cells of the following cell groups in ", sample.id[i], " sample are less than ", min.cells, " cells: ",toString(levels(object@idents)[cell.excludes.sample.i]), "!",'\n'))
        score.LR[cell.excludes.sample.i, , , i] <- 0
        score.LR[ ,cell.excludes.sample.i, , i] <- 0
      }
      #LR.nonzero.all <- c(LR.nonzero.all, LR.nonzero[apply(score.LR[ , , , i], 3, sum) != 0])
    }
    #LR.nonzero.jointOnly <- setdiff(LR.nonzero, unique(LR.nonzero.all))

    # get the excluded cell groups that are not observed in the merged data, which is very possible for rare populations
    cell.excludes.sample <- unique(cell.excludes.sample)
    if (length(cell.excludes.sample) > 0) {
      cell.excludes.sample <- setdiff(cell.excludes.sample, cell.excludes)
    }

    score.LR[score.LR > 0] <- 1 # binarize the interaction score
    score.LR.consitent <- array(0, dim = c(nlevels(group),nlevels(group),length(LR.nonzero)))
    LR.inconsitent <- c()
    for (jj in 1:length(LR.nonzero)) {
      score.LR.sum <- apply(score.LR[ , , jj, ], c(1,2), sum) # elements 2 and 1 means consistent and inconsistent interactions across samples, respectively.
      # set communication probability to be zero for inconsistent interactions across samples
      if (sum((score.LR.sum > 0) * (score.LR.sum < min.samples)) > 0) {
        #LR.inconsitent <- c(LR.inconsitent, LR.nonzero[jj])
        score.LR.consitent <- (score.LR.sum >= min.samples) * 1
        if (rare.keep == TRUE & length(cell.excludes.sample) > 0) {
          score.LR.consitent[cell.excludes.sample, ] <- 1
          score.LR.consitent[ ,cell.excludes.sample] <- 1
        }
        net$prob[ , , LR.nonzero[jj]] <- net$prob[ , , LR.nonzero[jj]] * score.LR.consitent
      }
    }
    num.interaction2 <- sum(net$prob > 0)
    pct.dicrease <- scales::percent((num.interaction1-num.interaction2)/num.interaction1, accuracy = .1)
    cat(paste0(pct.dicrease, " interactions are removed due to their inconsistence across ", min.samples, " samples!",'\n'))
  }

  object@net <- net
  return(object)
}


#' Identify all the significant interactions (L-R pairs) from some cell groups to other cell groups
#'
#' @param object CellChat object
#' @param from a vector giving the index or the name of source cell groups
#' @param to a corresponding vector giving the index or the name of target cell groups. Note: The length of 'from' and 'to' must be the same, giving the corresponding pair of cell groups for communication.
#' @param bidirection whether show the bidirectional communication, i.e., both 'from'->'to' and 'to'->'from'.
#' @param pair.only whether only return ligand-receptor pairs without pathway names and communication strength
#' @param pairLR.use0 ligand-receptor pairs to use; default is all the significant interactions
#' @param thresh threshold of the p-value for determining significant interaction
#'
#' @return
#' @export
#'
identifyEnrichedInteractions <- function(object, from, to, bidirection = FALSE, pair.only = TRUE, pairLR.use0 = NULL, thresh = 0.05){
  pairwiseLR <- object@net$pairwiseRank
  if (is.null(pairwiseLR)) {
    stop("The interactions between pairwise cell groups have not been extracted!
         Please first run `object <- rankNetPairwise(object)`")
  }
  group.names.all <- names(pairwiseLR)
  if (!is.numeric(from)) {
    from <- match(from, group.names.all)
    if (sum(is.na(from)) > 0) {
      message("Some input cell group names in 'from' do not exist!")
      from <- from[!is.na(from)]
    }
  }
  if (!is.numeric(to)) {
    to <- match(to, group.names.all)
    if (sum(is.na(to)) > 0) {
      message("Some input cell group names in 'to' do not exist!")
      to <- to[!is.na(to)]
    }
  }
  if (length(from) != length(to)) {
    stop("The length of 'from' and 'to' must be the same!")
  }
  if (bidirection) {
    from2 <- c(from, to)
    to <- c(to, from)
    from <- from2
  }
  if (is.null(pairLR.use0)) {
    k <- 0
    pairLR.use0 <- list()
    for (i in 1:length(from)){
      pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
      idx <- pairwiseLR_ij$pval < thresh
      if (length(idx) > 0) {
        k <- k +1
        pairLR.use0[[k]] <- pairwiseLR_ij[idx,]
      }
    }
    pairLR.use0 <- do.call(rbind, pairLR.use0)
  }

  k <- 0
  pval <- matrix(nrow = length(rownames(pairLR.use0)), ncol = length(from))
  prob <- pval
  group.names <- c()
  for (i in 1:length(from)) {
    k <- k+1
    pairwiseLR_ij <- pairwiseLR[[from[i]]][[to[i]]]
    pairwiseLR_ij <- pairwiseLR_ij[rownames(pairLR.use0),]
    pval_ij <- pairwiseLR_ij$pval
    prob_ij <- pairwiseLR_ij$prob
    pval_ij[pval_ij > 0.05] = 1
    pval_ij[pval_ij > 0.01 & pval_ij <= 0.05] = 2
    pval_ij[pval_ij <= 0.01] = 3
    prob_ij[pval_ij ==1] <- 0
    pval[,k] <- pval_ij
    prob[,k] <- prob_ij
    group.names <- c(group.names, paste(group.names.all[from[i]], group.names.all[to[i]], sep = " - "))
  }
  prob[which(prob == 0)] <- NA
  # remove rows that are entirely NA
  pval <- pval[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  pairLR.use0 <- pairLR.use0[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  prob <- prob[rowSums(is.na(prob)) != ncol(prob), ,drop = FALSE]
  if (pair.only) {
    pairLR.use0 <- dplyr::select(pairLR.use0, ligand, receptor)
  }
  return(pairLR.use0)
}


#' Compute the region distance based on the spatial locations of each splot/cell of the spatial transcriptomics
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param meta a data frame including at least two columns named `group` and `samples`. `meta$group` is a factor vector defining the regions/labels of each cell/spot. `meta$samples` is a factor vector defining the sample labels of each dataset.
#' @param interaction.range The maximum interaction/diffusion range of ligands. This hard threshold is used to filter out the connections between spatially distant regions
#' @param ratio a numerical vector giving the conversion factor when converting spatial coordinates from Pixels or other units to Micrometers (i.e.,Microns).
#'
#' For example, setting `ratio = 0.18` indicates that 1 pixel equals 0.18um in the coordinates.
#' For 10X visium, it is the ratio of the theoretical spot size (i.e., 65um) over the number of pixels that span the diameter of a theoretical spot size in the full-resolution image (i.e., 'spot.size.fullres' in the 'scalefactors_json.json' file).
#' @param tol a numerical vector giving the tolerance factor to increase the robustness when comparing the center-to-center distance against the `interaction.range`. This can be the half value of cell/spot size in the unit of um.
#'
#' For example, for 10X visium, `tol` can be set as `65/2`; for slide-seq, `tol` can be set as `10/2`.
#' If the cell/spot size is not known, we provide a function `computeCellDistance` to compute the center-to-center distance. `tol` can be the the half value of the minimum center-to-center distance.
#' @param k.min the minimum number of interacting cell pairs required for defining adjacent cell groups
#' @param contact.dependent Whether determining spatially proximal cell groups based on either the contact.range or the k-nearest neighbors (knn). By default `contact.dependent = TRUE` when inferring contact-dependent and juxtacrine signaling (including ECM-Receptor and Cell-Cell Contact signaling classified in CellChatDB$interaction$annotation).
#' If only focusing on `Secreted Signaling`, the `contact.dependent` will be automatically set as FALSE except for `contact.dependent.forced = TRUE`.
#' @param contact.range The interaction range (Unit: microns) to restrict the contact-dependent signaling.
#' For spatial transcriptomics in a single-cell resolution, `contact.range` is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent and juxtacrine signaling can only happens when the two cells are contact to each other.
#'
#' Typically, `contact.range = 10`, which is a typical human cell size. However, for low-resolution spatial data such as 10X visium, it should be the cell center-to-center distance (i.e., `contact.range = 100` for visium data). The function `computeCellDistance` can compute the center-to-center distance.
#'
#' @param contact.knn.k Number of neighbors to restrict the contact-dependent signaling within the neatest neighbors. By default, CellChat uses `contact.range` to restrict the contact-dependent signaling; however, users can also provide a value of `contact.knn.k`, in order to determine spatially proximal cell groups based on the k-nearest neighbors (knn).
#' For 10X visium, contact.knn.k = 6. For other spatial technologies, this value may be hard to determine because the sequenced cells/spots are usually not regularly arranged.
#' @param do.symmetric Whether converting the adjacent matrix into symmetric one when determining spatially proximal cell groups. Default is TRUE, indicating that if adj(i,j) or adj(j,i) is zero, then both are zeros.
#'
#' @importFrom BiocNeighbors queryKNN AnnoyParam
#' @return A list including a square matrix giving the pairwise region distances and an adjacent matrix indicating physically contacting cell groups based on either the contact.range or the k-nearest neighbors
#'
#' @export
computeRegionDistance <- function(coordinates, meta,
                                  interaction.range = NULL, ratio = NULL, tol = NULL, k.min = 10,
                                  contact.dependent = TRUE, contact.range = NULL, contact.knn.k = NULL, do.symmetric = TRUE
) {
  trim <- 0.1
  FunMean <- function(x) mean(x, trim = trim, na.rm = TRUE) # This is used for computing the average distance between two cell groups
  group <- meta$group
  numCluster <- nlevels(group)
  level.use <- levels(group)
  level.use <- level.use[level.use %in% unique(group)]
  samples <- meta$samples
  samples.use <- levels(samples)
  d.spatial <- array(NaN, dim = c(numCluster,numCluster,length(samples.use)))
  adj.spatial <- array(0, dim = c(numCluster,numCluster,length(samples.use)))
  adj.contact <- array(0, dim = c(numCluster,numCluster,length(samples.use)))
  adj.contact.knn <- array(0, dim = c(numCluster,numCluster,length(samples.use)))

  if (contact.dependent == TRUE & !is.null(contact.knn.k)) {
    ## find the k-nearest neighbors for each single cell
    # my.knn <- FNN::get.knn(coordinates, k = contact.knn.k)
    # nn.ranked <- my.knn$nn.index # this is a matrix with the size of nCell * contact.knn.k
    nn.ranked <- matrix(NA, nrow = nrow(coordinates), ncol = contact.knn.k)
    for (k in 1:length(samples.use)) {
      idx.k <- which(samples == samples.use[k])
      my.knn <- suppressWarnings(BiocNeighbors::findKNN(coordinates[idx.k, ], k = contact.knn.k, BNPARAM = BiocNeighbors::AnnoyParam(), get.index = TRUE))
      nn.ranked[idx.k, ] <- my.knn$index # this is a matrix with the size of nCell * contact.knn.k
    }
    k.min.contact <- k.min
  } else {
    nn.ranked <- matrix(1, nrow = nrow(coordinates), ncol = 1)
    k.min.contact <- -1 # this produces adj.contact.knn with all elements being 1
  }
  if (contact.dependent == TRUE) {
    if (is.null(contact.range) & is.null(contact.knn.k)) {
      stop("Please check the documentation of `computeCommunProb` and provide the value of either `contact.range` or `contact.knn.k`")
    }
  } else {
    contact.range <- 10000 # this produces adj.contact with all elements being 1
  }

  for (k in 1:length(samples.use)) {
    idx.k <- samples == samples.use[k]
    for (i in 1:numCluster) {
      for (j in 1:numCluster) {
        idx.i <- which((group == level.use[i]) & idx.k)
        idx.j <- which((group == level.use[j]) & idx.k)
        if (length(idx.i) == 0 | length(idx.j) == 0) {
          next # if one cell group is missing in one sample, just goes to next loop
        }
        data.spatial.i <- coordinates[idx.i, , drop = FALSE]
        data.spatial.j <- coordinates[idx.j, , drop = FALSE]
        # for each point in the i-th cell group, find its 1-nearest neighbor in the j-th cell group
        #qout <- suppressWarnings(BiocNeighbors::queryKNN(data.spatial.j, data.spatial.i, k = 1, BNPARAM = BiocNeighbors::KmknnParam(), get.index = TRUE))
        qout <- suppressWarnings(BiocNeighbors::queryKNN(data.spatial.j, data.spatial.i, k = 1, BNPARAM = BiocNeighbors::AnnoyParam(), get.index = TRUE))
        # qout$index is an one column matrix with length being `length(idx.i)`, which is the index of the 1-nearest neighbor in the j-th cell group defined by `idx.j`
        # qout$distance is an one column matrix with length being `length(idx.i)`, which is the distance to the 1-nearest neighbor in the j-th cell group defined by `idx.j`

        # conver the calculated distance into the distance in micrometers
        qout$distance <- qout$distance*ratio[k]
        # long-range distance
        idx <- qout$distance - interaction.range < tol[k]
        adj.spatial[i,j,k] <- (length(unique(qout$index[idx])) >= k.min) * 1
        # short-range distance based on contact.range
        idx2 <- qout$distance - contact.range < tol[k]
        adj.contact[i,j,k] <- (length(unique(qout$index[idx2])) >= k.min) * 1
        # short-range distance based on knn
        knn.i <- unique(as.vector(nn.ranked[idx.i, ]))
        #adj.contact.knn[i,j,k] <- (length(intersect(knn.i, idx.j)) >= k.min.contact) * 1
        adj.contact.knn[i,j,k] <- (length(intersect(knn.i, unique(qout$index[idx]))) >= k.min.contact) * 1 # knn within the long-range distance
        # computing the average distance between two cell groups
        d.spatial[i,j,k] <- FunMean(qout$distance) # since distances are positive values, different ways for computing the mean have little effects.

      }
    }
  }

  # merged spatial information from different samples
  d.spatial <- apply(d.spatial, c(1,2), function(x) mean(x, na.rm = TRUE))
  adj.spatial <- apply(adj.spatial, c(1,2), mean)
  adj.contact <- apply(adj.contact, c(1,2), mean)
  adj.contact.knn <- apply(adj.contact.knn, c(1,2), mean)
  # for multi-samples analysis, the following is needed
  adj.spatial[adj.spatial > 0] <- 1
  adj.contact[adj.contact > 0] <- 1
  adj.contact.knn[adj.contact.knn > 0] <- 1

  # make these adjacent matrix as symmetric
  if (do.symmetric) {
    adj.spatial <- adj.spatial * t(adj.spatial) # if one is zero, then both are zeros.
    adj.contact <- adj.contact * t(adj.contact) # if one is zero, then both are zeros.
    adj.contact.knn <- adj.contact.knn * t(adj.contact.knn) # if one is zero, then both are zeros.
  }
  d.spatial <- (d.spatial + t(d.spatial))/2

  # filter out the spatially distant cell groups
  adj.spatial[adj.spatial == 0] <- NaN
  d.spatial <- d.spatial * adj.spatial

  rownames(d.spatial) <- levels(group); colnames(d.spatial) <- levels(group)

  if (length(contact.knn.k) > 0) {
    adj.contact = adj.contact.knn
  }
  res <- list(d.spatial = d.spatial, adj.contact = adj.contact)
  return(res)

}

#' Compute cell-cell distance based on the spatial coordinates
#'
#' @param coordinates a data matrix in which each row gives the spatial locations/coordinates of each cell/spot
#' @param interaction.range The maximum interaction/diffusion range of ligands. This hard threshold is used to filter out the connections between spatially distant cells
#' @param ratio The conversion factor when converting spatial coordinates from Pixels or other units to Micrometers (i.e.,Microns).
#'
#' For example, setting `ratio = 0.18` indicates that 1 pixel equals 0.18um in the coordinates.
#' For 10X visium, it is the ratio of the theoretical spot size (i.e., 65um) over the number of pixels that span the diameter of a theoretical spot size in the full-resolution image (i.e., 'spot.size.fullres' in the 'scalefactors_json.json' file).
#' @param tol The tolerance factor to increase the robustness when comparing the center-to-center distance against the `interaction.range`. This can be the half value of cell/spot size in the unit of um.
#'
#' For example, for 10X visium, `tol` can be set as `65/2`; for slide-seq, `tol` can be set as `10/2`.
#' If the cell/spot size is not known, we provide a function `computeCellDistance` to compute the center-to-center distance. `tol` can be the the half value of the minimum center-to-center distance.
#'
#' @return an object of class "dist" giving the pairwise cell-cell distance
#' @export
#'
computeCellDistance <- function(coordinates, interaction.range = NULL, ratio = NULL, tol = NULL){
  if (ncol(coordinates) == 2) {
    colnames(coordinates) <- c("x_cent","y_cent")
  } else {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
  }

  d.spatial <- stats::dist(coordinates)
  if (!is.null(ratio)) {
    d.spatial <- d.spatial*ratio
  }

  if(!is.null(interaction.range) & !is.null(tol)){
    message("\n Apply a predefined spatial distance threshold based on the interaction length...")
    d.spatial[d.spatial > (interaction.range + tol)] <- NaN
  }
  return(d.spatial)
}


