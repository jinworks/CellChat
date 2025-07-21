#' Show the description of CellChatDB databse
#'
#' @param CellChatDB CellChatDB databse
#' @param nrow the number of rows in the plot
#' @importFrom dplyr group_by summarise n %>%
#'
#' @return
#' @export
#'
showDatabaseCategory <- function(CellChatDB, nrow = 1) {
  interaction_input <- CellChatDB$interaction
  geneIfo <- CellChatDB$geneInfo
  df <- interaction_input %>% group_by(annotation) %>% summarise(value=n())
  #df$group <- factor(df$annotation, levels = unique(df$annotation))
  df$group <- factor(df$annotation, levels = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling"))
  gg1 <- pieChart(df)
  binary <- (interaction_input$ligand %in% geneIfo$Symbol) & (interaction_input$receptor %in% geneIfo$Symbol)
  df <- data.frame(group = rep("Heterodimers", dim(interaction_input)[1]),stringsAsFactors = FALSE)
  df$group[binary] <- rep("Others",sum(binary),1)
  df <- df %>% group_by(group) %>% summarise(value=n())
  df$group <- factor(df$group, levels = c("Heterodimers","Others"))
  gg2 <- pieChart(df)

  kegg <- grepl("KEGG", interaction_input$evidence)
  df <- data.frame(group = rep("Literature", dim(interaction_input)[1]),stringsAsFactors = FALSE)
  df$group[kegg] <- rep("KEGG",sum(kegg),1)
  df <- df %>% group_by(group) %>% summarise(value=n())
  df$group <- factor(df$group, levels = c("KEGG","Literature"))
  gg3 <- pieChart(df)

  gg <- cowplot::plot_grid(gg1, gg2, gg3, nrow = nrow, align = "h", rel_widths = c(1, 1,1))
  return(gg)
}


#' Plot pie chart
#'
#' @param df a dataframe
#' @param label.size a character
#' @param color.use the name of the variable in CellChatDB interaction_input
#' @param title the title of plot
#' @import ggplot2
#' @importFrom scales percent
#' @importFrom dplyr arrange desc mutate
#' @importFrom ggrepel geom_text_repel
#' @return
#' @export
#'
pieChart <- function(df, label.size = 2.5, color.use = NULL, title = "") {
  df %>% arrange(dplyr::desc(value)) %>%
    mutate(prop = scales::percent(value/sum(value))) -> df

  gg <- ggplot(df, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)+theme_void() +
    ggrepel::geom_text_repel(aes(label = prop), size= label.size, show.legend = F, position = position_stack(vjust=0.5))
  #  ggrepel::geom_text_repel(aes(label = prop), size= label.size, show.legend = F, nudge_x = 0)
  gg <- gg + theme(legend.position="bottom", legend.direction = "vertical")

  if(!is.null(color.use)) {
    gg <- gg + scale_fill_manual(values=color.use)
    # gg <- gg + scale_color_manual(color.use)
  }

  if (!is.null(title)) {
    gg <- gg + guides(fill = guide_legend(title = title))
  }
  gg
}


#' Subset the ligand-receptor interactions for given specific signals in CellChatDB
#'
#' @param signaling a character vector
#' @param pairLR.use a dataframe containing ligand-receptor interactions
#' @param key the keyword to match
#' @param matching.exact whether perform exact matching
#' @param pair.only whether only return ligand-receptor pairs without cofactors
#' @importFrom future.apply future_sapply
#' @importFrom dplyr select
#' @return
#' @export
searchPair <- function(signaling = c(), pairLR.use, key = c("pathway_name","ligand"), matching.exact = FALSE, pair.only = TRUE) {
  key <- match.arg(key)
  pairLR = future.apply::future_sapply(
    X = 1:length(signaling),
    FUN = function(x) {
      if (!matching.exact) {
        index <- grep(signaling[x], pairLR.use[[key]])
      } else {
        index <- which(pairLR.use[[key]] %in% signaling[x])
      }
      if (length(index) > 0) {
        if (pair.only) {
          pairLR <- dplyr::select(pairLR.use[index, ], interaction_name, pathway_name, ligand, receptor)
        } else {
          pairLR <- pairLR.use[index, ]
        }
        return(pairLR)
      } else {
        stop(cat(paste("Cannot find ", signaling[x], ".", "Please input a correct name!"),'\n'))
      }
    }
  )
  if (pair.only) {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[c(4*i-3, 4*i-2, 4*i-1, 4*i)]), ncol=4, byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]][1:4]
    rownames(pairLR) <- pairLR[,1]
  } else {
    pairLR0 <- vector("list", length(signaling))
    for (i in 1:length(signaling)) {
      pairLR0[[i]] <- matrix(unlist(pairLR[(i*ncol(pairLR.use)-(ncol(pairLR.use)-1)):(i*ncol(pairLR.use))]), ncol=ncol(pairLR.use), byrow=F)
    }
    pairLR <- do.call(rbind, pairLR0)
    dimnames(pairLR)[[2]] <- dimnames(pairLR.use)[[2]]
    rownames(pairLR) <- pairLR[,1]
  }
  return(as.data.frame(pairLR, stringsAsFactors = FALSE))
}

#' Subset CellChatDB databse by only including interactions of interest
#'
#' @param CellChatDB CellChatDB databse
#' @param search a character vector, which is a subset of c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling"); Setting search = NULL & non_protein = FALSE will return all signaling except for "Non-protein Signaling".
#'
#' When `key` is a vector, the `search` should be a list with the size being `length(key)`, where each element is a character vector.
#' @param key a character vector and each element should be one of the column names of the interaction_input from CellChatDB.
#' @param non_protein whether to use the non-protein signaling for CellChat analysis. By default, non_protein = FALSE because most of non-protein signaling are the special synaptic signaling interactions that can only be used when inferring neuron-neuron communication.
#'
#' @return
#' @export
#'
subsetDB <- function(CellChatDB, search = c(), key = "annotation", non_protein = FALSE) {
  interaction_input <- CellChatDB$interaction
  if (is.null(search) & non_protein == FALSE & any(key == "annotation")) {
    search <- c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")
  } else if (is.null(search) & non_protein == TRUE & any(key == "annotation")) {
    search <- c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling")
  }

  if ("Non-protein Signaling" %in% unlist(search)) {
    non_protein = TRUE
    message("The non-protein signaling is now included for CellChat analysis, which is usually used for neuron-neuron and metabolic communication!")
  }
  if (non_protein == FALSE) {
    interaction_input <- subset(interaction_input, annotation != "Non-protein Signaling")
  }
  if (all(key %in% colnames(interaction_input)) == FALSE) {
    stop("Each element of the `key` should be one of the column names of the interaction_input from CellChatDB")
  }
  if (length(key) == 1) {
    interaction_input <- interaction_input[interaction_input[[key]] %in% search, ]
  } else {
    if (!is.list(search)) {
      stop("When `key` is a vector, the `search` should be a list. ")
    }
    idx.use <- TRUE
    for (i in 1:length(key)) {
      idx.use <- idx.use & (interaction_input[[key[i]]] %in% search[[i]])
    }
    interaction_input <- interaction_input[idx.use, , drop = FALSE]
  }

  CellChatDB$interaction <- interaction_input
  return(CellChatDB)
}



#' Extract the genes involved in CellChatDB
#'
#' @param CellChatDB CellChatDB databse used in the analysis
#'
#' @return
#' @export
#' @importFrom dplyr select
#'
extractGene <- function(CellChatDB) {
  interaction_input <- CellChatDB$interaction
  complex_input <- CellChatDB$complex
  cofactor_input <- CellChatDB$cofactor
  geneIfo <- CellChatDB$geneInfo
  # check whether all gene names in complex_input and cofactor_input are official gene symbol in geneIfo
  checkGeneSymbol(geneSet = unlist(complex_input), geneIfo)
  checkGeneSymbol(geneSet = unlist(cofactor_input), geneIfo)

  geneL <- unique(interaction_input$ligand)
  geneR <- unique(interaction_input$receptor)
  geneLR <- c(geneL, geneR)
  checkGeneSymbol(geneSet = geneLR[geneLR %in% rownames(complex_input) == "FALSE"], geneIfo)

  geneL <- extractGeneSubset(geneL, complex_input, geneIfo)
  geneR <- extractGeneSubset(geneR, complex_input, geneIfo)
  geneLR <- c(geneL, geneR)

  cofactor <- c(interaction_input$agonist, interaction_input$antagonist, interaction_input$co_A_receptor, interaction_input$co_I_receptor)
  cofactor <- unique(cofactor[cofactor != ""])
  cofactorsubunits <- select(cofactor_input[match(cofactor, rownames(cofactor_input), nomatch=0),], starts_with("cofactor"))
  cofactorsubunitsV <- unlist(cofactorsubunits)
  geneCofactor <- unique(cofactorsubunitsV[cofactorsubunitsV != ""])

  gene.use <- unique(c(geneLR, geneCofactor))
  return(gene.use)

}


#' Extract the gene name
#'
#' @param geneSet gene set
#' @param complex_input complex in CellChatDB databse
#' @param geneIfo official gene symbol
#'
#' @return
#' @importFrom dplyr select starts_with
#' @export
extractGeneSubset <- function(geneSet, complex_input, geneIfo) {
  complex <- geneSet[which(geneSet %in% geneIfo$Symbol == "FALSE")]
  geneSet <- intersect(geneSet, geneIfo$Symbol)
  complexsubunits <- dplyr::select(complex_input[match(complex, rownames(complex_input), nomatch=0),], starts_with("subunit"))
  complex <- intersect(complex, rownames(complexsubunits))
  complexsubunitsV <- unlist(complexsubunits)
  complexsubunitsV <- unique(complexsubunitsV[complexsubunitsV != ""])
  geneSet <- unique(c(geneSet, complexsubunitsV))
  return(geneSet)
}


#' Extract the signaling gene names from ligand-receptor pairs
#'
#' @param pairLR data frame must contain columns named `ligand` and `receptor`
#' @param object a CellChat object
#' @param complex_input complex in CellChatDB databse
#' @param geneInfo official gene symbol
#' @param combined whether combining the ligand genes and receptor genes
#'
#' @return
#' @export
extractGeneSubsetFromPair <- function(pairLR, object = NULL, complex_input = NULL, geneInfo = NULL, combined = TRUE) {
  if (!all(c("ligand", "receptor") %in% colnames(pairLR))) {
    stop("The input data frame must contain columns named `ligand` and `receptor`")
  }
  if (is.null(object)) {
    if (is.null(complex_input) | is.null(geneInfo)) {
      stop("Either `object` or `complex_input` and `geneInfo` should be provided!")
    } else {
      complex <- complex_input
    }
  } else {
    complex <- object@DB$complex
    geneInfo <- object@DB$geneInfo
  }
  geneL <- unique(pairLR$ligand)
  geneR <- unique(pairLR$receptor)
  geneL <- extractGeneSubset(geneL, complex, geneInfo)
  geneR <- extractGeneSubset(geneR, complex, geneInfo)
  geneLR <- c(geneL, geneR)
  if (combined) {
    return(geneLR)
  } else {
    return(list(geneL = geneL, geneR = geneR))
  }
}



#' check the official Gene Symbol
#'
#' @param geneSet gene set to check
#' @param geneIfo official Gene Symbol
#' @return
#' @export
#'
checkGeneSymbol <- function(geneSet, geneIfo) {
  geneSet <- unique(geneSet[geneSet != ""])
  genes_notOfficial <- geneSet[geneSet %in% geneIfo$Symbol == "FALSE"]
  if (length(genes_notOfficial) > 0) {
    cat("Issue identified!! Please check the official Gene Symbol of the following genes: ", "\n", genes_notOfficial, "\n")
  }
  return(FALSE)
}

#' Extract L-R pairs associated with a given gene set
#'
#' @param geneSet a vector of genes
#' @param db one of the CellChatDB databases (e.g., CellChatDB.human, CellChatDB.mouse...)
#' @export
#'
extractLRfromGenes <- function(geneSet, db) {
  interaction_input <- db$interaction
  complex_input <- db$complex
  geneIfo <- db$geneInfo
  geneSet1 <- intersect(geneSet, geneIfo$Symbol)
  idx1 <- which(interaction_input$ligand %in% geneSet1)
  idx2 <- which(interaction_input$receptor %in% geneSet1)
  idx <- unique(c(idx1, idx2)); idx <- setdiff(idx,0)
  LR.use <- interaction_input[idx,,drop = FALSE]
  genes.use <- extractGeneSubsetFromPair(LR.use, complex_input = complex_input, geneInfo = geneIfo)
  return(list(LR.use = LR.use, genes.use=genes.use))
}


#' Update CellChatDB by integrating new L-R pairs from other resources or adding more information
#'
#' @param db a data frame of the customized ligand-receptor database with at least two columns named as `ligand` and `receptor`. We highly suggest users to provide a column of pathway information named `pathway_name` associated with each L-R pair.
#' Other optional columns include `interaction_name` and `interaction_name_2`. The default columns of CellChatDB can be checked via `colnames(CellChatDB.human$interaction)`.
#' @param gene_info a data frame with at least one column named as `Symbol`. "When setting gene_info = NULL, the input `species_target` should be provided: either `human` or `mouse`.
#' @param other_info a list consisting of other information including a dataframe named as `complex` and a dataframe named as `cofactor`. This additional information is not necessary. If other_info is provided, the `complex` and `cofactor` are dataframes with defined rownames.
#' @param gene_info_columnNew a data frame with at least two columns named as `Symbol` and `AntibodyName`, which will add a new column named `AntibodyName` into `db$geneInfo`.
#' @param trim.pathway whether to delete the interactions with missing pathway names when the column `pathway_name` is provided in `db`.
#' @param merged whether merging the input database with the existing CellChatDB. setting merged = TRUE, the input `species_target` should be provided: either `human` or `mouse`.
#' @param species_target the target species for output: either `human` or `mouse`.
#' @return a list consisting of the customized L-R database for further CellChat analysis
#' @export
#'
#' @examples
#'\dontrun{
#' # integrating new L-R pairs from other resources or utilizing a custom database `db.user`
#' db.new <- updateCellChatDB(db = db.user, gene_info = gene_info)
#' db.new <- updateCellChatDB(db = db.user, gene_info = NULL, species_target = "human")
#' # Alternatively, users can integrate the customized L-R pairs into the built-in CellChatDB
#' db.new <- updateCellChatDB(db = db.user, merged = TRUE, species_target = "human")
#' # Add new columns (e.g., AntibodyName) into gene_info
#' db.new.human <- updateCellChatDB(db = CellChatDB.human$interaction, gene_info = CellChatDB.human$geneInfo, other_info=list(complex = CellChatDB.human$complex, cofactor = CellChatDB.human$cofactor),gene_info_columnNew = gene_info_columnNew)
#'
#' # Users can now use this new database in CellChat analysis
#' cellchat@DB <- db.new
#'}
updateCellChatDB <- function(db, gene_info = NULL, other_info = NULL, gene_info_columnNew = NULL, trim.pathway = FALSE, merged = FALSE, species_target = NULL) {
  db <- dplyr::mutate(db, across(everything(), as.character))
  if (all(c("ligand","receptor") %in% colnames(db)) == FALSE) {
    stop("The input `db` must contain at least two columns named as ligand,receptor")
  }
  if (all(c("pathway_name") %in% colnames(db)) == FALSE) {
    warning("The pathway_name associated with each L-R pair is not provided in `db`. We suggest to provide this information so that the versatile functionalities of CellChat can be fully used! \n")
    db$pathway_name <- rep("", nrow(db))
  } else {
    pathway.missing <- which(db$pathway_name == "")
    if (length(pathway.missing) > 0) {
      if (trim.pathway) {
        cat(paste0("The pathway names of ", length(pathway.missing) ," interactions are missing and the corresponding interactions are now deleted. \n"))
        db <- db[-pathway.missing, , drop = FALSE]
      } else {
        warning(paste0("The pathway names of ", length(pathway.missing) ," interactions are missing and it may cause error in the downstream analysis. Setting `trim.pathway = TRUE` to avoid such possible errors. \n"))
      }
    }
  }
  if (all(c("interaction_name") %in% colnames(db)) == FALSE) {
    db$interaction_name <- paste0(toupper(db$ligand), "_", toupper(db$receptor))
  }
  if (all(c("interaction_name_2") %in% colnames(db)) == FALSE) {
    db$interaction_name_2 <- paste0(db$ligand, " - ", db$receptor)
  }
  if ("agonist" %in% colnames(db) == FALSE) {
    db$agonist <- rep("", nrow(db))
  }
  if ("antagonist" %in% colnames(db) == FALSE) {
    db$antagonist <- rep("", nrow(db))
  }
  if ("co_A_receptor" %in% colnames(db) == FALSE) {
    db$co_A_receptor <- rep("", nrow(db))
  }
  if ("co_I_receptor" %in% colnames(db) == FALSE) {
    db$co_I_receptor <- rep("", nrow(db))
  }
  ## construct database
  idx.remove <- duplicated(db$interaction_name)
  if (sum(idx.remove) > 0) {
    warning(paste0(sum(idx.remove), " duplicated interaction_names are identified and the corresponding interactions are now deleted. \n"))
    db <- db[-which(idx.remove), ]
  }

  # build the interaction file
  interaction_input <- db
  rownames(interaction_input) <- interaction_input$interaction_name
  cols.default <- c("interaction_name","pathway_name","ligand","receptor","agonist","antagonist","co_A_receptor","co_I_receptor","annotation","interaction_name_2")
  cols.common <- intersect(cols.default,colnames(interaction_input))
  cols.specific <- setdiff(colnames(interaction_input), cols.default)
  interaction_input <- dplyr::select(interaction_input, c(cols.common, cols.specific))

  # build the complex file
  if (!is.null(other_info)) {
    if ("complex" %in% names(other_info) == TRUE) {
      complex_input <- other_info$complex
      if (all(colnames(complex_input) %in% paste0("subunit_", seq_len(100))) == FALSE) {
        stop("The colnames of the input `other_info$complex` should be `subunit_1`,`subunit_2`,...")
      }
    } else {
      complex_input <- data.frame()
    }
    # build the cofactor file
    if ("cofactor" %in% names(other_info) == TRUE) {
      cofactor_input <- other_info$cofactor
      if (all(colnames(cofactor_input) %in% paste0("cofactor", seq_len(100))) == FALSE) {
        stop("The colnames of the input `other_info$cofactor` should be `cofactor1`,`cofactor2`,...")
      }
    } else {
      cofactor_input <- data.frame()
    }
  } else {
    complex_input <- data.frame()
    cofactor_input <- data.frame()
  }

  # build the geneInfo file
  if (!is.null(gene_info)) {
    if ("Symbol" %in% colnames(gene_info) == FALSE) {
      stop("The input `gene_info` must contain at least one column named as `Symbol`")
    }
  } else {
    if (is.null(species_target)) {
      stop("When setting gene_info = NULL, the input `species_target` should be provided: either `human` or `mouse`. ")
    }
    if (species_target == "human") {
      gene_info <- CellChatDB.human$geneInfo
    } else if (species_target == "mouse") {
      gene_info <- CellChatDB.mouse$geneInfo
    }
  }
  geneInfo_input <- gene_info

  if (merged == TRUE) {
    if (is.null(species_target)) {
      stop("When setting merged = TRUE, the input `species_target` should be provided: either `human` or `mouse`. ")
    }
    if (species_target == "human") {
      db.cellchat <- CellChatDB.human
      cat("Starting to merge the input database with CellChatDB.human... \n")
    } else if (species_target == "mouse") {
      db.cellchat <- CellChatDB.mouse
      cat("Starting to merge the input database with CellChatDB.mouse... \n")
    }

    # build the interaction file
    interaction_input.cellchat <- db.cellchat$interaction
    interaction_input.cellchat$source.merged <- "CellChatDB"
    interaction_input$source.merged <- "User"
    cols.common <- intersect(colnames(interaction_input), colnames(interaction_input.cellchat))
    interaction_input <- interaction_input[, cols.common]
    interaction_input.cellchat <- interaction_input.cellchat[, cols.common]
    interaction_input.merged <- rbind(interaction_input.cellchat, interaction_input)
    idx.remove <- duplicated(interaction_input.merged$interaction_name)
    if (sum(idx.remove) > 0) {
      interaction_input.merged <- interaction_input.merged[-which(idx.remove), ]
    }

    # build the complex file
    complex_input.cellchat <- db.cellchat$complex
    num.subunit <- max(ncol(complex_input), ncol(complex_input.cellchat))
    if (ncol(complex_input) < num.subunit) {
      temp <- data.frame(rep("", nrow(complex_input)))
      complex_input <- cbind(complex_input, as.data.frame(do.call(cbind, rep(temp, num.subunit-ncol(complex_input)))))
      colnames(complex_input) <- paste0("subunit_", seq_len(num.subunit))
    }
    if (ncol(complex_input.cellchat) < num.subunit) {
      temp <- data.frame(rep("", nrow(complex_input.cellchat)))
      complex_input.cellchat <- cbind(complex_input.cellchat, as.data.frame(do.call(cbind, rep(temp, num.subunit-ncol(complex_input.cellchat)))))
      colnames(complex_input.cellchat) <- paste0("subunit_", seq_len(num.subunit))
    }
    complex_input.merged <- rbind(complex_input.cellchat, complex_input)
    idx.remove <- duplicated(rownames(complex_input.merged))
    if (sum(idx.remove) > 0) {
      complex_input.merged <- complex_input.merged[-which(idx.remove), ]
    }

    # build the cofactor file
    cofactor_input.cellchat <- db.cellchat$cofactor
    num.subunit <- max(ncol(cofactor_input), ncol(cofactor_input.cellchat))
    if (ncol(cofactor_input) < num.subunit) {
      temp <- data.frame(rep("", nrow(cofactor_input)))
      cofactor_input <- cbind(cofactor_input, as.data.frame(do.call(cbind, rep(temp, num.subunit-ncol(cofactor_input)))))
      colnames(cofactor_input) <- paste0("cofactor", seq_len(num.subunit))
    }
    if (ncol(cofactor_input.cellchat) < num.subunit) {
      temp <- data.frame(rep("", nrow(cofactor_input.cellchat)))
      cofactor_input.cellchat <- cbind(cofactor_input.cellchat, as.data.frame(do.call(cbind, rep(temp, num.subunit-ncol(cofactor_input.cellchat)))))
      colnames(cofactor_input.cellchat) <- paste0("cofactor", seq_len(num.subunit))
    }
    cofactor_input.merged <- rbind(cofactor_input.cellchat, cofactor_input)
    idx.remove <- duplicated(rownames(cofactor_input.merged))
    if (sum(idx.remove) > 0) {
      cofactor_input.merged <- cofactor_input.merged[-which(idx.remove), ]
    }

    interaction_input <- interaction_input.merged
    complex_input <- complex_input.merged
    cofactor_input <- cofactor_input.merged
  }

  if (!is.null(gene_info_columnNew)) {
    checkGeneSymbol(gene_info_columnNew$Symbol, geneInfo_input)
    idx <- match(gene_info_columnNew$Symbol, geneInfo_input$Symbol)
    geneInfo_input$AntibodyName <- NA
    geneInfo_input$AntibodyName[idx[!is.na(idx)]] <- gene_info_columnNew$AntibodyName[!is.na(idx)]
  }
  db.new <- list()
  db.new$interaction <- interaction_input
  db.new$complex <- complex_input
  db.new$cofactor <- cofactor_input
  db.new$geneInfo <- geneInfo_input

  return(db.new)
}
