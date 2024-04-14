#' Generate a Shiny App for interactive exploration of CellChat's outputs
#'
#' @param object CellChat object
#' @param ... Other parameters of `shinyApp` function from shiny R package
#' @return A Shiny app object on the basis of one CellChat object
#' @export
#' @importFrom stringr str_split_1
# #' @importFrom plotly subplot plot_ly ggplotly add_markers highlight highlight_key plotlyOutput layout
# #' @importFrom bsicons bs_icon
#' @import shiny bslib
#'
runCellChatApp <- function(object,...) {
  # ##########################################################################
  # set some global options
  # ##########################################################################
  options(stringsAsFactors = FALSE)

  # ##########################################################################
  # some useful elements for ui.R
  # ##########################################################################
  choices_cell_groups <-levels(object@idents)
  names(choices_cell_groups) <- levels(object@idents)

  choices_pathways <- object@netP$pathways
  names(choices_pathways) <- object@netP$pathways

  # all signaling gene names
  choices_gene_names <- CellChat::extractGene(object@DB)
  # all ligand-receptor pair names
  #choices_pairLR_use <- object@DB$interaction$interaction_name
  if ("LRs" %in% names(object@net)) {
    choices_pairLR_use <- object@net$LRs
  } else {
    thresh = 0.05
    prob <- object@net$prob
    prob[object@net$pval > thresh] <- 0
    LR <- dimnames(prob)[[3]]
    LR.sig <- LR[apply(prob, 3, sum) != 0]
    choices_pairLR_use <- LR.sig
  }


  # Palettes (sequential)
  choices_palettes_sequential <- stringr::str_split_1("Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu, YlOrBr, YlOrRd",", ")
  names(choices_palettes_sequential) <- choices_palettes_sequential
  choices_palettes_diverging <- stringr::str_split_1("BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral",", ")
  names(choices_palettes_diverging) <- choices_palettes_diverging

  # ##########################################################################
  # interactive visualization
  # ##########################################################################

  # interactive Heatmap
  # [Colors (ggplot2)](http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
  plotly_netVisual_heatmap <- function(obj_heatmap,palette.heatmap,direction.heatmap=1) {
    gg_heatmap <- obj_heatmap@matrix %>%
      as.data.frame() %>%
      mutate(row = rownames(.)) %>%
      tidyr::pivot_longer(
        data = .,
        cols = colnames(.)[-length(colnames(.))],
        names_to = "column",
        values_to = "value"
      ) %>%
      ggplot() +
      geom_tile(aes(row, column, fill = value),
                width = 0.95,
                height = 0.95) +
      # guides(fill=guide_legend(title=obj_heatmap@row_title))+
      labs(title = '',
           x = '',
           y = obj_heatmap@row_title,
           # I can't set the direction of the legend title, I thick it's a bug
           # fill = obj_heatmap@column_title,
      ) +
      scale_fill_distiller(
        palette = palette.heatmap,
        na.value = 'white',
        direction = direction.heatmap,
      ) +
      theme_minimal()+
      theme(axis.title.y = element_text(size = 14))

    # ggplot transpose the matrix, so we need use colSums to calc the 'rowSums'
    # of the matrix
    gg_right <- obj_heatmap@matrix %>%
      colSums(abs(.)) %>%
      tibble(row_sum = ., sources_name = names(.)) %>%
      ggplot() +
      geom_bar(aes(x = sources_name, y = row_sum, fill = sources_name),
               stat = 'identity') +
      labs(title = '',
           x = '',
           y = '',) +
      guides(fill = FALSE) +
      scale_fill_brewer(palette = "Set1", direction = 1) +
      theme_minimal() +
      coord_flip()

    gg_top <- obj_heatmap@matrix %>%
      rowSums(abs(.)) %>%
      tibble(col_sum = ., sources_name = names(.)) %>%
      ggplot() +
      # use fill to set the columns' colors
      geom_bar(aes(x = sources_name, y = col_sum, fill = sources_name),
               stat = 'identity') +
      labs(title = obj_heatmap@column_title,
           x = '',
           y = '',) +
      guides(fill = FALSE)+
      scale_fill_brewer(palette = "Set1", direction = 1) +
      # theme() function should be used behind the theme_*()
      theme_minimal()+
      theme(plot.title = element_text(hjust = 0.5,size = 14))

    return(plotly::subplot(
      gg_top,
      plotly::plotly_empty(),
      gg_heatmap,
      gg_right,
      nrows = 2,
      heights = c(0.2, 0.8),
      widths = c(0.8, 0.2),
      margin = 0,
      shareX = TRUE,
      shareY = TRUE,
      titleX = TRUE,
      titleY = TRUE
    )
    )
  }

  # interactive DimPlot
  plotly_DimPlot <- function (object,
                              color.use = NULL,
                              group.by = NULL,
                              sample.use = NULL,
                              reduction = NULL,
                              sources.use = NULL,
                              targets.use = NULL,
                              idents.use = NULL,
                              alpha = 1,
                              title.name = NULL,
                              point.size = 1)
  {
    if (is.null(group.by)) {
      labels <- object@idents
    } else {
      labels = object@meta[, group.by]
      labels <- factor(labels)
    }
    if (length(names(object@dr)) == 0) {
      stop("Please check `addReduction` to add a new reduced space into `object@dr`. \n")
    }
    if (!is.null(reduction)) {
      coords <- object@dr[[reduction]]
    } else {
      if ("umap" %in% names(object@dr)) {
        coords <- object@dr$umap
      } else if ("tsne" %in% names(object@dr)){
        coords <- object@dr$tsne
      } else {
        stop(paste0("The `object@dr` contains the following reduced space: ", toString(names(object@dr)), ". Please specify the dimensionality reduction to use. \n"))
      }
    }
    coordinates <- as.data.frame(coords)
    samples <- object@meta$samples
    if (ncol(coordinates) >= 2) {
      coordinates <- coordinates[, c(1,2)]
      colnames(coordinates) <- c("x_cent","y_cent")
      if (length(unique(samples)) > 1) {
        if (is.null(sample.use)) {
          stop("`sample.use` should be provided for visualizing signaling on each individual sample.")
        } else if (sample.use %in% unique(samples)) {
          coordinates = coordinates[samples == sample.use, ]
          labels = labels[samples == sample.use]
        } else {
          stop("Please check the input `sample.use`, which should be the element in `meta$samples`.")
        }
      }
      # temp_coordinates = coordinates
      # coordinates[,1] = temp_coordinates[,2]
      # coordinates[,2] = temp_coordinates[,1]
    } else {
      stop("Please check the input 'object@dr' and make sure it has at least two columns.")
    }



    cells.level <- levels(labels)
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      cell.use <- !(labels %in% idents.use)
      labels[cell.use] <- NA
      cells.level <- cells.level[cells.level %in% idents.use]
      labels <- factor(labels, levels = cells.level)
    }
    if (is.null(sources.use) & is.null(targets.use)) {
      if (is.null(color.use)) {
        color.use <- scPalette(nlevels(labels))
      }
    }
    else {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      group <- rep("Others", length(labels))
      group[(labels %in% sources.use)] <- sources.use
      group[(labels %in% targets.use)] <- targets.use
      group = factor(group, levels = c(sources.use, targets.use,
                                       "Others"))
      if (is.null(color.use)) {
        color.use.all <- scPalette(nlevels(labels))
        color.use <- color.use.all[match(c(sources.use,
                                           targets.use), levels(labels))]
        color.use[nlevels(group)] <- "grey90"
      }
      labels <- group
    }
    # print(color.use)->return a color vector
    coordinates$cell_labels <- labels

    py <- plotly::highlight_key(coordinates,~cell_labels) %>%
      plotly::plot_ly(x = ~x_cent, y = ~y_cent,marker = list(size = point.size)) %>%
      plotly::add_markers(color=~cell_labels,alpha=alpha,colors=color.use) %>%
      plotly::layout(
        title = title.name,
        yaxis = list(
          title = "",
          #autorange = "reversed",
          showgrid = FALSE,
          ticks = "",
          # ticktext = "",
          tickvals = "",
          zeroline = FALSE,
          showline = FALSE
        ),
        xaxis = list(
          title = "",
          showgrid = FALSE,
          ticks = "",
          # ticktext = "",
          tickvals = "",
          zeroline = FALSE,
          showline = FALSE
        )
      ) %>%
      plotly::highlight(on = "plotly_click",
                        off = "plotly_relayout")

    return(py)
  }

  # interactive FeaturePlot
  # https://plotly.com/r/subplots/
  plotly_FeaturePlot <- function (object,
                                  features = NULL,
                                  signaling = NULL,
                                  pairLR.use = NULL,
                                  sample.use = NULL,
                                  reduction = NULL,
                                  enriched.only = TRUE,
                                  thresh = 0.05,
                                  do.group = TRUE,
                                  color.heatmap = "Reds",
                                  n.colors = 8,
                                  direction = -1,
                                  do.binary = FALSE,
                                  cutoff = NULL,
                                  color.use = NULL,
                                  alpha = 1,
                                  point.size = 0.8,
                                  legend.size = 3,
                                  legend.text.size = 8,
                                  shape.by = 16,
                                  plot_nrows = 1,
                                  show.legend = TRUE,
                                  show.legend.combined = FALSE){
    if (!is.null(reduction)) {
      coords <- object@dr[[reduction]]
    } else {
      if ("umap" %in% names(object@dr)) {
        coords <- object@dr$umap
      } else if ("tsne" %in% names(object@dr)){
        coords <- object@dr$tsne
      } else {
        stop("Please make sure `object@dr` contains a low-dimensional space of the data and specify the dimensionality reduction to use.")
      }
    }

    samples <- object@meta$samples
    cell_labels <- object@idents
    data <- as.matrix(object@data)
    meta <- object@meta
    coords <- as.data.frame(coords)
    if (ncol(coords) >= 2) {
      coords <- coords[, c(1,2)]
      colnames(coords) <- c("x_cent","y_cent")
      if (length(unique(samples)) > 1) {
        if (is.null(sample.use)) {
          stop("`sample.use` should be provided for visualizing signaling on each individual sample.")
        } else if (sample.use %in% unique(samples)) {
          coords = coords[samples == sample.use, ]
          meta = meta[samples == sample.use, ]
          data = data[, samples == sample.use]
        } else {
          stop("Please check the input `sample.use`, which should be the element in `meta$samples`.")
        }
      }
    } else {
      stop("Please check the input 'object@dr' and make sure it has at least two columns.")
    }

    # add idents info

    if (length(color.heatmap) == 1) {
      colormap <- tryCatch({
        RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
      }, error = function(e) {
        (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
      })
      if (direction == -1) {
        colormap <- rev(colormap)
      }
      colormap <- colorRampPalette(colormap)(99)
      colormap[1] <- "#E5E5E5"
    }
    else {
      colormap <- color.heatmap
    }
    if (is.null(features) &
        is.null(signaling) & is.null(pairLR.use)) {
      stop("Please input either features, signaling or pairLR.use.")
    }
    if (!is.null(features) & !is.null(signaling)) {
      stop("Please don't input features or signaling simultaneously.")
    }
    if (!is.null(features) & !is.null(pairLR.use)) {
      stop("Please don't input features or pairLR.use simultaneously.")
    }
    if (!is.null(signaling) & !is.null(pairLR.use)) {
      stop("Please don't input signaling or pairLR.use simultaneously.")
    }
    df <- data.frame(x = coords[, 1], y = coords[, 2],
                     cell_labels = cell_labels)


    if (!do.binary) {
      if (!is.null(signaling)) {
        res <- extractEnrichedLR(
          object,
          signaling = signaling,
          geneLR.return = TRUE,
          enriched.only = enriched.only,
          thresh = thresh
        )
        feature.use <- res$geneLR
      }
      else if (!is.null(pairLR.use)) {
        if (is.character(pairLR.use)) {
          pairLR.use <- data.frame(interaction_name = pairLR.use)
        }
        if (enriched.only) {
          if (do.group) {
            object@net$prob[object@net$pval > thresh] <- 0
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob)[[3]]]
            prob <- object@net$prob[, , pairLR.use.name,
                                    drop = FALSE]
            prob.sum <- apply(prob > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          else {
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob.cell)[[3]]]
            prob.cell <- object@net$prob.cell[, , pairLR.use.name,
                                              drop = FALSE]
            prob.sum <- apply(prob.cell > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          if (length(pairLR.use$interaction_name) == 0) {
            stop(
              paste0(
                "There is no significant communication related with the input `pairLR.use`. Set `enriched.only = FALSE` to show non-significant signaling."
              )
            )
          }
        }
        LR.pair <- object@LR$LRsig[pairLR.use$interaction_name,
                                   c("ligand", "receptor")]
        geneL <- unique(LR.pair$ligand)
        geneR <- unique(LR.pair$receptor)
        geneL <- extractGeneSubset(geneL, object@DB$complex,
                                   object@DB$geneInfo)
        geneR <- extractGeneSubset(geneR, object@DB$complex,
                                   object@DB$geneInfo)
        feature.use <- c(geneL, geneR)
      }
      else {
        feature.use <- features
      }
      if (length(intersect(feature.use, rownames(data))) >
          0) {
        feature.use <- feature.use[feature.use %in% rownames(data)]
        data.use <- data[feature.use, , drop = FALSE]
      }
      else if (length(intersect(feature.use, colnames(meta))) >
               0) {
        feature.use <- feature.use[feature.use %in% colnames(meta)]
        data.use <- t(meta[, feature.use, drop = FALSE])
      }
      else {
        stop("Please check your input! ")
      }
      if (!is.null(cutoff)) {
        cat("Applying a cutoff of ", cutoff, "to the values...",
            "\n")
        data.use[data.use <= cutoff] <- 0
      }


      numFeature = length(feature.use)
      gg <- vector("list", numFeature)

      plot_ncols <- ceiling(numFeature / plot_nrows)
      plot_width <- 1 / plot_ncols
      plot_height <- 1 / plot_nrows
      annotations_pos <- vector("list", 0)
      for (i in 1:plot_nrows) {
        for (j in 1:plot_ncols) {
          x <- (j - 0.5) * plot_width
          y <- 1-(i - 0.95) * plot_height
          annotations_pos[[length(annotations_pos)+1]] <- list(
            x=x,
            y=y
          )
        }

      }
      annotations <- vector("list", numFeature)

      for (i in seq_len(numFeature)) {
        feature.name <- feature.use[i]
        df$feature.data <- data.use[i,]
        g <-
          ggplot(data = df, aes(x, y)) + geom_point(
            aes(colour = feature.data,cell_labels=cell_labels),
            alpha = alpha,
            size = point.size,
            shape = shape.by
          ) +
          scale_colour_gradientn(
            colours = colormap,
            guide = guide_colorbar(
              title = NULL,
              ticks = T,
              label = T,
              barwidth = 0.5
            ),
            na.value = "grey90"
          ) +
          theme(legend.position = "right") + theme(
            legend.title = element_blank(),
            legend.text = element_text(size = legend.text.size),
            legend.key.size = unit(0.15, "inches")
          ) + ggtitle(feature.name) +
          theme(plot.title = element_text(
            hjust = 0.5,
            vjust = 0,
            size = 10
          )) + theme(
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
          ) +
          xlab(NULL) + ylab(NULL) + theme(legend.key = element_blank())
        # g <- g + coord_fixed() +
        #   scale_y_reverse()

        gg[[i]] <- g %>% plotly::ggplotly(height = 400)

        annotations[[i]] <- list(
          x=annotations_pos[[i]]$x,
          y=annotations_pos[[i]]$y,
          text = feature.name,
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "top",
          showarrow = FALSE,
          font = list(size = 16))
      }
      # print(annotations_pos)
      if (numFeature > 1) {
        gg <- plotly::subplot(gg, nrows = plot_nrows,margin = 0.02,shareX = FALSE,shareY = FALSE) %>%
          plotly::layout(title = '',annotations = annotations)
      }
      else {
        gg <- plotly::ggplotly(gg[[1]])
      }
    }
    else {
      # do.binary
      set_individual_legend <- function(plt) {
        # plt is a plotly plot obj
        plt_build <- plotly::plotly_build(plt)

        # get the num of traces
        len_legend <- length(plt_build$x$data)

        for (i in 1:len_legend) {
          # set legendgroup
          plt_build$x$data[[i]]$legendgroup <- feature.name
          # set legendtitle
          plt_build$x$data[[i]]$legendgrouptitle <- list(text=feature.name,font=list(size=12))
        }
        return(plt_build)
      }
      if (is.null(color.use)) {
        color.use <- ggPalette(4)
        color.use[4] <- "grey90"
      }
      color.use1 = color.use
      if (!is.null(signaling)) {
        res <- extractEnrichedLR(
          object,
          signaling = signaling,
          enriched.only = enriched.only,
          thresh = thresh
        )
        LR.pair <- object@LR$LRsig[res$interaction_name,
                                   c("ligand", "receptor")]
      }
      else if (!is.null(pairLR.use)) {
        if (is.character(pairLR.use)) {
          pairLR.use <- data.frame(interaction_name = pairLR.use)
        }
        if (enriched.only) {
          if (do.group) {
            object@net$prob[object@net$pval > thresh] <- 0
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob)[[3]]]
            prob <- object@net$prob[, , pairLR.use.name,
                                    drop = FALSE]
            prob.sum <- apply(prob > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          else {
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob.cell)[[3]]]
            prob.cell <- object@net$prob.cell[, , pairLR.use.name,
                                              drop = FALSE]
            prob.sum <- apply(prob.cell > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          if (length(pairLR.use$interaction_name) == 0) {
            stop(
              paste0(
                "There is no significant communication related with the input `pairLR.use`. Set `enriched.only = FALSE` to show non-significant signaling."
              )
            )
          }
        }
        LR.pair <- object@LR$LRsig[pairLR.use$interaction_name,
                                   c("ligand", "receptor")]
      }
      else {
        stop("Please input either `pairLR.use` or `signaling` for `binary` mode!")
      }
      geneL <- as.character(LR.pair$ligand)
      geneR <- as.character(LR.pair$receptor)
      complex_input <- object@DB$complex
      dataL <- computeExpr_LR(geneL, data, complex_input)
      dataR <- computeExpr_LR(geneR, data, complex_input)
      rownames(dataL) <- geneL
      rownames(dataR) <- geneR

      feature.use <- rownames(LR.pair)
      numFeature = nrow(LR.pair)

      if (is.null(cutoff)) {
        stop("A `cutoff` must be provided when plotting expression in binary mode! ")
      }
      gg <- vector("list", numFeature)

      # set subplot title pos
      plot_ncols <- ceiling(numFeature / plot_nrows)
      plot_width <- 1 / plot_ncols
      plot_height <- 1 / plot_nrows
      annotations_pos <- vector("list", 0)
      for (i in 1:plot_nrows) {
        for (j in 1:plot_ncols) {
          x <- (j - 0.5) * plot_width
          y <- 1-(i - 1) * plot_height
          annotations_pos[[length(annotations_pos)+1]] <- list(
            x=x,
            y=y
          )
        }

      }
      annotations <- vector("list", numFeature)

      for (i in seq_len(numFeature)) {
        feature.name <- feature.use[i]
        idx1 = dataL[i,] > cutoff
        idx2 = dataR[i,] > cutoff
        idx3 = idx1 & idx2
        group = rep("None", ncol(dataL))
        group[idx1] = geneL[i]
        group[idx2] = geneR[i]
        group[idx3] = "Both"
        group = factor(group, levels = c(geneL[i], geneR[i],
                                         "Both", "None"))
        color.use <- color.use1
        names(color.use) <- c(geneL[i], geneR[i], "Both",
                              "None")
        if (length(setdiff(levels(group), unique(group))) >
            0) {
          color.use <- color.use[names(color.use) %in% unique(group)]
          group = droplevels(group, exclude = setdiff(levels(group),
                                                      unique(group)))
        }
        df$feature.data <- group
        g <-
          ggplot(data = df, aes(x, y)) + geom_point(
            aes(colour = feature.data,cell_labels=cell_labels),
            alpha = alpha,
            size = point.size,
            shape = shape.by
          ) +
          scale_color_manual(values = color.use, na.value = "grey90") +
          theme(legend.position = "right") + theme(
            legend.title = element_blank(),
            legend.text = element_text(size = legend.text.size),
            legend.key.size = unit(0.15, "inches")
          ) + guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(feature.name) + theme(plot.title = element_text(
            hjust = 0.5,
            vjust = 0,
            size = 10
          )) + theme(
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
          ) +
          xlab(NULL) + ylab(NULL) + theme(legend.key = element_blank())
        # g <- g + coord_fixed() +
        #   scale_y_reverse()

        # cat(feature.name)
        gg[[i]] <- g %>% plotly::ggplotly(
          type = 'scatter',
          mode='markers+text',
        ) %>% set_individual_legend()

        annotations[[i]] <- list(
          x=annotations_pos[[i]]$x,
          y=annotations_pos[[i]]$y,
          text = feature.name,
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "top",
          showarrow = FALSE,
          font = list(size = 16))
      }

      if (numFeature > 1) {
        gg <- plotly::subplot(gg, nrows = plot_nrows,margin = 0.02,shareX = FALSE,shareY = FALSE) %>%
          plotly::layout(title = '',
                         annotations = annotations,
                         legend = list(tracegroupgap = 10,title=list(text=''))
          )
      }
      else {
        gg <- plotly::ggplotly(gg[[1]],
                               type = 'scatter',
                               mode = 'markers') %>%
          plotly::layout(legend = list(title = list(text = '')))
      }
    }
    return(gg)
  }

  # interactive spatialDimPlot
  plotly_spatialDimPlot <- function (object,
                              color.use = NULL,
                              group.by = NULL,
                              sample.use = NULL,
                              sources.use = NULL,
                              targets.use = NULL,
                              idents.use = NULL,
                              alpha = 1,
                              title.name = NULL,
                              point.size = 1)
  {
    if (is.null(group.by)) {
      labels <- object@idents
    } else {
      labels = object@meta[, group.by]
      labels <- factor(labels)
    }

    coordinates <- as.data.frame(object@images$coordinates)
    samples <- object@meta$samples
    if (ncol(coordinates) == 2) {
      colnames(coordinates) <- c("x_cent","y_cent")
      if (length(unique(samples)) > 1) {
        if (is.null(sample.use)) {
          stop("`sample.use` should be provided for visualizing signaling on each individual sample.")
        } else if (sample.use %in% unique(samples)) {
          coordinates = coordinates[samples == sample.use, ]
          labels = labels[samples == sample.use]
        } else {
          stop("Please check the input `sample.use`, which should be the element in `meta$samples`.")
        }
      }
      temp_coordinates = coordinates
      coordinates[,1] = temp_coordinates[,2]
      coordinates[,2] = temp_coordinates[,1]
    } else {
      stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
    }



    cells.level <- levels(labels)
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      cell.use <- !(labels %in% idents.use)
      labels[cell.use] <- NA
      cells.level <- cells.level[cells.level %in% idents.use]
      labels <- factor(labels, levels = cells.level)
    }
    if (is.null(sources.use) & is.null(targets.use)) {
      if (is.null(color.use)) {
        color.use <- scPalette(nlevels(labels))
      }
    }
    else {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      group <- rep("Others", length(labels))
      group[(labels %in% sources.use)] <- sources.use
      group[(labels %in% targets.use)] <- targets.use
      group = factor(group, levels = c(sources.use, targets.use,
                                       "Others"))
      if (is.null(color.use)) {
        color.use.all <- scPalette(nlevels(labels))
        color.use <- color.use.all[match(c(sources.use,
                                           targets.use), levels(labels))]
        color.use[nlevels(group)] <- "grey90"
      }
      labels <- group
    }
    # print(color.use)->return a color vector
    coordinates$cell_labels <- labels

    py <- plotly::highlight_key(coordinates,~cell_labels) %>%
      plotly::plot_ly(x = ~x_cent, y = ~y_cent,marker = list(size = point.size)) %>%
      plotly::add_markers(color=~cell_labels,alpha=alpha,colors=color.use) %>%
      plotly::layout(
        title = title.name,
        yaxis = list(
          autorange = "reversed",
          title = "",
          showgrid = FALSE,
          ticks = "",
          # ticktext = "",
          tickvals = "",
          showline = FALSE
        ),
        xaxis = list(
          title = "",
          showgrid = FALSE,
          ticks = "",
          # ticktext = "",
          tickvals = "",
          showline = FALSE
        )
      ) %>%
      plotly::highlight(on = "plotly_click",
                off = "plotly_relayout")

    return(py)
  }

  # interactive spatialFeaturePlot
  # https://plotly.com/r/subplots/
  plotly_spatialFeaturePlot <- function (object,
                                  features = NULL,
                                  signaling = NULL,
                                  pairLR.use = NULL,
                                  sample.use = NULL,
                                  enriched.only = TRUE,
                                  thresh = 0.05,
                                  do.group = TRUE,
                                  color.heatmap = "Reds",
                                  n.colors = 8,
                                  direction = -1,
                                  do.binary = FALSE,
                                  cutoff = NULL,
                                  color.use = NULL,
                                  alpha = 1,
                                  point.size = 0.8,
                                  legend.size = 3,
                                  legend.text.size = 8,
                                  shape.by = 16,
                                  plot_nrows = 1,
                                  show.legend = TRUE,
                                  show.legend.combined = FALSE){
    coords <- as.data.frame(object@images$coordinates)
    samples <- object@meta$samples
    cell_labels <- object@idents
    data <- as.matrix(object@data)
    meta <- object@meta

    if (ncol(coords) == 2) {
      colnames(coords) <- c("x_cent","y_cent")
      if (length(unique(samples)) > 1) {
        if (is.null(sample.use)) {
          stop("`sample.use` should be provided for visualizing signaling on each individual sample.")
        } else if (sample.use %in% unique(samples)) {
          coords = coords[samples == sample.use, ]
          meta = meta[samples == sample.use, ]
          data = data[, samples == sample.use]
        } else {
          stop("Please check the input `sample.use`, which should be the element in `meta$samples`.")
        }
      }
      temp_coords = coords
      coords[,1] = temp_coords[,2]
      coords[,2] = temp_coords[,1]
    } else {
      stop("Please check the input 'coordinates' and make sure it is a two column matrix.")
    }

    # add idents info

    if (length(color.heatmap) == 1) {
      colormap <- tryCatch({
        RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
      }, error = function(e) {
        (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
      })
      if (direction == -1) {
        colormap <- rev(colormap)
      }
      colormap <- colorRampPalette(colormap)(99)
      colormap[1] <- "#E5E5E5"
    }
    else {
      colormap <- color.heatmap
    }
    if (is.null(features) &
        is.null(signaling) & is.null(pairLR.use)) {
      stop("Please input either features, signaling or pairLR.use.")
    }
    if (!is.null(features) & !is.null(signaling)) {
      stop("Please don't input features or signaling simultaneously.")
    }
    if (!is.null(features) & !is.null(pairLR.use)) {
      stop("Please don't input features or pairLR.use simultaneously.")
    }
    if (!is.null(signaling) & !is.null(pairLR.use)) {
      stop("Please don't input signaling or pairLR.use simultaneously.")
    }
    df <- data.frame(x = coords[, 1], y = coords[, 2],
                     cell_labels = cell_labels)


    if (!do.binary) {
      if (!is.null(signaling)) {
        res <- extractEnrichedLR(
          object,
          signaling = signaling,
          geneLR.return = TRUE,
          enriched.only = enriched.only,
          thresh = thresh
        )
        feature.use <- res$geneLR
      }
      else if (!is.null(pairLR.use)) {
        if (is.character(pairLR.use)) {
          pairLR.use <- data.frame(interaction_name = pairLR.use)
        }
        if (enriched.only) {
          if (do.group) {
            object@net$prob[object@net$pval > thresh] <- 0
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob)[[3]]]
            prob <- object@net$prob[, , pairLR.use.name,
                                    drop = FALSE]
            prob.sum <- apply(prob > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          else {
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob.cell)[[3]]]
            prob.cell <- object@net$prob.cell[, , pairLR.use.name,
                                              drop = FALSE]
            prob.sum <- apply(prob.cell > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          if (length(pairLR.use$interaction_name) == 0) {
            stop(
              paste0(
                "There is no significant communication related with the input `pairLR.use`. Set `enriched.only = FALSE` to show non-significant signaling."
              )
            )
          }
        }
        LR.pair <- object@LR$LRsig[pairLR.use$interaction_name,
                                   c("ligand", "receptor")]
        geneL <- unique(LR.pair$ligand)
        geneR <- unique(LR.pair$receptor)
        geneL <- extractGeneSubset(geneL, object@DB$complex,
                                   object@DB$geneInfo)
        geneR <- extractGeneSubset(geneR, object@DB$complex,
                                   object@DB$geneInfo)
        feature.use <- c(geneL, geneR)
      }
      else {
        feature.use <- features
      }
      if (length(intersect(feature.use, rownames(data))) >
          0) {
        feature.use <- feature.use[feature.use %in% rownames(data)]
        data.use <- data[feature.use, , drop = FALSE]
      }
      else if (length(intersect(feature.use, colnames(meta))) >
               0) {
        feature.use <- feature.use[feature.use %in% colnames(meta)]
        data.use <- t(meta[, feature.use, drop = FALSE])
      }
      else {
        stop("Please check your input! ")
      }
      if (!is.null(cutoff)) {
        cat("Applying a cutoff of ", cutoff, "to the values...",
            "\n")
        data.use[data.use <= cutoff] <- 0
      }


      numFeature = length(feature.use)
      gg <- vector("list", numFeature)

      plot_ncols <- ceiling(numFeature / plot_nrows)
      plot_width <- 1 / plot_ncols
      plot_height <- 1 / plot_nrows
      annotations_pos <- vector("list", 0)
      for (i in 1:plot_nrows) {
        for (j in 1:plot_ncols) {
          x <- (j - 0.5) * plot_width
          y <- 1-(i - 0.95) * plot_height
          annotations_pos[[length(annotations_pos)+1]] <- list(
            x=x,
            y=y
          )
        }

      }
      annotations <- vector("list", numFeature)

      for (i in seq_len(numFeature)) {
        feature.name <- feature.use[i]
        df$feature.data <- data.use[i,]
        g <-
          ggplot(data = df, aes(x, y)) + geom_point(
            aes(colour = feature.data,cell_labels=cell_labels),
            alpha = alpha,
            size = point.size,
            shape = shape.by
          ) +
          scale_colour_gradientn(
            colours = colormap,
            guide = guide_colorbar(
              title = NULL,
              ticks = T,
              label = T,
              barwidth = 0.5
            ),
            na.value = "grey90"
          ) +
          theme(legend.position = "right") + theme(
            legend.title = element_blank(),
            legend.text = element_text(size = legend.text.size),
            legend.key.size = unit(0.15, "inches")
          ) + ggtitle(feature.name) +
          theme(plot.title = element_text(
            hjust = 0.5,
            vjust = 0,
            size = 10
          )) + theme(
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
          ) +
          xlab(NULL) + ylab(NULL) + theme(legend.key = element_blank())
        g <- g + coord_fixed() +
          scale_y_reverse()

        gg[[i]] <- g %>% plotly::ggplotly(height = 400)

        annotations[[i]] <- list(
          x=annotations_pos[[i]]$x,
          y=annotations_pos[[i]]$y,
          text = feature.name,
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "top",
          showarrow = FALSE,
          font = list(size = 16))
      }
      # print(annotations_pos)
      if (numFeature > 1) {
        gg <- plotly::subplot(gg, nrows = plot_nrows,margin = 0.02,shareX = FALSE,shareY = FALSE) %>%
          plotly::layout(title = '',annotations = annotations)
      }
      else {
        gg <- plotly::ggplotly(gg[[1]])
      }
    }
    else {
      # do.binary
      set_individual_legend <- function(plt) {
        # plt is a plotly plot obj
        plt_build <- plotly::plotly_build(plt)

        # get the num of traces
        len_legend <- length(plt_build$x$data)

        for (i in 1:len_legend) {
          # set legendgroup
          plt_build$x$data[[i]]$legendgroup <- feature.name
          # set legendtitle
          plt_build$x$data[[i]]$legendgrouptitle <- list(text=feature.name,font=list(size=12))
        }
        return(plt_build)
      }
      if (is.null(color.use)) {
        color.use <- ggPalette(4)
        color.use[4] <- "grey90"
      }
      color.use1 = color.use
      if (!is.null(signaling)) {
        res <- extractEnrichedLR(
          object,
          signaling = signaling,
          enriched.only = enriched.only,
          thresh = thresh
        )
        LR.pair <- object@LR$LRsig[res$interaction_name,
                                   c("ligand", "receptor")]
      }
      else if (!is.null(pairLR.use)) {
        if (is.character(pairLR.use)) {
          pairLR.use <- data.frame(interaction_name = pairLR.use)
        }
        if (enriched.only) {
          if (do.group) {
            object@net$prob[object@net$pval > thresh] <- 0
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob)[[3]]]
            prob <- object@net$prob[, , pairLR.use.name,
                                    drop = FALSE]
            prob.sum <- apply(prob > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          else {
            pairLR.use.name <-
              pairLR.use$interaction_name[pairLR.use$interaction_name %in%
                                            dimnames(object@net$prob.cell)[[3]]]
            prob.cell <- object@net$prob.cell[, , pairLR.use.name,
                                              drop = FALSE]
            prob.sum <- apply(prob.cell > 0, 3, sum)
            names(prob.sum) <- pairLR.use.name
            signaling.includes <- names(prob.sum)[prob.sum >
                                                    0]
            pairLR.use <- pairLR.use[pairLR.use$interaction_name %in%
                                       signaling.includes, , drop = FALSE]
          }
          if (length(pairLR.use$interaction_name) == 0) {
            stop(
              paste0(
                "There is no significant communication related with the input `pairLR.use`. Set `enriched.only = FALSE` to show non-significant signaling."
              )
            )
          }
        }
        LR.pair <- object@LR$LRsig[pairLR.use$interaction_name,
                                   c("ligand", "receptor")]
      }
      else {
        stop("Please input either `pairLR.use` or `signaling` for `binary` mode!")
      }
      geneL <- as.character(LR.pair$ligand)
      geneR <- as.character(LR.pair$receptor)
      complex_input <- object@DB$complex
      dataL <- computeExpr_LR(geneL, data, complex_input)
      dataR <- computeExpr_LR(geneR, data, complex_input)
      rownames(dataL) <- geneL
      rownames(dataR) <- geneR

      feature.use <- rownames(LR.pair)
      numFeature = nrow(LR.pair)

      if (is.null(cutoff)) {
        stop("A `cutoff` must be provided when plotting expression in binary mode! ")
      }
      gg <- vector("list", numFeature)

      # set subplot title pos
      plot_ncols <- ceiling(numFeature / plot_nrows)
      plot_width <- 1 / plot_ncols
      plot_height <- 1 / plot_nrows
      annotations_pos <- vector("list", 0)
      for (i in 1:plot_nrows) {
        for (j in 1:plot_ncols) {
          x <- (j - 0.5) * plot_width
          y <- 1-(i - 1) * plot_height
          annotations_pos[[length(annotations_pos)+1]] <- list(
            x=x,
            y=y
          )
        }

      }
      annotations <- vector("list", numFeature)

      for (i in seq_len(numFeature)) {
        feature.name <- feature.use[i]
        idx1 = dataL[i,] > cutoff
        idx2 = dataR[i,] > cutoff
        idx3 = idx1 & idx2
        group = rep("None", ncol(dataL))
        group[idx1] = geneL[i]
        group[idx2] = geneR[i]
        group[idx3] = "Both"
        group = factor(group, levels = c(geneL[i], geneR[i],
                                         "Both", "None"))
        color.use <- color.use1
        names(color.use) <- c(geneL[i], geneR[i], "Both",
                              "None")
        if (length(setdiff(levels(group), unique(group))) >
            0) {
          color.use <- color.use[names(color.use) %in% unique(group)]
          group = droplevels(group, exclude = setdiff(levels(group),
                                                      unique(group)))
        }
        df$feature.data <- group
        g <-
          ggplot(data = df, aes(x, y)) + geom_point(
            aes(colour = feature.data,cell_labels=cell_labels),
            alpha = alpha,
            size = point.size,
            shape = shape.by
          ) +
          scale_color_manual(values = color.use, na.value = "grey90") +
          theme(legend.position = "right") + theme(
            legend.title = element_blank(),
            legend.text = element_text(size = legend.text.size),
            legend.key.size = unit(0.15, "inches")
          ) + guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(feature.name) + theme(plot.title = element_text(
            hjust = 0.5,
            vjust = 0,
            size = 10
          )) + theme(
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
          ) +
          xlab(NULL) + ylab(NULL) + theme(legend.key = element_blank())
        g <- g + coord_fixed() +
          scale_y_reverse()

        # cat(feature.name)
        gg[[i]] <- g %>% plotly::ggplotly(
          type = 'scatter',
          mode='markers+text',
        ) %>% set_individual_legend()

        annotations[[i]] <- list(
          x=annotations_pos[[i]]$x,
          y=annotations_pos[[i]]$y,
          text = feature.name,
          xref = "paper",
          yref = "paper",
          xanchor = "center",
          yanchor = "top",
          showarrow = FALSE,
          font = list(size = 16))
      }

      if (numFeature > 1) {
        gg <- plotly::subplot(gg, nrows = plot_nrows,margin = 0.02,shareX = FALSE,shareY = FALSE) %>%
          plotly::layout(title = '',
                 annotations = annotations,
                 legend = list(tracegroupgap = 10,title=list(text=''))
          )
      }
      else {
        gg <- plotly::ggplotly(gg[[1]],
                               type = 'scatter',
                               mode = 'markers') %>%
          plotly::layout(legend = list(title = list(text = '')))
      }
    }
    return(gg)
  }



  # ##########################################################################
  # Shiny App's UI
  # ##########################################################################
  ui <- fluidPage(
    theme =  bslib::bs_theme(version = 5),
    # ##########################################################################
    # meta info of the HTML pages
    # ##########################################################################
    tags$head(
      # title
      tags$title("Interactive CellChat Explorer"),
      # icon
      tags$link(rel = "shortcut icon", type = "image/x-icon", href = "favicon.ico"),
      tags$link(rel="stylesheet",href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.5/font/bootstrap-icons.css"),
    ),
    tags$body(
      # ##########################################################################
      # logo and title of the website
      # ##########################################################################
      tags$nav(class="navbar navbar-light bg-light",
               div(class="container-fluid justify-content-center",
                   tags$a(
                     class="navbar-brand",href="http://www.cellchat.org/",
                     img(src="https://s2.loli.net/2023/08/08/2qjSoRACDtHByOY.png",class="d-inline",alt="",height="30"),
                     tags$p("Interactive CellChat Explorer",class="fs-1 d-inline")
                   )

               )),
      # ##########################################################################
      # 1.Basic exploration of spatial-resolved gene expression
      # ##########################################################################

      # Visualize cell groups and signaling expression
      h3(tags$i(class="bi bi-1-square-fill"),
         "Visualize cell groups and signaling expression",class="h3"),
      bslib::card(
        bslib::card_header(
          h6(tags$i(class="bi bi-bookmark"),
             "Dim Plot",class="h6")),
        layout_sidebar(
          sidebar = accordion(
            accordion_panel(
              "Numerical",
              icon = bsicons::bs_icon("sliders"),
              sliderInput(
                "dimplot_point_size",
                label = "Point size",
                min = 3,
                max = 8,
                step = 0.5,
                value = 3
              ),
              sliderInput(
                "dimplot_alpha",
                label = "Alpha",
                min = 0,
                max = 1,
                step = 0.2,
                value = 1
              ),
            )
          ),
          div(class="d-flex justify-content-center",
              plotly::plotlyOutput(outputId = "DimPlot",
                           width = 664,height = 498)
          )

        ),

      ),

      # gene expression distribution
      # https://shiny.posit.co/r/gallery/widgets/datatables-options/
      # https://shiny.posit.co/r/gallery/widgets/selectize-examples/
      navset_card_tab(
        title = h6(tags$i(class="bi bi-bookmark-dash"),
                   "Feature Plot",class="h6"),
        sidebar = NULL,
        # content
        nav_panel(
          title = "use gene names",
          layout_sidebar(
            sidebar = accordion(
              accordion_panel(
                "Select",
                icon = bsicons::bs_icon("menu-app"),
                selectizeInput(
                  inputId = 'selectize_gene_names',
                  label = 'Gene Names',
                  choices = NULL,
                  multiple = TRUE,
                  # options = list(maxItems = 4)
                ),
                numericInput(
                  "nrows_feature_plot1",
                  label = "nrows",
                  min = 1,
                  step = 1,
                  value = 1,
                ),
              ),

              accordion_panel(
                title = "Color",
                icon = tags$i(class="bi bi-palette-fill"),
                selectInput(
                  "direction_feature_plot1",
                  label = "direction",
                  choices = list("1"=1,"-1"=-1),
                  selected = 1,
                  multiple = F
                ),
                selectInput(
                  "palette_feature_plot1",
                  label = "palette",
                  choices = c(choices_palettes_diverging,choices_palettes_sequential),
                  selected = "Reds",
                  multiple = F
                ),
              ),

              accordion_panel(
                title = "Numerical",
                icon = bsicons::bs_icon("sliders"),
                numericInput(
                  "cut.off_feature_plot1",
                  label = "cut.off",
                  min = 0,
                  step = 0.1,
                  value = 0,
                ),
                sliderInput(
                  "point.size_feature_plot1",
                  label = "point.size",
                  min = 0,
                  max = 2,
                  step = 0.1,
                  value = 0.4
                )
              )

            ),
            # nav content
            div(class="d-flex justify-content-center",
                plotly::plotlyOutput(outputId = "gene_expression_distribution",width = 664,height = 498),
            ),
          )
        ),
        nav_panel(
          title = "use L-R pairs",
          layout_sidebar(
            sidebar = accordion(
              accordion_panel(
                "Select",
                icon = bsicons::bs_icon("menu-app"),
                selectizeInput(
                  inputId = 'selectize_pairLR_use',
                  label = 'pairLR_use',
                  choices = NULL,
                  multiple = T
                ),
                numericInput(
                  "nrows_feature_plot2",
                  label = "nrows",
                  min = 1,
                  step = 1,
                  value = 1,
                ),
                checkboxInput(
                  "do.binary_feature_plot",
                  label = "do.binary",
                  value = TRUE),
              ),
              accordion_panel(
                title = "Color",
                icon = tags$i(class="bi bi-palette-fill"),
                selectInput(
                  "direction_feature_plot2",
                  label = "direction",
                  choices = list("1"=1,"-1"=-1),
                  selected = 1,
                  multiple = F
                ),
                selectInput(
                  "palette_feature_plot2",
                  label = "palette",
                  choices = c(choices_palettes_diverging,choices_palettes_sequential),
                  selected = "Reds",
                  multiple = F
                ),
              ),

              accordion_panel(
                title = "Numerical",
                icon = bsicons::bs_icon("sliders"),
                numericInput(
                  "cut.off_feature_plot2",
                  label = "cut.off",
                  min = 0,
                  step = 0.1,
                  value = 0,
                ),
                sliderInput(
                  "point.size_feature_plot2",
                  label = "point.size",
                  min = 0,
                  max = 2,
                  step = 0.1,
                  value = 0.4
                )
              )
            ),
            # nav content
            div(class="d-flex justify-content-center",
                plotly::plotlyOutput(outputId = "gene_expression_distribution2",width = 664,height = 498)
            )
          ),
        )
      ),

      # ##########################################################################
      # 2.Examine signaling between cell groups
      # ##########################################################################
      h2(tags$i(class="bi bi-2-square-fill"),
         "Examine signaling between cell groups"),
      navset_card_tab(
        title = NULL,
        sidebar = NULL,
        nav_panel("Heatmap",
                  layout_sidebar(
                    sidebar = accordion(
                      h6("The number of interactions/interaction strength between any two cell groups",
                         class="h6"),
                      hr(),

                      accordion_panel(
                        "Select",
                        icon = bsicons::bs_icon("menu-app"),
                        selectInput(
                          "measure_heatmap",
                          label = "measurement",
                          choices = list("count" = "count", "weight" = "weight"),
                          selected = "count"
                        ),
                        selectInput(
                          "palette_heatmap",
                          label = "palette (sequential)",
                          choices = choices_palettes_sequential,
                          selected = "Blues"
                        ),
                        # Sets the order of colours in the scale. If 1, the default, colours are as output by RColorBrewer::brewer.pal(). If -1, the order of colours is reversed.
                        selectInput(
                          "direction_heatmap",
                          label = "direction",
                          choices = list(
                            "1"=1,
                            "-1"=-1
                          ),
                          selected = 1,
                        )

                        # refer to: https://ggplot2.tidyverse.org/reference/scale_brewer.html
                      ),
                    ),

                    # content
                    div(class="d-flex justify-content-center",
                        plotly::plotlyOutput(outputId = "netVisual_heatmap",width = 664,height = 498)
                    )
                  )
        ),

        # the enriched signaling among one selected pair of cell groups
        nav_panel("rankNet",
                  layout_sidebar(
                    sidebar = accordion(
                      h6("The enriched signaling",
                         class="h6"),
                      hr(),
                      accordion_panel(
                        "Select",
                        icon = bsicons::bs_icon("menu-app"),
                        selectInput(
                          "select1_cell_group",
                          label = "cell groups for sources.use",
                          choices = choices_cell_groups,
                          selected = choices_cell_groups[1],
                          multiple = TRUE
                        ),
                        selectInput(
                          "select2_cell_group",
                          label = "cell groups for targets.use",
                          choices = choices_cell_groups,
                          selected = choices_cell_groups[2],
                          multiple = TRUE
                        ),

                        # selectInput(
                        #   "measure_ranknet",
                        #   label = "measurement",
                        #   choices = list("count" = "count", "weight" = "weight"),
                        #   selected = "count"
                        # ),
                        selectInput(
                          "slot.name_ranknet",
                          label = "slot.name",
                          choices = list("net" = "net", "netP" = "netP"),
                          selected = "netP"
                        ),
                        # selectInput(
                        #   "palette_ranknet",
                        #   label = "palette (sequential)",
                        #   choices = choices_palettes_sequential,
                        #   selected = "Blues"
                        # ),
                      ),
                    ),

                    # content
                    plotly::plotlyOutput(outputId = "rankNet")
                  )
        ),
        nav_panel("Contribution Plot",
                  layout_sidebar(
                    sidebar = accordion(
                      h6("Contribution of each L-R pair to overall signaling",
                         class = "h6"),
                      hr(),
                      accordion_panel(
                        "Select",
                        icon = bsicons::bs_icon("menu-app"),
                        selectInput(
                          "pathway_contribution_plot",
                          label = "a pathway to show",
                          choices = choices_pathways,
                          selected = choices_pathways[1]
                        ),
                        selectInput(
                          "select3_cell_group",
                          label = "a cell group for sources.use",
                          choices = choices_cell_groups,
                          selected = choices_cell_groups[1]
                        ),
                        selectInput(
                          "select4_cell_group",
                          label = "a cell group for targets.use",
                          choices = choices_cell_groups,
                          selected = choices_cell_groups[2]
                        ),
                      ),

                      accordion_panel(
                        title = "Numerical",
                        icon = bsicons::bs_icon("sliders"),
                        sliderInput(
                          "font.size_contribution_plot",
                          label = "font.size",
                          min = 10,
                          max=30,
                          step = 5,
                          value = 20
                        )
                      )
                    ),

                    # content
                    plotOutput(outputId = "netAnalysis_contribution"),
                  )
        ),
      ),

      # Contribution of each L-R pair to overall signaling
      # width = 3.2inch, height = 1.5inch,
      # height might change dependent on dataset

      ## Examine individual signaling pathway
      ## (the following four plots will be appeared based on user's input)
      h2(tags$i(class="bi bi-3-square-fill"),
         "Examine individual signaling pathway"),
      navset_card_tab(
        title = h6("Plots",class="h6"),
        sidebar = accordion(
          selectizeInput(
            inputId = 'selectize_pathway',
            label = 'Select a pathway to show',
            choices = NULL,
            multiple = FALSE
          ),
          hr(),
          accordion_panel(
            title = "Circle plot",
            icon = tags$i(class="bi bi-circle-fill"),
            # edge.width.max = 5, vertex.size.max = 12, vertex.label.cex = 0.8
            sliderInput(
              "slider_Circle_plot_edge.width.max",
              label = "edge.width.max",
              min = 5,
              max = 15,
              value = 8,
              step = 1
            ),
            sliderInput(
              "slider_Circle_plot_vertex.size.max",
              label = "vertex.size.max",
              min = 8,
              max = 16,
              value = 12,
              step = 2),
            sliderInput(
              "slider_Circle_plot_vertex.label.cex",
              label = "vertex.label.cex",
              min = 1,
              max = 2,
              value = 1,
              step = 0.2
            ),
          ),
          accordion_panel(
            title = "Spatial plot",
            icon = tags$i(class="bi bi-layers-half"),
            # edge.width.max = 5, vertex.size.max = 1,
            # point.size = 2.5,
            # alpha.image = 0.2, vertex.label.cex = 5
            sliderInput(
              "slider_Spatial_plot_edge.width.max",
              label = "edge.width.max",
              min = 2,
              max = 8,
              value = 5,
              step = 1
            ),
            sliderInput(
              "slider_Spatial_plot_vertex.size.max",
              label = "vertex.size.max",
              min = 2,
              max = 8,
              value = 5,
              step = 1),
            sliderInput(
              "slider_Spatial_plot_vertex.label.cex",
              label = "vertex.label.cex",
              min = 5,
              max = 10,
              value = 8,
              step = 1
            ),

            sliderInput(
              "slider_Spatial_plot_point.size",
              label = "point.size",
              min = 1,
              max = 3,
              value = 2.4,
              step = 0.2
            ),
            sliderInput(
              "slider_Spatial_plot_alpha.image",
              label = "alpha.image",
              min = 0,
              max = 1,
              value = 0.2,
              step = 0.05
            ),
          ),
          accordion_panel(
            title = "Contribution of each L-R pair",
            icon = tags$i(class="bi bi-bar-chart-fill"),
          ),
        ),

        # nav tab
        nav_panel(
          title = "Circle plot",
          div(class="d-flex justify-content-center",
              plotOutput(outputId = "Circle_plot",
                         height = "780px",width = "580px")
          )
        ),
        nav_panel(
          title = "Spatial plot",
          div(class="d-flex justify-content-center",
              plotOutput(outputId = "Spatial_plot",
                         height = "780px",width = "580px")
          )
        ),
        nav_panel(
          title = "Contribution of each L-R pair",
          plotly::plotlyOutput(outputId = "LR_pair_contribution",
                       height = "900px")
        ),

      ),
      # body

    ),
    # page
  )
  # ##########################################################################
  # Shiny App's Server
  # ##########################################################################
  server <- function(input, output, session) {
    ############################################################################
    if (object@options$datatype == "RNA") {
      output$DimPlot <- plotly::renderPlotly({
        plotly_DimPlot(
          object,
          point.size = input$dimplot_point_size,
          alpha = input$dimplot_alpha,
        ) %>%
          plotly::config(toImageButtonOptions = list(
            format = "svg",
            filename = "DimPlot",
            width = 800,
            height = 600
          ))
      })
    } else {
      output$spatialDimPlot <- plotly::renderPlotly({
        plotly_spatialDimPlot(
          object,
          point.size = input$dimplot_point_size,
          alpha = input$dimplot_alpha,
        ) %>%
          plotly::config(toImageButtonOptions = list(
            format = "svg",
            filename = "spatialDimPlot",
            width = 800,
            height = 600
          ))
      })
    }


    observe({
      updateSelectizeInput(
        session,
        "selectize_gene_names",
        # selected = c("Wnt10a", "Fzd1", "Lrp6"),
        selected = choices_gene_names[1:2],
        # selected = c("Wnt10a", "Fzd1", "Lrp6","Ror2",
        #              "Nrp1","Nrp2","Bmpr2","Ret"),
        choices = choices_gene_names,
        server = TRUE
      )
    })
    # output$out6 <- renderPrint(input$selectize_gene_names)

    if (object@options$datatype == "RNA") {
      output$gene_expression_distribution <- plotly::renderPlotly(plotly_FeaturePlot(
        object,
        features = input$selectize_gene_names,
        plot_nrows = input$nrows_feature_plot1,
        point.size = input$point.size_feature_plot1,
        cutoff = input$cut.off_feature_plot1,
        color.heatmap = input$palette_feature_plot1,
        direction = input$direction_feature_plot1,
      ) %>%
        plotly::config(toImageButtonOptions = list(
          format = "svg",
          filename = "FeaturePlot (use gene names)",
          width = 600,
          height = 600
        ))
      )
    } else {
      output$gene_expression_distribution <- plotly::renderPlotly(plotly_spatialFeaturePlot(
        object,
        features = input$selectize_gene_names,
        plot_nrows = input$nrows_feature_plot1,
        point.size = input$point.size_feature_plot1,
        cutoff = input$cut.off_feature_plot1,
        color.heatmap = input$palette_feature_plot1,
        direction = input$direction_feature_plot1,
      ) %>%
        plotly::config(toImageButtonOptions = list(
          format = "svg",
          filename = "spatialFeaturePlot (use gene names)",
          width = 600,
          height = 600
        ))
      )
    }


    observe({
      updateSelectizeInput(
        session,
        "selectize_pairLR_use",
        selected = choices_pairLR_use[1],
        # selected = c("WNT10A_FZD1_LRP6","WNT10A_FZD10_LRP6","BMP2_BMPR1A_ACVR2A"),
        choices = choices_pairLR_use,
        server = TRUE
      )
    })
    # output$out7 <- renderPrint(input$selectize_pairLR_use)
    if (object@options$datatype == "RNA") {
      output$gene_expression_distribution2 <- plotly::renderPlotly({
        plotly_FeaturePlot(
          object,
          pairLR.use = input$selectize_pairLR_use,
          point.size = input$point.size_feature_plot2,
          do.binary = input$do.binary_feature_plot,
          cutoff = input$cut.off_feature_plot2,
          enriched.only = F,
          color.heatmap = input$palette_feature_plot2,
          direction = input$direction_feature_plot2,
          plot_nrows = as.numeric(input$nrows_feature_plot2)
        ) %>%
          plotly::config(toImageButtonOptions = list(
            format = "svg",
            filename = "FeaturePlot(use pairLRs)",
            width = 600,
            height = 600
          ))
      })
    } else {
      output$gene_expression_distribution2 <- plotly::renderPlotly({
        plotly_spatialFeaturePlot(
          object,
          pairLR.use = input$selectize_pairLR_use,
          point.size = input$point.size_feature_plot2,
          do.binary = input$do.binary_feature_plot,
          cutoff = input$cut.off_feature_plot2,
          enriched.only = F,
          color.heatmap = input$palette_feature_plot2,
          direction = input$direction_feature_plot2,
          plot_nrows = as.numeric(input$nrows_feature_plot2)
        ) %>%
          plotly::config(toImageButtonOptions = list(
            format = "svg",
            filename = "spatialFeaturePlot(use pairLRs)",
            width = 600,
            height = 600
          ))
      })
    }


    ############################################################################
    output$netVisual_heatmap <- plotly::renderPlotly({
      suppressWarnings({
        netVisual_heatmap(object,
                          measure = input$measure_heatmap,
        ) %>%
          plotly_netVisual_heatmap(
            palette.heatmap = input$palette_heatmap,
            direction.heatmap = input$direction_heatmap)
      })
    })

    output$rankNet <- plotly::renderPlotly({
      rankNet(
        object,
        mode = "single",
        measure = "weight",
        sources.use = input$select1_cell_group,
        targets.use = input$select2_cell_group,
        slot.name = input$slot.name_ranknet
      ) %>%
        plotly::ggplotly()
    })

    output$netAnalysis_contribution <- renderPlot({
      netAnalysis_contribution(
        object,
        signaling = input$pathway_contribution_plot,
        sources.use = input$select3_cell_group,
        targets.use = input$select4_cell_group,
        font.size = input$font.size_contribution_plot,
        font.size.title = input$font.size_contribution_plot,
      )
    },res = 96)
    ############################################################################
    observe({
      updateSelectizeInput(
        session,
        "selectize_pathway",
        selected = choices_pathways[1],
        choices = choices_pathways,
        server = TRUE
      )
    })
    output$Circle_plot <- renderPlot({
      netVisual_aggregate(
        object,
        signaling = input$selectize_pathway,
        layout = "circle",
        edge.width.max = input$slider_Circle_plot_edge.width.max,
        vertex.size.max = input$slider_Circle_plot_vertex.size.max,
        vertex.label.cex = input$slider_Circle_plot_vertex.label.cex
      )
    },res = 96)
    output$Spatial_plot <- renderPlot({
      netVisual_aggregate(
        object,
        signaling = input$selectize_pathway,
        layout = "spatial",
        edge.width.max = input$slider_Spatial_plot_edge.width.max,
        vertex.size.max = input$slider_Spatial_plot_vertex.size.max,
        vertex.label.cex = input$slider_Spatial_plot_vertex.label.cex,
        alpha.image = input$slider_Spatial_plot_alpha.image,
        point.size = input$slider_Spatial_plot_point.size,
      )
    })
    output$LR_pair_contribution <- plotly::renderPlotly({
      netAnalysis_contribution(
        object,
        signaling = input$selectize_pathway,
        font.size = 12,
        font.size.title = 14
      )
    })
    ############################################################################
  }


  # Running a Shiny app
  shinyApp(ui = ui, server = server,...)
}
