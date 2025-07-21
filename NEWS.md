# Version 2.1.2 (2024-2-6)
- Update `filterCommunication` to allow identifying consistent cell-cell communication across samples when mutiple samples/replicates/batches are merged as an input for CellChat analysis.
- Update `computeCommunProb` where the ECM-receptor signaling is now assumed as diffusible signaling when analyzing spatial transcrptomics.
- Update 'updateCellChat' to change `object@meta$slices` to `object@meta$samples` in order to identify consistent cell-cell communication across samples.
- Typos in CellChatDB were fixed.

## Updated functions with minor changes
Functions that have been updated for analyzing spatial transcriptomics and perform comparison analysis, including  `createCellChat`, 'updateCellChat', `computeCommunProb`,`filterCommunication`, `spatialFeaturePlot`, `spatialDimPlot`, `netVisual_spatial`, `netVisual_aggregate`, `netVisual`, `computeRegionDistance`, `StackedVlnPlot`.

# Version 2.1.1 (2023-12-12)
- Enable to flexibly infer contact-dependent or juxtacrine signaling from any type of spatial transcriptomics data by defining the `contact.range` in `computeCommunProb`. 
- Change `scale.factors` to `spatial factors`.

## Updated functions with minor changes
Functions that have been updated for analyzing spatial transcriptomics and perform comparison analysis, including  `createCellChat`, 'updateCellChat', `computeCommunProb`,`computeRegionDistance`, `netMappingDEG`, `computeEnrichmentScore`, `netVisual_bubble`, `identifyOverExpressedInteractions`, `identifyOverExpressedGenes`.

# Version 2.1.0 (2023-11-26)
- CellChat v2 now enables the [inferrence of cell-cell communication from multiple spatially resolved transcriptomics datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.html). Users should update the previously calculated individual CellChat object for spatial transcriptomics data analysis via `updateCellChat` function.
- We add [Frequently Asked Questions (FAQ) when analyzing spatially resolved transcriptomics datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html), particularly on how to apply CellChat to different technologies of spatial transcriptomics data, including sequencing-based and in-situ imaging-based readouts. In addition, we redefine the `scale.factors` for easier interpretation when applying other spatial technologies.

## Updated functions with changes
Functions that have been updated for analyzing spatial transcriptomics, including `CellChat-class`, `createCellChat`, 'updateCellChat', `computeCommunProb`,`computeRegionDistance`, `computeCellDistance`, `subsetDB`, `spatialFeaturePlot`, `spatialDimPlot`, `netVisual_spatial`, `netVisual_aggregate`, `netVisual`. 

# CellChat 2.0.0 (2023-11-01)

CellChat v2 is an updated version that includes
- inference of spatially proximal cell-cell communication between interacting cell groups from spatially resolved transcriptomics
- expanded database CellChatDB v2 by including more than 1000 protein and non-protein interactions such as metabolic and synaptic signaling, and by adding rich annotations to each ligand-receptor pair. A function named `updateCellChatDB` is also provided for easily updating CellChatDB. 
- new functionalities enabling easily interface with other computational tools for single-cell data analysis and cell-cell communication analysis
- interactive web browser function to allow exploration CellChat outputs of spatially proximal cell-cell communication

## Updated functions with minor changes
Functions that have been updated for analyzing spatial imaging data, such as `netVisual_circle`, `netVisual`, `netVisual_aggregate`, `netVisual_individual`, `pieChart`, `showDatabaseCategory`, `subsetDB`, `identifyOverExpressedGenes`, `computeRegionDistance`, `computeCommunProb`, ''

## New added functions
`runCellChatApp`, `updateCellChatDB`, `updateCCC_score`, `findEnrichedSignaling`,  `spatialFeaturePlot`, `spatialDimPlot`


# CellChat 1.6.0 (2022-11-12)

CellChat is now applicable to spatial imaging data. We showcase its application to 10X Visium data. When spatial locations of spots/cells are available, CellChat infers spatial-informed cell-cell communication between interacting cell groups. CellChat restricts cell-cell communication within the maximum interaction/diffusion length of molecules. 

A brief [tutorial](https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_imaging_data.html) for spatial imaging data analysis is available in the tutorial directory. CellChat's various functionality can be used for further data exploration, analysis, and visualization. 

We have redesigned the structure of CellChat object. When loading previous CellChat object (version < 1.6.0), please update the object via `updateCellChat`. 

## Updated functions with minor changes
Functions that have been updated for analyzing spatial imaging data, such as `CellChat-class`, `createCellChat`, `updateCellChat`, `computeCommunProb`, `netVisual`, `netVisual_aggregate`, `netVisual_individual`

## New added functions
`netVisual_spatial`, `computeRegionDistance`

# CellChat 1.5.0 (2022-08-06)
We have now presented our comparison framework for systematically detecting dysregulated cell-cell communication across biological conditions, and
then utilized it to study the aging-induced signaling changes during skin wound healing. Our results not only present general communication rules and signaling mechanisms in wound healing associated with aging, but also provide a paradigm for other researchers to study cell-cell communication in other contexts. Please check out [our paper (Vu#, Jin#, Sun# et al., Cell Reports, 2022)](https://doi.org/10.1016/j.celrep.2022.111155) for the detailed methods and applications.

# CellChat 1.4.0 (2022-05-07)
## Updated functions with minor changes
A number of functions have been updated with enhanced functionalities, such as `rankNet`, `netVisual_aggregate`, `netVisual_circle`, 'plotGeneExpression', `netVisual_embedding`.  

## New added functions
`computeEnrichmentScore`, 'netVisual_barplot', 'barPlot'

# CellChat 1.1.3 (2021-08-14)
## Updated functions with minor changes
`netAnalysis_computeCentrality`, `netVisual_embeddingPairwise`, `netAnalysis_signalingChanges_scatter`
NB: The method for computing the 'influencer' metric in the function `netAnalysis_computeCentrality` has been changed due to the issue that "lower" mode gives different results when re-ordering the matrix. The results changed, but it looks like the dominant patterns do not change, i.e., the top cell groups ranked based on this metrix retain the same.  

# CellChat 1.1.2 (2021-07-10)
## New added functions
* Add `netAnalysis_diff_signalingRole_scatter` for 2D visualization of differential signaling roles of each cell group when comparing mutiple datasets.

## Updated functions with minor changes
`netVisual`, `netVisual_aggregate`, `netVisual_individual`, `netVisual_hierarchy1`, `netVisual_hierarchy2`,`netVisual_circle`,`createCellChat`,`netAnalysis_computeCentrality` ,`netAnalysis_signalingRole_scatter`, `netAnalysis_signalingChanges_scatter`

## Important changes
* Add `thresh = 0.05` in `netAnalysis_computeCentrality` to only consider the significant interactions. This will slightly change (very likely quantitative instead of qualitative change) the results computed by previous version of CellChat. 

* In the updated `netAnalysis_computeCentrality`, we now also compute unweighted outdegree (i.e., the total number of outgoing links) and indegree (i.e., the total number of incoming links). 

* `netAnalysis_signalingRole_scatter` and `netAnalysis_signalingChanges_scatter` now support the comparison of the total number of outgoing and incoming links. 

* Change the default setting for visualizing cell-cell communication network: 1) using `circle` plot instead of `hierarchy`; 2) using the same node size instead of different size (setting `vertex.weight = NULL` will give different size as in previous version of CellChat).  


# CellChat 1.1.1 (2021-06-19)
## New added functions
* Add `updateClusterLabels` to update cell cluster labels without re-running the time-consuming function `computeCommunProb`.
## Updated functions with minor changes
`netAnalysis_signalingRole_scatter`, `setIdent`, `netVisual_diffInteraction`, `sketchData`, `identifyOverExpressedGenes`,`netEmbedding`,`runUMAP`

# CellChat 1.1.0 (2021-04-26)
## New added functions
* Add `subsetCellChat` to create a object using a portion of cells
* Add `computeAveExpr` to compute the average expression per cell group
* Add `sketchData` to downsample the single cell data for faster calculation
* Add `netAnalysis_signalingChanges_scatter` to identify the signaling changes associated with one cell group

## Important changes
* `computeCommunProb` now supports a faster calculation and displays a ProgressBar
* `computeCommunProbPathway` now returns the significant pathways that are ordered based on the total communication probabilities

## Publish a new release
* The first release was published as version 1.0.0 before updating to version 1.1.0 

# CellChat 1.0.0 (2021-02-17)

* CellChat paper is now officially published [(Jin et al., Nature Communications, 2021)](https://www.nature.com/articles/s41467-021-21246-9). Compared to the preprint, we have now experimentally validated CellChat's predictions on embryonic skin using RNAscope technique, applied CellChat to a human diseased skin dataset and updated many others. 
* We have now developed a [standalone CellChat Shiny App](https://github.com/sqjin/CellChatShiny) for interactive exploration of the cell-cell communication analyzed by CellChat. Want to share your results with your collaborators like biologists for further exploration? Try it out! 


# CellChat 0.5.0 (2021-01-05)

* Slight changes of CellChat object (Please update your previously calculated CellChat object via `updateCellChat()`)
* Enhanced documentation of functions and tutorials (use `help()` to check the documentation, e.g., `help(CellChat)`)
* New features for comparison analysis of multiple datasets
* Support for creating a new CellChat object from Seurat V3 or SingleCellExperiment object


