<p align="center">
  <img width="200"  src="https://github.com/jinworks/CellChat/blob/main/CellChat_Logo.png">
</p>

# CellChat v2
## Version 2.1.0
- CellChat v2 now enables the [inferrence of cell-cell communication from multiple spatially resolved transcriptomics datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.html). Users should update the previously calculated individual CellChat object for spatial transcriptomics data analysis via `updateCellChat` function.
- We add [Frequently Asked Questions (FAQ) when analyzing spatially resolved transcriptomics datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html), particularly on how to apply CellChat to different technologies of spatial transcriptomics data, including sequencing-based and in-situ imaging-based readouts. In addition, we redefine the `scale.factors` for easier interpretation when applying other spatial technologies. In version 2.1.1, we change `scale.factors` to `spatial.factors`, but users still can run the old CellChat object with `scale.factors`. Users can also update the old CellChat object.
- In version 2.1.2, we change `object@meta$slices` to `object@meta$samples` in order to allow the identification of consistent communication across samples using the updated function `filterCommunication`. Users need to update the old CellChat object using `updateCellChat`. The ECM-Receptor signaling is assumed as diffusible signaling by default when analyzing spatial transcriptomics. 

CellChat v2 is an updated version that includes
- inference of spatially proximal cell-cell communication between interacting cell groups from spatially resolved transcriptomics
- expanded database CellChatDB v2 by including more than 1000 protein and non-protein interactions (e.g. metabolic and synaptic signaling) with rich annotations. A function named `updateCellChatDB` is also provided for easily updating CellChatDB. 
- new functionalities enabling easily interface with other computational tools for single-cell data analysis and cell-cell communication analysis
- interactive web browser function to allow exploration of CellChat outputs of spatially proximal cell-cell communication

For the version history and detailed important changes, please see the [NEWS file](https://github.com/jinworks/CellChat/blob/master/NEWS.md).

A step-by-step protocol for cell-cell communication analysis using CellChat is now available at [Jin et al., Nature Protocols 2024](https://www.nature.com/articles/s41596-024-01045-4). Please kindly cite this paper when using CellChat version >= 1.5. We greatly appreciate the users' support and suggestions that make it possible for us to update CellChat since we published the first version in the year of 2021. 

## Capabilities
In addition to infer the intercellular communication from any given scRNA-seq data and spatially resolved transcriptomics data, CellChat provides functionality for further data exploration, analysis, and visualization. 

- It can quantitatively characterize and compare the inferred cell-cell communication networks using an integrated approach by combining social network analysis, pattern recognition, and manifold learning approaches.
- It provides an easy-to-use tool for extracting and visualizing high-order information of the inferred networks. For example, it allows ready prediction of major signaling inputs and outputs for all cell populations and how these populations and signals coordinate together for functions.
- It enables comparative analysis of cell-cell communication across different conditions and identification of altered signaling and cell populations.
- It provides several visualization outputs to facilitate intuitive user-guided data interpretation.

## Installation

CellChat R package can be easily installed from Github using devtools:  

```
devtools::install_github("jinworks/CellChat")
```
**Please make sure you have installed the correct version of `NMF` and `circlize` package**. See instruction below. 

### Installation of other dependencies
- Install [NMF (>= 0.23.0)](http://renozao.github.io/NMF/devel/PAGE-INSTALLATION.html) using `install.packages('NMF')`. Please check [here](https://github.com/sqjin/CellChat/issues/16) for other solutions if you encounter any issue. You might can set `Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)` if it throws R version error. 
- Install [circlize (>= 0.4.12)](https://github.com/jokergoo/circlize) using `devtools::install_github("jokergoo/circlize")` if you encounter any issue.
- Install [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")` if you encounter any issue.
- Install UMAP python pacakge for dimension reduction: ```pip install umap-learn```. Please check [here](https://github.com/lmcinnes/umap) if you encounter any issue. 

Some users might have issues when installing CellChat pacakge due to different operating systems and new R version. Please check the following solutions:

- **Installation on Mac OX with R > 3.6**: Please re-install [Xquartz](https://community.rstudio.com/t/imager-package-does-not-work-in-r-3-6-1/38119).
- **Installation on Windows, Linux and Centos**: Please check the solution for [Windows](https://github.com/jinworks/CellChat/issues/5) and [Linux](https://github.com/jinworks/CellChat/issues/131).  



## Tutorials
Please check the tutorial directory of the repo. Example datasets are publicly available at [figshare](https://figshare.com/projects/Example_data_for_cell-cell_communication_analysis_using_CellChat/157272). 
Please check the [Jin et al., Nature Protocols 2024](https://www.nature.com/articles/s41596-024-01045-4) for a comprehensive protocol of cell-cell communication analysis using CellChat. 


### Analysis of single-cell transcriptomics data
- [Full tutorial for CellChat analysis of a single dataset with detailed explanation of each function](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html)
- [Full tutorial for comparison analysis of multiple datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html)
- [Comparison analysis of multiple datasets with different cellular compositions](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions.html)

### Analysis of spatially resolved omics data
- [Brief tutorial for CellChat analysis of a single spatially resolved transcriptomic dataset](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.html)
- [Brief tutorial for CellChat analysis of multiple spatially resolved transcriptomic datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.html)
- [Brief tutorial for CellChat analysis of spatial multiomics data](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_spatial_multiomics_data.html)
- [Frequently Asked Questions when analyzing spatially resolved transcriptomics datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html)

### Additional utilities
- [Interface with other single-cell analysis toolkits (e.g., Seurat, SingleCellExperiment, Scanpy)](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Interface_with_other_single-cell_analysis_toolkits.html)
- [Tutorial for updating ligand-receptor database CellChatDB](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html)


## Web-based “CellChat Explorer” 

We build a user-friendly web-based “[CellChat Explorer](http://www.cellchat.org/)” that contains two major components:
- **Ligand-Receptor Interaction Explorer** that allows easy exploration of our novel ligand-receptor interaction database, a comprehensive recapitulation of known molecular compositions including multimeric complexes and co-factors. *Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both **human and mouse***. Of note, this Explorer currently only shows the original CellChatDB, but did not include the new interactions in CellChatDB v2. 
- **Cell-Cell Communication Atlas Explorer** that allows easy exploration of the cell-cell communication for any given scRNA-seq dataset that has been processed by our R toolkit CellChat.

We also developed an Interactive Web Browser that allows exploration of CellChat outputs of spatially proximal cell-cell communication using a built-in function `runCellChatApp`, and a [standalone CellChat Shiny App](https://github.com/sqjin/CellChatShiny) for the above Cell-Cell Communication Atlas Explorer. 


## Help, Suggestion and Contribution
If you have any question, comment or suggestion, please use github issue tracker to report coding related [issues](https://github.com/jinworks/CellChat/issues) of CellChat. 

### Before reporting an issue
- First **check the GitHub [issues](https://github.com/jinworks/CellChat/issues)** to see if the same or a similar issues has been reported and resolved. This relieves the developers from addressing the same issues and helps them focus on adding new features!
- The best way to figure out the issues is **running the sources codes** of the specific functions by yourself. This will also relieve the developers and helps them focus on the common issues! I am sorry, but I have to say I have no idea on many errors except that I can reproduce the issues. 
- Minimal and **reproducible example** are required when filing a GitHub issue. In certain cases, please share your CellChat object and related codes to reproduce the issues. 
- Users are encouraged to discuss issues and bugs using the github [issues](https://github.com/jinworks/CellChat/issues) instead of email exchanges.

### Contribution
CellChat is an open source software package and any contribution is highly appreciated! 

We use GitHub's [Pull Request](https://github.com/jinworks/CellChat/pulls) mechanism for reviewing and accepting submissions of any contribution. Issue a pull request on the GitHub website to request that we merge your branch's changes into CellChat's master branch. Be sure to include a description of your changes in the pull request, as well as any other information that will help the CellChat developers involved in reviewing your code. 

## System Requirements
- Hardware requirements: CellChat package requires only a standard computer with enough RAM to support the in-memory operations.

- Software requirements: This package is supported for macOS, Windows and Linux. The package has been tested on macOS: Ventura (version 13.5) and Windows 10. Dependencies of CellChat package are indicated in the Description file, and can be automatically installed when installing CellChat pacakge. CellChat can be installed on a normal computer within few mins.


# About CellChat and CellChatDB
CellChat is an R package designed for inference, analysis, and visualization of cell-cell communication from single-cell and spatially resolved transcriptomics. CellChat aims to enable users to identify and interpret cell-cell communication within an easily interpretable framework, with the emphasis of clear, attractive, and interpretable visualizations.  

CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in mutiple species, leading to a comprehensive recapitulation of known molecular interaction mechanisms including multi-subunit structure of ligand-receptor complexes and co-factors.

If you use CellChat or CellChatDB in your research, please considering citing our papers: 
- [Suoqin Jin et al., CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics, Nature Protocols 2024](https://www.nature.com/articles/s41596-024-01045-4) [CellChat v2] (Please kindly cite this paper when using CellChat version >= 1.5)
- [Suoqin Jin et al., Inference and analysis of cell-cell communication using CellChat, Nature Communications 2021](https://www.nature.com/articles/s41467-021-21246-9) [CellChat v1]


<p align="center">
  <a href="https://clustrmaps.com/site/1bpq2">
     <img width="200"  src="https://clustrmaps.com/map_v2.png?cl=ffffff&w=a&t=n&d=42WqeykSXznN_NSaBlpf6CtSXQxhqmIs6QusUsguFdY" />
   </a>
</p>
<p align="center">
  <a href="#">
     <img src="https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fjinworks%2FCellChat&labelColor=%233499cc&countColor=%2370c168" />
   </a>
</p>





