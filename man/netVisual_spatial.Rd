% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{netVisual_spatial}
\alias{netVisual_spatial}
\title{Spatial plot of cell-cell communication network}
\usage{
netVisual_spatial(
  net,
  coordinates,
  meta,
  sample.use = NULL,
  color.use = NULL,
  title.name = NULL,
  sources.use = NULL,
  targets.use = NULL,
  idents.use = NULL,
  remove.isolate = FALSE,
  remove.loop = TRUE,
  top = 1,
  weight.scale = FALSE,
  vertex.weight = 20,
  vertex.weight.max = NULL,
  vertex.size.max = NULL,
  vertex.label.cex = 5,
  vertex.label.color = "black",
  edge.weight.max = NULL,
  edge.width.max = 8,
  edge.curved = 0.2,
  alpha.edge = 0.6,
  arrow.angle = 5,
  arrow.size = 0.2,
  alpha.image = 0.15,
  point.size = 1.5,
  legend.size = 5
)
}
\arguments{
\item{net}{A weighted matrix representing the connections}

\item{coordinates}{a data matrix in which each row gives the spatial locations/coordinates of each cell/spot}

\item{meta}{a data frame with at least two columns named `labels` and `samples`.
`meta$labels` is a vector giving the group label of each cell/spot. `meta$samples` is a factor vector defining the sample labels of each dataset. The length should be the same as the number of rows in `coordinates`.}

\item{sample.use}{the sample used for visualization, which should be the element in `meta$samples`.}

\item{color.use}{Colors represent different cell groups}

\item{title.name}{the name of the title}

\item{sources.use}{a vector giving the index or the name of source cell groups}

\item{targets.use}{a vector giving the index or the name of target cell groups.}

\item{idents.use}{a vector giving the index or the name of cell groups of interest.}

\item{remove.isolate}{whether remove the isolate nodes in the communication network}

\item{remove.loop}{whether remove the self-loop in the communication network. Default: TRUE}

\item{top}{the fraction of interactions to show}

\item{weight.scale}{whether scale the weight}

\item{vertex.weight}{The weight of vertex: either a scale value or a vector}

\item{vertex.weight.max}{the maximum weight of vertex; defualt = max(vertex.weight)}

\item{vertex.size.max}{the maximum vertex size for visualization}

\item{vertex.label.cex}{The label size of vertex}

\item{vertex.label.color}{The color of label for vertex}

\item{edge.weight.max}{the maximum weight of edge; defualt = max(net)}

\item{edge.width.max}{The maximum edge width for visualization}

\item{edge.curved}{Specifies whether to draw curved edges, or not.
This can be a logical or a numeric vector or scalar.
First the vector is replicated to have the same length as the number of
edges in the graph. Then it is interpreted for each edge separately.
A numeric value specifies the curvature of the edge; zero curvature means
straight edges, negative values means the edge bends clockwise, positive
values the opposite. TRUE means curvature 0.5, FALSE means curvature zero}

\item{alpha.edge}{the transprency of edge}

\item{arrow.angle}{The width of arrows}

\item{arrow.size}{the size of arrow}

\item{alpha.image}{the transparency of individual spots}

\item{point.size}{the size of spots}

\item{legend.size}{the size of legend}
}
\value{
an object of ggplot
}
\description{
Autocrine interactions are omitted on this plot. Group centroids may be not accurate for some data due to complex geometry.
The width of edges represent the strength of the communication.
}
