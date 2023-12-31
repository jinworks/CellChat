% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{barplot_internal}
\alias{barplot_internal}
\title{Bar plot for dataframe}
\usage{
barplot_internal(
  df,
  x = "cellType",
  y = "value",
  fill = "condition",
  legend.title = NULL,
  width = 0.6,
  title.name = NULL,
  xlabel = NULL,
  ylabel = NULL,
  color.use = NULL,
  remove.xtick = FALSE,
  stat.add = FALSE,
  stat.method = "wilcox.test",
  percent.y = FALSE,
  label.x = 1.5,
  show.legend = TRUE,
  x.lab.rot = FALSE,
  size.text = 10
)
}
\arguments{
\item{df}{a dataframe}

\item{x}{Name of one column to show on the x-axis}

\item{y}{Name of one column to show on the y-axis}

\item{fill}{Name of one column to compare the values}

\item{legend.title}{Name of legend}

\item{width}{bar width}

\item{title.name}{Name of the main title}

\item{xlabel}{Name of x label}

\item{ylabel}{Name of y label}

\item{color.use}{defining the color of bar plot;}

\item{remove.xtick}{whether remove x tick}

\item{stat.add}{whether adding statistical test}

\item{stat.method, label.x}{parameters for ggpubr::stat_compare_means}

\item{percent.y}{whether showing y-values as percentage}

\item{show.legend}{Whether show the legend}

\item{x.lab.rot}{Whether rorate the xtick labels}

\item{size.text}{font size}
}
\value{
ggplot2 object
}
\description{
Bar plot for dataframe
}
