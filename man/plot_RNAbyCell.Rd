% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_RNAbyCell}
\alias{plot_RNAbyCell}
\title{Plot RNA distribution by cell type}
\usage{
plot_RNAbyCell(data, time, cols_plot, cell = NULL, facet = TRUE)
}
\arguments{
\item{data}{table with column data containing percent reads by RNA type}

\item{time}{time hpi to filter on}

\item{cols_plot}{column numbers identifying percent RNA captured for RNAs of interest}

\item{cell}{set to null, can provide cell type of interest}

\item{facet}{set to FALSE if providing specific cell type}
}
\description{
Plot RNA distribution by cell type
}
\examples{

plot_RNAbyCell(final_table, "12", 13:20, "B6", facet = FALSE)

plot_RNAbyCell(final_table, "12", 13:20)

}
