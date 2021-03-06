% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage.R
\name{normalize_table}
\alias{normalize_table}
\title{Normalize coverage table using the sum of aligned reads per library}
\usage{
normalize_table(library_sizes, table)
}
\arguments{
\item{x}{path to depth files containing all aligned reads in each libary}

\item{y}{complete table with position and count data for specific RNA of interest}
}
\value{
table with position, count, normalized count, and dinucleotide data
}
\description{
Normalize coverage table using the sum of aligned reads per library
}
\examples{
complete_table(path, table)


}
