% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_plot.R
\name{model_plot}
\alias{model_plot}
\title{Outlier visualization}
\usage{
model_plot(
  df,
  output.type = c("single", "combo"),
  labels = NULL,
  label.size = NULL,
  write.out = c("yes", "no"),
  num.col = NULL,
  h = NULL,
  w = NULL,
  output.dir = NULL,
  file.name = NULL
)
}
\arguments{
\item{df}{The dataframe generated with \code{random_frag_model()}}

\item{output.type}{If a grid of plots is desired define as \code{"combo"}, if individual plots are desired define as \code{"single"}.}

\item{labels}{The labels to be assigned to each panel. Use \code{"AUTO"} to automatically assign alphabetical labels starting with "A". For custom labels, enter them as a character vector.}

\item{label.size}{The size of the labels on the combined plot output.}

\item{write.out}{Determines whether or not to save the generated plot. To save the plots, define as \code{"yes"}.}

\item{num.col}{The desired number of columns in the cowplot output}

\item{h}{The desired height of the generated plot}

\item{w}{The desired width of the generated plot}

\item{output.dir}{The directory in which the files will be output. If the working directory is the desired directory, enter \code{NULL}.}

\item{file.name}{The name of the output file. Enter as a character string.}
}
\description{
Generates separate plots of genic integration frequency vs. gene length for each sample being analyzed.
Depicts model boundaries within gray shaded areas and colors upper and lower outliers blue and red, respectively.
}
\examples{
model_plot()

}
