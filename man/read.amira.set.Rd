% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.amira.set.R
\name{read.amira.set}
\alias{read.amira.set}
\title{read.amira.set}
\usage{
read.amira.set(name.file, nland)
}
\arguments{
\item{name.file}{character: path of a landmarkAscii file}

\item{nland}{numeric: number of landmark sampled in Amira, if is set on "auto" it will be automatically recognized}
}
\value{
array.set numeric: a kx3x1 array with landmark coordinates
}
\description{
This function converts a landmarkAscii file set in a kx3x1 array
}
\author{
Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
}
