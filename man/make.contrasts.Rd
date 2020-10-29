\name{make.contrasts}
\alias{make.contrasts}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Calculate contrasts required by the function}
\description{
  Calculates the contrasts that \code{hwde} uses in fitting models
  for genotypic disequilibrium.  At present these are limited to
  first order interaction terms between the distinct loci,
  possibly with distinct values foe each different level of a
  grouping factor.
}
\usage{
make.contrasts(data = data[, loci], allele.chars = letters)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{Data frame, with one column per locus, and as
    many rows as there are observations. Genotypes are represented
    using a two-letter code, e.g., AA, Aa, aa}
  \item{allele.chars}{By default the letters "a", "b", ...,
    one for each locus, are used.}
}
\details{
  Any pair of letters can be used for the two alleles.  The setting
  of \code{allele.chars} determines the coding used in the model
  formulae.  This function is called by \code{hwde}}
\value{
  \item{contrasts.df}{Data frame whose columns hold the contrasts that
    will be explanatory variables for the Poisson regression model.}
  \item{oset}{Offset vector for the GLM model.}
\item{list.columns}{List in which each element is a text string
that holds the codes for a sequence of model terms.}
}
\author{J.H. Maindonald}

\seealso{\code{identify.genotypes}, called by \code{make.contrasts}}

\examples{
make.contrasts(data=c("AA","Aa","aa"))
}

\keyword{models}
