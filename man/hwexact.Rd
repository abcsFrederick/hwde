\name{hwexact}

\alias{hwexact}
\alias{snphwe}

\title{Exact test of HWE}

\description{This function implements an exact SNP test of
  Hardy-Weinberg equilibrium. It is a wrapper for \code{snphwe}}

\usage{hwexact(obs.hom1, obs.hets, obs.hom2, data, loci, observed)}

\arguments{
  \item{obs.hom1}{Observed number of individuals homozygous for the common allele.}
  \item{obs.hets}{Observed number of heterozygous individuals.}
  \item{obs.hom2}{Observed number of individuals homozygous for the rare allele.}
 \item{data}{If not \code{NULL}, specifies a data frame holding the data.}
 \item{loci}{If \code{data} is not \code{NULL}, column that holds some
   combination of symbols of the form \code{AA}, \code{Aa} and \code{aa}.}
 \item{observed}{Column that holds the frequencies, if \code{data} is not
   \code{NULL}.}
}

\value{The exact p-value testing departure from Hardy-Weinberg
  Equilibrium, conditional on the observed relative numbers of alleles
  of the two types}

\references{Wigginton, J.E.; Cutler, D.J. & Abecasis, G.R.; A Note on
  exact tests of Hardy-Weinberg equilibrium; American Journal of Human
  Genetics, 2005, 76, 887-93}

\author{Randall Johnson -- adapted from code by Wigginton et al}

\examples{hwexact(68, 28, 4)}

\keyword{univar}
\keyword{htest}
