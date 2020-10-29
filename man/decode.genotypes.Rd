\name{decode.genotypes}
\alias{decode.genotypes}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Creates contrasts that relate to a single locus}
\description{
Vectors ma, maa and oset (offset) are returned.
}
\usage{
decode.genotypes(genotype)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{genotype}{The vector \code{genotype} holds two-letter codes
  for the three genotypes.  For example, the values may be AA, Aa and aa.}
}
\value{
  \item{oset}{The offset values that would be appropriate, in the
    multiplicative version of the model, if there was just this one
  locus.}
  \item{ma}{A contrast for the Hardy-Weinberg model, at this locus.}
  \item{maa}{A contrast that measures departure from the Hardy-Weinberg
    model, at this locus.}
  \item{types}{A vector of length three whose elements are the
    two-letter codes used for the three genotypes.}
}
\author{J.H. Maindonald}
\note{Called by \code{make.contrasts}}
\seealso{\code{hwde}}

\examples{
decode.genotypes(rep(c("AA","Aa","aa"),2))
}

\keyword{models}
