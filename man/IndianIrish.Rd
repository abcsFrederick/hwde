\name{IndianIrish}
\alias{IndianIrish}
\non_function{}
\title{Observed genotype frequencies at MN and S loci, for 2 populations}
\usage{data(IndianIrish)}
\description{
The \code{IndianIrish} data frame has 18 rows and 4 columns.
The data are genotype frequencies for two locations, for
Xavante Indian and Irish populations respectively
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Population}{Factor with levels:
    \code{Indian} and \code{Irish} }
    \item{locus1}{Factor with levels:
    \code{MM}, \code{MN} and \code{NN}}
    \item{locus2}{Factor with levels:
    \code{SS}, \code{Ss} and \code{ss}}
    \item{Observed}{a numeric vector giving the frequency for each
    category of the tale}
  }
}
\source{
  Mourant et al (1977) and Huttley and Wilson (2000).
}
\references{
  1. Huttley, G.A. and Wilson, S.R. 2000. Testing for concordant
  equilibria between population samples. \emph{Genetics} \bold{156}, 2127-2135. \cr
  2. Mourant, A.E., Kopec, A.C. and Domaniewska-Sobczak, K. 1976.
\emph{The Distribution of the Human Blood Groups and Other Polymorphisms.}
Oxford University Press.\cr
  3. Weir, B.S. 1996.  \emph{Genetic Data Analysis II.} Sinauer.
}
\seealso{\code{\link{hwde}}}

\examples{
data(IndianIrish)
hwde(data=IndianIrish)
}
\keyword{datasets}
