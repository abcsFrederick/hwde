\name{hwde}
\alias{hwde}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fit relevant models, and test for various types of departure from
  Hardy-Weinberg equilibrium.  Allows only 2 alleles per locus.  The
  number of loci is arbitrary.
}
\description{
Fits models for genotypic disequilibria, as described in Huttley and
Wilson (2000), Weir (1996) and Weir and Wilson (1986).
}
\usage{
hwde(data = hwde::IndianIrish, gp = "Population", termlist = NULL, refmodel = NULL,
loci = paste("locus", 1:(dim(data)[2] - 1), sep = ""), observed = "Observed",
keep.models = FALSE, aovtable.print = TRUE, group.terms = TRUE, allele.chars = letters)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{Must have a column of frequencies, by default called
    \code{Observed}, and one or more columns giving genotype
    information,
    in the form AA, Aa, aa. or e.g., MM, MN, NN.  (The choice of letters
    is arbitrary.)  Additionally, there may be a column that gives
    information on groupings, by default called \code{Population}.
  }
  \item{gp}{Gives the name of the column, if any, that has information
    on groups within the data.}
  \item{termlist}{Use to specify a user-defined sequence of models.
  See the vignette \bold{hwde.pdf} or \bold{hwde.html}}
  \item{refmodel}{For each model in \code{termlist}, specifies a
    reference model, which will be updated to include the additional terms.}
  \item{loci}{Gives name(s) of columns that hold information on
  genotypes. By default, these are taken to be \code{locus1},
  \code{locus2}, etc.}
  \item{observed}{Name (by default \code{Observed}) of the column that
holds the frequenceies.}
  \item{keep.models}{Should a list be returned that holds the full
    sequence of models that were fitted?}
  \item{aovtable.print}{Should the anova table be printed?}
  \item{group.terms}{Should model terms be grouped according to
    hierarchy,
  for the anova table?}
\item{allele.chars}{A sequence of letters used to code for the loci.
By default \emph{a}, \emph{b}, \emph{c}, ... are used}
}
\details{
  See the document \bold{hwde.pdf} or \bold{hwde.html} for details.
  See the references (below) for information on the interpretation of
  model parameters.
}
\value{
  \item{anovatab}{anova (analysis of deviance) table}
  \item{data.df}{Data, and contrasts used in fitting the various
    models.}
  \item{aovtab.terms}{This string holds, for each model that is
    fitted. the terms that have appeared in the model formula.
    The text strings for the distinct models are concatenated.}
  \item{models}{Optionally, this holds the complete sequence of
    qmodel objects that were fitted}
}
\references{
  1. Huttley, G.A. and Wilson, S.R. 2000. Testing for concordant
  equilibria between population samples. \emph{Genetics} \bold{156}:2127-2135. \cr
  2. Weir, B.S. 1996.  \emph{Genetic Data Analysis II.} Sinauer.\cr
  3. Weir, B.S. and Wilson, S.R. 1986.  Log-linear models for
linked loci. \emph{Biometrics} \bold{42}:665-670.
}
\author{J.H. Maindonald}

\seealso{\code{\link{make.contrasts}}, \code{\link{decode.genotypes}} }

\examples{
data(IndianIrish)
hwde(data=IndianIrish)
data(mendelABC)
hwde(data=mendelABC, loci=c("seedshape", "cotylcolor", "coatcolor"))
}
\keyword{models}
