#' Observed genotype frequencies at MN and S loci, for 2 populations
#' 
#' The \code{IndianIrish} data frame has 18 rows and 4 columns. The data are genotype frequencies for two locations, for Xavante Indian and Irish populations respectively.
#' 
#' @format This data frame contains the following columns:
#' \describe{
#'   \item{Population}{Factor with levels: \code{Indian} and \code{Irish} }
#'   \item{locus1}{Factor with levels: \code{MM}, \code{MN} and \code{NN}}
#'   \item{locus2}{Factor with levels: \code{SS}, \code{Ss} and \code{ss}}
#'   \item{Observed}{a numeric vector giving the frequency for each category of the tale}
#' }
#' 
#' @source Mourant et al (1977) and Huttley and Wilson (2000).
#' 
#' @references 
#' \itemize{
#'   \item Huttley, G.A. and Wilson, S.R. 2000. Testing for concordant equilibria between population samples. \emph{Genetics} \bold{156}, 2127-2135.
#'   \item Mourant, A.E., Kopec, A.C. and Domaniewska-Sobczak, K. 1976. \emph{The Distribution of the Human Blood Groups and Other Polymorphisms.} Oxford University Press.
#'   \item Weir, B.S. 1996.  \emph{Genetic Data Analysis II.} Sinauer.
#' }
#' 
#' @seealso \code{\link{hwde}}
#' 
#' @examples 
#' data(IndianIrish)
#' hwde(data=IndianIrish)
"IndianIrish"


#' Mendel's F2 trifactorial data for seed shape (A: round or wrinkled), cotyledon color (B: albumen yellow or green), and seed coat color (C: grey-brown or white)
#' 
#' The \code{mendel3} data frame has 27 rows and 4 columns. Data are from Mendel (1886), and are reproduced in Fisher (1936) and Weir (1996).
#' 
#' @format This data frame contains the following columns:
#' \describe{
#'   \item{seedshape}{Factor with levels: \code{AA}, \code{Aa} and \code{aa}}
#'   \item{cotylcolor}{Factor with levels: \code{BB}, \code{Bb} and \code{bb}}
#'   \item{coatcolor}{Factor with levels: \code{CC}, \code{Cc} and \code{cc}}
#'   \item{Observed}{a numeric vector that holds the frequencies.}
#'   }
#' 
#' @details 
#' The data are reviewed in detail in Fisher (1936). For a brief discussion, and references to work that revisits Fisher's conclusions, see Weir (1996).
#' 
#' @source  Data are from Mendel (1886), and are reproduced in Fisher (1936) and Weir (1996).
#' 
#' @references 
#' \itemize{
#'   \item Fisher, R.A. 1936. Has Mendel's work been rediscovered? \emph{Annals of Science} \bold{1}:115-137.
#'   \item Mendel, G. 1886. Versuche uber Pflanzen-Hybriden. Verhandlugen des Naturforschenden Vereines in Brunn \bold{4}:3-47. (An English translation. with annotations, is available from \url{http://www.esp.org/foundations/genetics/classical/gm-65.pdf} NB also the English translation by Royal Horticultural Society of London, reprinted in Peters, J.A. 1959. \emph{Classic Papers in Genetics.} Prentice-Hall.)
#'   \item Weir, B.S. 1996.  \emph{Genetic Data Analysis II.} Sinauer.
#' }
#' 
#' @examples 
#' ## Lay table out as a 3D array, as in Fisher (1936)
#' abc <- aperm(array(mendelABC[,4], dim=c(3,3,3)), c(1,3,2))
#' dimnames(abc) <- list(B=c('BB','Bb','bb'), 
#'                       A=c('AA','Aa','aa'),
#'                       C=c('CC','Cc','cc'))
#' abc
#' 
#' ## Fit Hardy-Weinberg disequilibium model
#' hwde(mendelABC, loci=c("seedshape","cotylcolor","coatcolor"))
"mendelABC"