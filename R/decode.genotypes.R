#' Creates contrasts that relate to a single locus
#' 
#' Vectors ma, maa and oset (offset) are returned.
#' 
#' @param genotype The vector \code{genotype} holds two-letter codes for the three genotypes. For example, the values may be AA, Aa and aa.
#' 
#' @return
#' \itemize{
#'     \item{oset}{The offset values that would be appropriate, in the multiplicative version of the model, if there was just this one locus.}
#'     \item{ma}{A contrast for the Hardy-Weinberg model, at this locus.}
#'     \item{maa}{A contrast that measures departure from the Hardy-Weinberg model, at this locus.}
#'     \item{types}{A vector of length three whose elements are the two-letter codes used for the three genotypes.}
#' }
#' 
#' @author J.H. Maindonald
#' 
#' @note Called by \code{make.contrasts}
#' 
#' @seealso \code{hwde}
#' 
#' @examples 
#' decode.genotypes(rep(c("AA","Aa","aa"),2))
#' 
#' @export
    
"decode.genotypes" <-
function (genotype) 
{
    gname <- deparse(substitute(genotype))
    genotype <- factor(genotype)
    gtypes <- levels(genotype)
    charswide <- nchar(gtypes)
    nabc <- length(gtypes)
    if ((nabc != 3) | any(charswide != 2)) 
        stop(paste("Illegal codes", paste(gtypes, collapse = ", "), 
            "at locus", gname))
    ch1 <- substring(gtypes, 1, 1)
    ch2 <- substring(gtypes, 2, 2)
    AA <- gtypes[ch2 == ch1][1]
    Aa <- gtypes[ch1 != ch2]
    aa <- gtypes[ch2 == ch1][2]
    ord <- match(c(AA, Aa, aa), gtypes)
    id <- c(2, 1, 0)
    names(id) <- c(AA, Aa, aa)
    n <- length(genotype)
    idloc <- c(1, 2, 1)
    names(idloc) <- c(AA, Aa, aa)
    oset <- idloc[as.character(genotype)]
    ma <- id[as.character(genotype)]
    idaa <- c(1, 0, 0)
    names(idaa) <- c(AA, Aa, aa)
    maa <- idaa[as.character(genotype)]
    list(oset = oset, ma = ma, maa = maa, types = c(AA, Aa, aa))
}
