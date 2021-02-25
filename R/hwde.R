#' Fit relevant models, and test for various types of departure from Hardy-Weinberg equilibrium. Allows only 2 alleles per locus. The number of loci is arbitrary.
#'
#' Fits models for genotypic disequilibria, as described in Huttley and Wilson (2000), Weir (1996) and Weir and Wilson (1986).
#' 
#' @param data Must have a column of frequencies, by default called \code{Observed}, and one or more columns giving genotype information, in the form AA, Aa, aa. or e.g., MM, MN, NN. (The choice of letters is arbitrary.) Additionally, there may be a column that gives information on groupings, by default called \code{Population}.
#' @param gp Gives the name of the column, if any, that has information on groups within the data.
#' @param termlist Use to specify a user-defined sequence of models. See the vignette \bold{hwde.pdf} or \bold{hwde.html}
#' @param refmodel For each model in \code{termlist}, specifies a reference model, which will be updated to include the additional terms.
#' @param loci Gives name(s) of columns that hold information on genotypes. By default, these are taken to be \code{locus1}, \code{locus2}, etc.
#' @param observed Name (by default \code{Observed}) of the column that holds the frequenceies.
#' @param keep.models Should a list be returned that holds the full sequence of models that were fitted?
#' @param aovtable.print Should the anova table be printed?
#' @param group.terms Should model terms be grouped according to hierarchy, for the anova table?
#' @param allele.chars A sequence of letters used to code for the loci. By default \emph{a}, \emph{b}, \emph{c}, ... are used
#' 
#' @details See the document \bold{hwde.pdf} or \bold{hwde.html} for details. See the references (below) for information on the interpretation of model parameters.
#' 
#' @return 
#' \itemize{
#'   \item{anovatab}{anova (analysis of deviance) table}
#'   \item{data.df}{Data, and contrasts used in fitting the various models.}
#'   \item{aovtab.terms}{This string holds, for each model that is fitted. the terms that have appeared in the model formula. The text strings for the distinct models are concatenated.}
#'   \item{models}{Optionally, this holds the complete sequence of qmodel objects that were fitted}
#' }
#' 
#' @references 
#' \itemize{
#'   \item Huttley, G.A. and Wilson, S.R. 2000. Testing for concordant equilibria between population samples. \emph{Genetics} \bold{156}:2127-2135.
#'   \item Weir, B.S. 1996.  \emph{Genetic Data Analysis II.} Sinauer.
#'   \item 3. Weir, B.S. and Wilson, S.R. 1986.  Log-linear models for linked loci. \emph{Biometrics} \bold{42}:665-670.
#' }
#' 
#' @author J.H. Maindonald
#' 
#' @seealso \code{\link{make.contrasts}}, \code{\link{decode.genotypes}}
#' 
#' @examples 
#' data(IndianIrish)
#' hwde(data=IndianIrish)
#' data(mendelABC)
#' hwde(data=mendelABC, loci=c("seedshape", "cotylcolor", "coatcolor"))
#' 
#' @export
#' @importFrom stats formula glm poisson update

hwde <-
function (data = hwde::IndianIrish, gp = "Population", termlist = NULL,
    refmodel = NULL, loci = paste("locus", 1:(dim(data)[2] -
        1), sep = ""), observed = "Observed", keep.models = FALSE,
    aovtable.print = TRUE, group.terms = TRUE, allele.chars = letters)
{
   tmaker <- function(trms = c("ma", "mb"), currbase = 1, group.terms = TRUE,
        pref = "") {
        if (length(trms) == 1) {
            refmodel <- currbase
            plus <- trms
        }
        else if (group.terms) {
            plus <- paste("(", paste(trms, collapse = "+"), ")",
                sep = "")
            refmodel <- currbase
        }
        else {
            plus <- paste(pref, c(trms, paste("(", paste(trms,
                collapse = "+"), ")", sep = "")), sep = "")
            refmodel <- rep(currbase, length(plus))
        }
        list(plus, refmodel)
    }
    nloci <- length(loci)
    colnames <- names(data)
    for (i in 1:length(loci)) if (!(loci[i] %in% colnames)) {
        nloci <- i - 1
        loci <- loci[1:nloci]
        break
    }
    obs <- data[, observed]
    nobs <- length(obs)
    gploc <- loci
    contrasts.info <- make.contrasts(data = data[, loci], allele.chars = allele.chars)
    contr.df <- contrasts.info$contrasts.df
    oset <- contrasts.info$oset
    list.columns <- contrasts.info$list.columns
    if (gp %in% colnames)
        data.df <- cbind(data.frame(obs = obs, gp = data[, gp]),
            contr.df, oset=oset)
    else data.df <- cbind(data.frame(obs = obs), contr.df,
                          oset=oset)
    gpslot <- as.numeric("gp" %in% names(data.df))
    if (is.null(termlist)) {
        addterms <- NULL
        newbase <- 1
        refmodel <- NULL
        for (colvec in list.columns) {
            trms <- tmaker(colvec, currbase = newbase, group.terms = group.terms)
            addterms <- c(addterms, trms[[1]])
            newbase <- newbase + length(trms[[1]])
            refmodel <- c(refmodel, trms[[2]])
        }
        if (gpslot) {
            addtermsg <- paste("gp:", addterms, sep = "")
            addterms <- c("gp", addterms, addtermsg)
            refmodel <- c(1, refmodel + 1, refmodel + newbase)
        }
        addterms <- paste("+", addterms, sep = "")
    }
    else {
        addterms <- termlist
    }
    form1 <- formula(paste("obs", "~", "1"))
    m1 <- glm(form1, family = poisson, offset = log(oset), data = data.df)
    n <- length(addterms)
    firstchar <- rep(" ", n)
    currmod <- m1
    mlist <- list(m1)
    for (i in 1:n) {
        modi <- paste("m", i + 1, sep = "")
        updi <- formula(paste(".~.", addterms[i], sep = ""))
        newmod <- update(currmod, updi)
        if ((i < n) & (refmodel[i + 1] > refmodel[i])) {
            if (!group.terms)
                firstchar[i] <- "r"
            else firstchar <- " "
            currmod <- newmod
        }
        assign(paste("m", i + 1, sep = ""), newmod)
        mlist <- c(mlist, list(newmod))
    }
    anovatab <- do.call("anova", mlist)
    nam <- c("m", paste(firstchar, addterms, sep = ""))
    if (!group.terms)
        nam[1] <- paste("r", nam[1], sep = "")
    anovatab[2:(n + 1), "Deviance"] <- anovatab[refmodel, "Resid. Dev"] -
        anovatab[2:(n + 1), "Resid. Dev"]
    anovatab[2:(n + 1), "Df"] <- anovatab[refmodel, "Resid. Df"] -
        anovatab[2:(n + 1), "Resid. Df"]
    aovtab.terms <- attributes(anovatab)$heading[2]
    row.names(anovatab) <- c(refmodel[1], addterms)
    attributes(anovatab)$heading <- NULL
    if (aovtable.print == TRUE) {
        print("Analysis of Deviance Table")
        print(anovatab, rowlab = nam)
    }
    if (!keep.models)
        invisible(list(anovatab = anovatab, data.df = data.df,
            aovtab.terms = aovtab.terms))
    else invisible(list(anovatab = anovatab, data.df = data.df,
        aovtab.terms = aovtab.terms, models = mlist))
}
