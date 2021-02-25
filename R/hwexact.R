#' Exact test of HWE
#' 
#' This function implements an exact SNP test of Hardy-Weinberg equilibrium.
#' 
#' @param obs.hom1 Observed number of individuals homozygous for the common allele.
#' @param obs.hets Observed number of heterozygous individuals.
#' @param obs.hom2 Observed number of individuals homozygous for the rare allele.
#' @param data If not \code{NULL}, specifies a data frame holding the data.
#' @param loci If \code{data} is not \code{NULL}, column that holds some combination of symbols of the form \code{AA}, \code{Aa} and \code{aa}.
#' @param observed Column that holds the frequencies, if \code{data} is not \code{NULL}.
#' 
#' @return The exact p-value testing departure from Hardy-Weinberg Equilibrium, conditional on the observed relative numbers of alleles of the two types.
#' 
#' @references Wigginton, J.E.; Cutler, D.J. & Abecasis, G.R.; A Note on exact tests of Hardy-Weinberg equilibrium; American Journal of Human
#' 
#' @author Randall Johnson -- adapted from code by Wigginton et al Genetics, 2005, 76, 887-93
#' 
#' @examples 
#' hwexact(68, 28, 4)
#' 
#' data(IndianIrish)
#' hwexact(data = IndianIrish, loci = 'locus1', observed = 'Observed')
#' 
#' @export

hwexact <- function(obs.hom1=68, obs.hets=28, obs.hom2=4, data=NULL, loci="locus1", observed="Observed")
{
  if(is.data.frame(obs.hom1))data=obs.hom1
  if(!is.null(data)){
    loc <- data[, loci]
    types <- decode.genotypes(loc)$types
    obs.hom1 <- data[loc==types[1], observed]
    obs.hom2 <- data[loc==types[3], observed]
    obs.hets <- data[loc==types[2], observed]
  }
  if(length(obs.hom1) != length(obs.hets) | length(obs.hets) != length(obs.hom2))
    stop('length of obs.hom1, obs.hets, and obs.hom2 must be the same')

  return(unlist(lapply(1:length(obs.hom1), function(i){snphwe(obs.hom1[i], obs.hets[i], obs.hom2[i])})))
}



snphwe <- function(obs.hom1, obs.hets, obs.hom2)
{
  if(obs.hom1 < 0 || obs.hom2 < 0 || obs.hets < 0)
    stop('Negative count(s) are not valid')
  
  # total number of genotypes
  N <- obs.hom1 + obs.hom2 + obs.hets
  
  # rare homozygotes, common homozygotes
  obs.homr <- min(obs.hom1, obs.hom2)
  obs.homc <- max(obs.hom1, obs.hom2)
  
  # number of rare allele copies
  rare  <- obs.homr * 2 + obs.hets
  
  # Initialize probability array
  probs <- rep(0, 1 + rare)
  
  # Find midpoint of the distribution
  mid <- floor(rare * ( 2 * N - rare) / (2 * N))
  if( (mid %% 2) != (rare %% 2) )
    mid <- mid + 1
  
  probs[mid + 1] <- 1.0
  mysum <- 1.0
  
  # Calculate probablities from midpoint down 
  curr.hets <- mid
  curr.homr <- (rare - mid) / 2
  curr.homc <- N - curr.hets - curr.homr
  
  while( curr.hets >=  2)
  {
    probs[curr.hets - 1]  <- probs[curr.hets + 1] * curr.hets * (curr.hets - 1.0) / (4.0 * (curr.homr + 1.0)  * (curr.homc + 1.0))
    mysum <- mysum + probs[curr.hets - 1]
    
    # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
    curr.hets <- curr.hets - 2
    curr.homr <- curr.homr + 1
    curr.homc <- curr.homc + 1
  }    
  
  # Calculate probabilities from midpoint up
  curr.hets <- mid
  curr.homr <- (rare - mid) / 2
  curr.homc <- N - curr.hets - curr.homr
  
  while( curr.hets <= rare - 2)
  {
    probs[curr.hets + 3] <- probs[curr.hets + 1] * 4.0 * curr.homr * curr.homc / ((curr.hets + 2.0) * (curr.hets + 1.0))
    mysum <- mysum + probs[curr.hets + 3]
    
    # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
    curr.hets <- curr.hets + 2
    curr.homr <- curr.homr - 1
    curr.homc <- curr.homc - 1
  }    
  
  # P-value calculation
  target <- probs[obs.hets + 1]
  p <- min(1.0, sum(probs[probs <= target])/ mysum)
  
  #plo <- min(1.0, sum(probs[1:obs.hets + 1]) / mysum)
  
  #phi <- min(1.0, sum(probs[obs.hets + 1: rare + 1]) / mysum)
  
  return(p)
}
