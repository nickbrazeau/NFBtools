#' MIDmaker
#' no export
#' #' @noRd

MIDmaker <- function(length, homopolymercutoff){

  MID <- paste0(c(sample(x=c("C","G"), size=1), sample(x=c("C","G", "A", "T"), size=length-1, replace = T)), collapse="") # can't start with A or T

    while(max(unlist(parallel::mclapply(c("A", "C", "T", "G"), function(x){longestConsecutive(MID, x)}))) >= homopolymercutoff | # avoid homopolymers
          length(findPalindromes(BString(MID), min.armlength=2)) != 0 # avoid palindromes
          ){
      MID <- paste0(c(sample(x=c("C","G"), size=1), sample(x=c("C","G", "A", "T"), size=length-1, replace = T)), collapse="")
      }

  return(MID)

}




#' @title MID Primer Design
#'
#' @author Nick Brazeau
#' @description This will produces MIDs
#' @details This is a brute-force approach with recursion as the number of new MIDs proposed is based on "how close we are to solving the problem" with the current sit of MIDs
#'
#' @param design_fasta is expected to have the forward, reverse, and target gene in it. Must have `forward` and `reverse` in the contig names for the forward and reverse primers
#' @param MID2targetmismatchesAllowed is the number of base-pair mismatches allowed in the MID before it would bind to ANY place on the target gene


MIDPrimerFinder <- function(design_fasta = NULL,
                            MID2MIDmatchesAllowed = NULL,
                            MIDlength = NULL,
                            MIDhomopolymerallowance = NULL,
                            MID2targetmismatchesAllowed = NULL,
                            MIDnum = NULL
                            ){

  design_fasta <- Biostrings::readDNAStringSet(filepath = design_fasta, format = "fasta", use.names = T)

  ## Housekeeping
  fw <- design_fasta[ grepl("forward", tolower(names(design_fasta))) ]
  rv <- Biostrings::reverseComplement(design_fasta[ grepl("reverse", tolower(names(design_fasta))) ])
  tg <- design_fasta[ ! tolower(names(design_fasta)) %in% tolower(c(names(fw), names(rv))) ]

  # error handling

  hits <- lapply(list(fw, rv), function(x) {
    Biostrings::vmatchPattern(as.character(x), tg) # coercion here is fine since dnastrinset of length 1 ... becomes a named character vector
  })


  if( any(sapply(hits, function(x){ elementNROWS(x) }) == 0 ) ){
    stop("Your primers do not bind to your target")
  }

  if( any(sapply(hits, function(x){ elementNROWS(x) }) > 1 ) ){
    stop("Your primers have more than one binding site on your target")
  }

  if( all(sapply(c(MID2MIDmatchesAllowed, MIDlength, MIDhomopolymerallowance, MID2targetmismatchesAllowed, MIDnum), class) != "numeric") ){
    stop("Your MID specifications all must be of class numeric")
  }

  if( !all( c(MID2MIDmatchesAllowed, MIDlength, MIDhomopolymerallowance, MID2targetmismatchesAllowed, MIDnum) %in% sapply(c(MID2MIDmatchesAllowed, MIDlength, MIDhomopolymerallowance, MID2targetmismatchesAllowed, MIDnum), floor)) ){
    stop(" Your MID specifications must be whole numbers (no decimals)")
  }

  # housekeeping II
  # this is a very strange S4 object. I assume this will not always be backwards compatible...
  primerinfo <-
    tibble( primer = c("Forward", "Reverse"),
          start = c(
            hits[[1]]@ends[[1]] - hits[[1]]@width0 + 1,
            hits[[2]]@ends[[1]] - hits[[2]]@width0 + 1
            ),
          ends = c(
            hits[[1]]@ends[[1]],
            hits[[2]]@ends[[1]]
            ),
          width = c(
            hits[[1]]@width0,
            hits[[2]]@width0
            )
          )

  # START
  finalMIDs <- NULL # init N
  currMIDs <- NULL # init
  n <- MIDnum # init N
  int <- 1 # housekeeping


  while( length(finalMIDs) < n ){

    propMIDs <- replicate(n, MIDmaker(MIDlength, MIDhomopolymerallowance)) # init
    propMIDs <- c(propMIDs, currMIDs)


    while( any( TRUE == c( as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(propMIDs))) %in% as.vector(Biostrings::DNAStringSet(propMIDs))) ) | # catch any reverse complements -- repropose
           sum(sapply(propMIDs, function(x){Biostrings::vcountPattern(x, tg, max.mismatch = MID2targetmismatchesAllowed) })) > 0 | # catch any close matches to target
           any( TRUE %in%  duplicated(propMIDs) )
           ){
        # find culprits for rev comp
        revcount <- sum(Biostrings::reverseComplement(Biostrings::DNAStringSet(propMIDs)) %in% Biostrings::DNAStringSet(propMIDs))
        if(revcount > 0){
          propMIDs[ Biostrings::reverseComplement(Biostrings::DNAStringSet(propMIDs)) %in% Biostrings::DNAStringSet(propMIDs) ] <- replicate(revcount, MIDmaker(MIDlength, MIDhomopolymerallowance))
        }

        # find culprits for overlap targets
        ot_vect <- which( sapply(propMIDs, function(x){Biostrings::vcountPattern(x, tg, max.mismatch = MID2targetmismatchesAllowed) }) > 0 )
        if(length(ot_vect) > 0){
        propMIDs[ot_vect] <- replicate(length(ot_vect), MIDmaker(MIDlength, MIDhomopolymerallowance))
        }

        # find culprits for duplicates
        propMIDs[duplicated(propMIDs)] <-  replicate(sum(duplicated(propMIDs)), MIDmaker(MIDlength, MIDhomopolymerallowance))


        # make sure to return a vector
        propMIDs <- unlist(propMIDs) # if replicate has n=0, it returns a list class which will coerce our atomic vector into a list. Don't want this

      }





    if(int%%10 == 0){cat(paste0("Iteration  ", int, "\n ------- \n"))}
    # get pairs and pair.list
    pairs <- t(combn(propMIDs, m=2)) # choose 2
    pairs_list <- split(pairs, seq(nrow(pairs)))
    pairs_list_ret <- parallel::mclapply(pairs_list, function(x){stringdist::stringdistmatrix(x[1], x[2], method = c("hamming"))})
    pairs_list_ret <- tibble(primer = c(pairs[,1], pairs[,2]), Hdist = rep(unlist(pairs_list_ret),2)) # TODO sloppy

    MIDret <- pairs_list_ret %>%
      dplyr::group_by(primer) %>%
      dplyr::summarise(Hdist=min(Hdist))

    # update MID proposals
    if(min(MIDret$Hdist) <= MIDlength-MID2MIDmatchesAllowed){
      # Use distances as weights to fine tune the number of MIDs to keep
      keep <- floor( runif( n = 1, min = MIDnum*sum(MIDret$Hdist)/(nrow(MIDret)*MIDlength), max =  MIDnum - 1 ) ) # - 1 at end, so all propose
      # Standardize to length and then select those MIDs that have the best distances so far
      currMIDs <- sample(MIDret$primer, size = keep, prob = MIDret$Hdist/MIDlength)
      # update N for propsoal
      n <- MIDnum - length(currMIDs)

      # housekeeping
      int <- int + 1

    } else {

      finalMIDs <- MIDret$primer

    }

  }

  ret <- list(finalMIDs = finalMIDs,
              primerinfo = primerinfo)
  return(ret)


}
