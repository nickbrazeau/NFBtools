#' @title calcCovperc
#'
#'
#' not exported

calcCovperc <- function(depth, ttlbp, lvl){
  ret <- sum(depth > lvl)/ttlbp
  names(ret) <- paste0(lvl, "x")
  return(ret)
}


#' @title bedtoolsgenomecov2bamCov
#'
#'
#' @export

require(tidyverse)

function(gencovdir = NULL, lvls = c(1,5,seq(10,150, by=10))){

  gencovfiles <- dir(gencovdir, full.names = T)

  retlist <- lapply(gencovfiles, function(file){
    dat <- read_tsv(file, col_names = F)
    colnames(dat) <- c("CHROM", "POS", "depth")
    datchromlist <- split(x = dat, f=factor(dat$CHROM))
    ttlbp <- sum(unlist(lapply(datchromlist, function(x) return(max(x$POS)))))

    ## genome coverage summary
    genomecov <- summary(dat)

    ## Percentage of Genome Covered by levels
    genomcovperc <- sapply(lvls, calcCovperc, depth=dat$depth, ttlbp=ttlbp)
    genomcovperc <- data.frame(covdepth = names(genomcovperc),
                               val = genomcovperc)

    ## chrom coverage summary
    chromlistcov <- sapply(datchromlist, function(x){
      chromlistcov <- summary(x)
        }
      )


    ## Percentage of Chromosome Covered by levels
    chromlistcovperc <- lapply(datchromlist, function(df){
      sapply(lvls, calcCovperc, depth=df$depth, ttlbp=ttlbp)

      }
    )

    chromlistcovperc <- do.call("rbind", chromlistcovperc)
    chromlistcovperc <- cbind.data.frame(smpl = gsub("long.cov", "", basename(file)),
                                         CHROM = row.names(chromlistcovperc), chromlistcovperc)



    retlist <- list(genomecov = cbind.data.frame(smpl = gsub("long.cov", "", basename(file)), genomecov),
                    genomcovperc = genomcovperc,
                    chromlistcov = cbind.data.frame(smpl = gsub("long.cov", "", basename(file)), chromlistcov),
                    chromlistcovperc = chromlistcovperc)
    return(retlist)
    # end of function for lapply loop
    }
  # end of lapply for files
  )
# end of function overall
}
