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


bedtoolsgenomecov2bamCov <- function(gencovdir = NULL, lvls = c(1,5,10,25,50,75,100), bpwindow=1e4){

  require(tidyverse)
  require(zoo)
  require(parallel)

  gencovfiles <- dir(gencovdir, full.names = T)

  retlist_all <- lapply(gencovfiles, function(file){
    dat <- read_tsv(file, col_names = F)
    colnames(dat) <- c("CHROM", "POS", "depth")
    datchromlist <- split(x = dat, f=factor(dat$CHROM))
    ttlbp <- sum(unlist(lapply(datchromlist, function(x) return(max(x$POS)))))


    ## genome coordinates
    genomcoords <- dat %>%
      dplyr::group_by(CHROM) %>%
      dplyr::summarise(smpl = gsub(".long.cov", "", basename(file)), start=min(POS), end=max(POS))
    genomcoords <- rbind.data.frame(genomcoords, c("genome",
                                                   gsub(".long.cov", "", basename(file)),
                                                   1,
                                                   sum(unlist(lapply(datchromlist, function(x) return(max(x$POS)))))))


    # window cov
    windowfunctioncalculator <- function(df){
    start <- seq(1, max(df$POS), bpwindow)
    end = lead(start)
    end[is.na(end)] <- max(df$POS)
    windowcov <- data.frame(start=start, end=end, meancov=NA)
    for(i in 1:nrow(windowcov)){

      windowcov$meancov[i] <- mean( df$depth[ windowcov$start[i]:windowcov$end[i] ] ) # dat has every line indexed, so this is OK

    }
    CHROM <- levels(factor(df$CHROM))
    windowcov <- cbind.data.frame(CHROM=CHROM, windowcov)
      return(windowcov)
  }

  windowcov <- parallel::mclapply(datchromlist, windowfunctioncalculator)
  windowcov <- do.call("rbind.data.frame", windowcov)

    ## genome summary
    genomsummary <- dat %>%
      dplyr::select(POS) %>%
      dplyr::summarise(smpl = gsub(".long.cov", "", basename(file)),
                       n=n(),
                       min=min(POS),
                       q25 = quantile(POS, prob=0.25),
                       median = median(POS),
                       mean = mean(POS),
                       q75 = quantile(POS, prob=0.75),
                       max = max(POS)
      )

    ## Percentage of Genome Covered by levels
    genomcovperc <- sapply(lvls, calcCovperc, depth=dat$depth, ttlbp=ttlbp)
    genomcovperc <- data.frame(smpl = gsub(".long.cov", "", basename(file)),
                               covdepth = names(genomcovperc),
                               val = genomcovperc)
    genomcovperc$covdepth <- factor(genomcovperc$covdepth, levels=paste0(lvls,"x"))

    ## chrom coverage summary
    chromcovsummary <- dat %>%
      group_by(CHROM) %>%
      dplyr::select(CHROM, POS) %>%
      dplyr::summarise(smpl = gsub(".long.cov", "", basename(file)),
                       n=n(),
                       min=min(POS),
                       q25 = quantile(POS, prob=0.25),
                       median = median(POS),
                       mean = mean(POS),
                       q75 = quantile(POS, prob=0.75),
                       max = max(POS)
      )


    ## Percentage of Chromosome Covered by levels
    chromlistcovperc <- parallel::mclapply(datchromlist, function(df){
      sapply(lvls, calcCovperc, depth=df$depth, ttlbp=ttlbp)

    }
    )

    chromlistcovperc <- do.call("rbind", chromlistcovperc)
    chromlistcovperc <- cbind.data.frame(smpl = gsub(".long.cov", "", basename(file)),
                                         CHROM = row.names(chromlistcovperc), chromlistcovperc)



    retlist <- list(windowcov = windowcov,
                    genomcoords = genomcoords,
                    genomsummary = genomsummary,
                    genomcovperc = genomcovperc,
                    chromcovsummary = chromcovsummary,
                    chromlistcovperc = chromlistcovperc)


    return(retlist)
    # end of function for lapply loop
  }
  # end of lapply for files
  )
  # end of function overall
  class(retlist_all) <- "bamCov"
  names(retlist_all) <- gsub(".long.cov", "", basename(dir(gencovdir, full.names = T))) # add names
  return(retlist_all)
}





#' @title bamCov2OverallPercCov
#'
#'
#' @export
#'

bamCov2OverallPercCov <- function(input = NULL){

  require(tidyverse)
  require(RColorBrewer)

  stop(class(input) != "bamCov", "Input must be of class bamCov. See the bedtoolsgenomecov2bamCov function.")


  dat <- do.call("rbind.data.frame", lapply(input, function(x){x[["genomcovperc"]]}))

  # Thanks to Jelena-bioinf & Megatron for this cool trick -- https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  n <- length(levels((factor(dat$covdepth))))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))), n, replace = F)

  plotObj <- dat %>%
    dplyr::mutate(covdepth_f = factor(covdepth)) %>%
    ggplot() +
    geom_dotplot(aes(x=smpl, y=val, fill=covdepth_f), binaxis='y', stackdir='center') +
    geom_vline(aes(xintercept = as.numeric(smpl)), linetype="dashed", colour = "#bdbdbd") +
    scale_fill_manual("x-Fold \nCoverage", values = col_vector) +
    ggtitle("Genomic Coverage by Sample") +
    xlab("Samples") + ylab("Proportion of Genome Covered") +
    theme(plot.title = element_text(hjust = 0.5, family="Arial", size=17, face = "bold"),
          axis.text.x = element_text(family="Arial", size=13, face = "bold", angle=90, vjust=0.5),
          axis.title.x = element_text(family="Arial", size=15, face = "bold"),
          axis.text.y = element_text(family="Arial", size=13, face = "bold", hjust=0.5),
          axis.title.y = element_text(family="Arial", size=15, face = "bold"),
          legend.title = element_text(family="Arial", size=12, face = "bold"),
          legend.text = element_text(family="Arial", size=10, face = "bold"),
          legend.title.align = 0.5,
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size=2),
          axis.ticks = element_blank()
          )


  return(plotObj)

}






#' @title bamCov2SmplGenomCov
#'
#'
#' @export
#'
#'



bamCov2SmplGenomCov <- function(input = NULL, window = 1e4){
  stop(class(input) != "bamCov", "Input must be of class bamCov")

  sapply(input, "[[", 2)


  dat <- do.call("rbind.data.frame", )



}


