#' inherits bamcovobj
#' @noRd

is.bamCovObj <- function(x){inherits(x, "bamCovObj")}

#' calcCovperc
#' @noRd


calcCovperc <- function(depth, ttlbp, lvl){
  ret <- sum(depth >= lvl)/ttlbp
  names(ret) <- paste0(lvl, "x")
  return(ret)
}


#' @title bedtoolsgenomecov2bamCov
#'
#'
#' @export


bedtoolsgenomecov2bamCov <- function(gencovdir = NULL, lvls = c(1,5,10,25,50,75,100), bpwindow=1e4){

  gencovfiles <- dir(gencovdir, full.names = T)

  retlist_all <- parallel::mclapply(gencovfiles, function(file){
    dat <- readr::read_tsv(file, col_names = F)
    colnames(dat) <- c("CHROM", "POS", "depth")
    datchromlist <- split(x = dat, f=factor(dat$CHROM))
    ttlbp <- sum(unlist(lapply(datchromlist, function(x) return(max(x$POS)))))
    if(ttlbp != nrow(dat)){
      stop("Total Base-pairs calculated does not equal rows from Bedtools GenomeCov. Error in how file was generated from Bedtools.")
    }


    ## genome coordinates
    genomcoords <- dat %>%
      dplyr::group_by(CHROM) %>%
      dplyr::summarise(smpl = gsub(".long.cov", "", basename(file)), start=min(POS), end=max(POS))
    genomcoords <- rbind.data.frame(genomcoords, data.frame(CHROM = "genome",
                                                   smpl = gsub(".long.cov", "", basename(file)),
                                                   start = 1,
                                                   end = sum(unlist(lapply(datchromlist, function(x) return(max(x$POS)))))))


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
    windowcov <- cbind.data.frame(smpl = gsub(".long.cov", "", basename(file)), CHROM=CHROM, windowcov)
      return(windowcov)
  }

  windowcov <- parallel::mclapply(datchromlist, windowfunctioncalculator)
  windowcov <- do.call("rbind.data.frame", windowcov)

    ## genome summary
    genomsummarydepth <- dat %>%
      dplyr::select(depth) %>%
      dplyr::summarise(smpl = gsub(".long.cov", "", basename(file)),
                       n=n(),
                       min=min(depth),
                       q25 = quantile(depth, prob=0.25),
                       median = median(depth),
                       mean = mean(depth),
                       q75 = quantile(depth, prob=0.75),
                       max = max(depth)
      )

    ## Percentage of Genome Covered by levels
    genomcovperc <- sapply(lvls, calcCovperc, depth=dat$depth, ttlbp=ttlbp)
    genomcovperc <- data.frame(smpl = gsub(".long.cov", "", basename(file)),
                               covdepth = factor(names(genomcovperc), levels=paste0(lvls,"x")),
                               val = genomcovperc)

    ## chrom coverage summary
    chromcovdepthsummary <- dat %>%
      group_by(CHROM) %>%
      dplyr::select(CHROM, depth) %>%
      dplyr::summarise(smpl = gsub(".long.cov", "", basename(file)),
                       n=n(),
                       min=min(depth),
                       q25 = quantile(depth, prob=0.25),
                       median = median(depth),
                       mean = mean(depth),
                       q75 = quantile(depth, prob=0.75),
                       max = max(depth)
      )

    ## Percentage of Chromosome Covered by levels
    chromlistcovperc <- parallel::mclapply(datchromlist,
                                           function(df){
                                             chrombp = nrow(df)
                                             ret <- sapply(lvls, calcCovperc, depth=df$depth, ttlbp=chrombp)
                                             return(ret)

    }
    )

    chromlistcovperc <- do.call("rbind", chromlistcovperc)
    chromlistcovperc <- cbind.data.frame(smpl = gsub(".long.cov", "", basename(file)),
                                         CHROM = row.names(chromlistcovperc), chromlistcovperc)

    chromlistcovperc <- chromlistcovperc %>%
      tidyr::gather(key="covdepth", value="val", 3:ncol(chromlistcovperc)) %>%
      dplyr::mutate(covdepth=factor(covdepth, levels= paste0(lvls, "x")))





    retlist <- list(windowcov = windowcov,
                    genomcoords = genomcoords,
                    genomsummarydepth = genomsummarydepth,
                    genomcovperc = genomcovperc,
                    chromcovsummarydepth = chromcovdepthsummary,
                    chromlistcovperc = chromlistcovperc)


    # end of function for lapply loop
  }
  # end of lapply for files
  )
  # end of function overall
  class(retlist_all) <- append(class(retlist_all), "bamCovObj")
  names(retlist_all) <- gsub(".long.cov", "", basename(dir(gencovdir, full.names = T))) # add names
  return(retlist_all)
}




#' @title bamCov2OverallPercCov
#'
#'
#' @export
#'

bamCov2OverallPercCov <- function(input = NULL){

  if(!is.bamCovObj(input)){
    stop("Input must be of class bamCovObj. See the bedtoolsgenomecov2bamCov function.")
  }

  # extract
  dat <- do.call("rbind.data.frame", lapply(input, function(x){x[["genomcovperc"]]}))

  # tidy & plot
  plotObj <- dat %>%
    ggplot() +
    geom_dotplot(aes(x=smpl, y=val, fill=covdepth), binaxis='y', stackdir='center', alpha=0.5) +
    geom_vline(aes(xintercept = as.numeric(smpl)), linetype="dashed", colour = "#bdbdbd") +
    scale_fill_manual("x-Fold \nCoverage", values = brewer.pal(n = length(levels((factor(dat$covdepth)))), name = "Set2")) +
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






#' @title bamCov2SmplPercCov
#'
#'
#' @export
#'
#'


bamCov2SmplPercCov <- function(chromlistcovpercdf = NULL){

  # tidy and plot

  plotObj <- chromlistcovpercdf %>%
    ggplot() +
    geom_dotplot(aes(x=CHROM, y=val, fill=covdepth), binaxis='y', stackdir='center', alpha=0.5) +
    geom_vline(aes(xintercept = as.numeric(CHROM)), linetype="dashed", colour = "#bdbdbd") +
    scale_fill_manual("x-Fold \nCoverage", values = brewer.pal(n = length(levels((factor(chromlistcovpercdf$covdepth)))), name = "Set2")) +
    ggtitle(paste("Coverage by Chromosome for Sample", chromlistcovpercdf$smpl[1])) +
    xlab("Chromosome") + ylab("Proportion of Chromosome Covered") +
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





#' @title bamCov2SmplChromCov
#'
#'
#' @export
#'
#'


bamCov2SmplChromCov <- function(genomcoordsdf = NULL, windowcovdf = NULL){
  genomcoordsdf <- genomcoordsdf %>%
    dplyr::filter(CHROM != "genome") %>%
    dplyr::arrange(end) %>%
    dplyr::mutate(CHROM_num = seq(1:length(CHROM)))  # hacky on purpose

  chromnumdf <- genomcoordsdf %>%
    dplyr::select(CHROM, CHROM_num)

  windowcovdf <- windowcovdf %>%
    dplyr::left_join(x=., y=chromnumdf, by="CHROM")

  plotObj <- ggplot() +
    geom_rect(data=genomcoordsdf, aes(xmin=start, xmax=end, ymin=CHROM_num-0.25, ymax=CHROM_num+0.25), fill="#d9d9d9", color="#d9d9d9") +
    geom_rect(data=windowcovdf, aes(xmin=start, xmax=end, ymin=CHROM_num-0.2, ymax=CHROM_num+0.2, fill=meancov)) +
    scale_fill_gradientn("Mean \n Coverage", colours = c("#313695", "#ffffbf", "#d53e4f")) +
    scale_y_continuous(breaks = 1:nrow(genomcoordsdf), labels=chromnumdf$CHROM) +
    ggtitle(paste("Mean Coverage with a ", windowcovdf$end[1] - windowcovdf$start[1], " base-pair Sliding Window \n by Chromosome for Sample:", windowcovdf$smpl[1])) +
    xlab("Chromosome Position") + ylab("Chromosome") +
    theme(plot.title = element_text(hjust = 0.5, family="Arial", size=17, face = "bold"),
          axis.text.x = element_text(family="Arial", size=13, face = "bold", angle=90, vjust=0.5),
          axis.title.x = element_text(family="Arial", size=15, face = "bold"),
          axis.text.y = element_text(family="Arial", size=13, face = "bold", hjust=0.5),
          axis.title.y = element_text(family="Arial", size=15, face = "bold"),
          legend.title = element_text(family="Arial", size=12, face = "bold"),
          legend.text = element_text(family="Arial", size=10, face = "bold"),
          legend.title.align = 0.5,
          legend.box.just = "center",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size=2),
          axis.ticks = element_blank()
    )

  return(plotObj)

}



#' @title bamCov2SmplRaster
#'
#'
#' @export
#'
#'

bamCov2SmplRaster <- function(input=NULL){

  if(!is.bamCovObj(input)){
    stop("Input must be of class bamCovObj. See the bedtoolsgenomecov2bamCov function.")
  }

  chromcovsummarydepthlist <- lapply(input, function(x){x[["chromcovsummarydepth"]]})
  chromlistcovperclist <- lapply(input, function(x){x[["chromlistcovperc"]]})
  genomcoordslist <- lapply(input, function(x){x[["genomcoords"]]})
  windowcovlist <- lapply(input, function(x){x[["windowcov"]]})

  if(FALSE %in% c(names(windowcovlist) == names(chromlistcovperclist) & names(chromlistcovperclist) == names(chromcovsummarydepthlist))){
    stop("There is an error in the input file. See the bedtoolsgenomecov2bamCov function.")
  }

  smpls <- names(windowcovlist)

  SmplPercCovlist <- lapply(chromlistcovperclist, bamCov2SmplPercCov)
  SmplChromCovlist <- mapply(bamCov2SmplChromCov, genomcoordsdf = genomcoordslist, windowcovdf = windowcovlist, SIMPLIFY = F)



# relist and make each sample its own grob of 3 plot objs
  grobsrelist <- lapply(smpls, function(smpl){
    # plot 1
      p1 <- SmplPercCovlist[[smpl]]
    # summary table (i.e. plot 2)        #https://stackoverflow.com/questions/11774703/adding-text-to-a-grid-table-plot/11775211#11775211
      p2 <- gridExtra::tableGrob(chromcovsummarydepthlist[[smpl]], rows=NULL)
      title <- grid::textGrob("Summary of Coverage Depth by Chromosome", gp=grid::gpar(fontfamily="Arial", fontsize=16, fontface="bold"))
      padding <- unit(5,"mm")
      p2 <- gtable::gtable_add_rows(
        p2,
        heights = grid::grobHeight(title) + padding,
        pos = 0)
      p2 <- gtable::gtable_add_grob(
        p2,
        list(title),
        1, 1, 1, ncol(p2))
    # plot 3
      p3 <- SmplChromCovlist[[smpl]]
    # plot 4 white space - https://stackoverflow.com/questions/21529926/arrange-ggplots-together-in-custom-ratios-and-spacing
    blank <- grid::rectGrob(gp=grid::gpar(col="white")) # make a white spacer grob
    ret <- list(p1,
                p2,
                p3,
                blank)
    return(ret)

  })


# plot the grobs
  lapply(grobsrelist, function(grobs){
    grid::grid.newpage()
    gridExtra::grid.arrange(grobs = grobs, layout_matrix = rbind(c(rep(1,10), NA, rep(3,10)),
                                                                 c(rep(1,10), NA, rep(3,10)),
                                                                 c(rep(1,10), NA, rep(3,10)),
                                                                 c(rep(2,10), NA, rep(3,10)),
                                                                 c(rep(2,10), NA, rep(3,10)),
                                                                 c(rep(2,10), NA, rep(3,10)),
                                                                 c(rep(4,21))

    )
    )

  })



}










