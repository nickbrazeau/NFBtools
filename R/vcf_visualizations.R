
#------------------------------------------------
#' @title Plot vcf outputs collectively
#'
#' @description Takes a \code{vcf} or \code{vcfRobject} and plots the GT levels
#'
#' @param
#'
#' @export

vcfR2GTCov <- function(vcfRobject=NULL, biallelic=T){
  require(vcfR)
  require(tidyverse)
  # plot vcf chrom as a list object
  if(is.null(vcfRobject)){
    stop("Must provide an input file")
  }
  if(! class(vcfRobject) %in% "vcfR"){
    stop("vcfRobject must contain  class vcfR")
  }

  vcf <- vcfRobject # legacy

  if( !identical(vcfRobject, vcfR::extract.indels(vcf[vcfR::is.biallelic(vcf)], return.indels = F)) ){
    stop("Your file does not appear to be a vcf with only biallelic SNPs")
  }


  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  ploidy <- stringr::str_extract(vcf@meta[grepl("ploidy", vcf@meta)], "ploidy=[0-9]")
  ploidy <- stringr::str_split(t(ploidy), "=", simplify = T)
  ploidy <- as.numeric(ploidy[1,2])

  if(ploidy == 2){

    snpmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
    snpmatrix[snpmatrix == "0/0"] <- 0
    snpmatrix[snpmatrix == "0/1"] <- 1
    snpmatrix[snpmatrix == "1/1"] <- 2
    snpmatrix <- apply(snpmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
  } else if(ploidy == 1){
    snpmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=T)
    snpmatrix[snpmatrix == 0] <- 0
    snpmatrix[snpmatrix == 1] <- 2
  } else {
    stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by polyIBD")
  }

  CHROM <- vcfR::getCHROM(vcf)
  POS <- vcfR::getPOS(vcf)
  snpdf <- cbind.data.frame(CHROM, POS, snpmatrix)
  snpdf <- snpdf %>%
    tidyr::gather(key="sample", value="GT", 3:ncol(snpdf)) %>%
    dplyr::mutate(GTfact = factor(GT, levels=c(0,1,2), labels=c("Ref", "Het", "Alt")))

  snpdflist <- split(snpdf, f=snpdf$CHROM)

  lagPOS <- function(df){
    df <- df %>%
      dplyr::mutate(start = lag(POS)) %>%
      dplyr::mutate(end = POS) %>%
      dplyr::mutate(sample = factor(sample)) %>%
      dplyr::mutate(sampleindex = as.numeric(factor(sample)))  # Will pull out levels
    return(df)
  }

  snpdflist <- lapply(snpdflist, lagPOS)

  snpdf <- do.call("rbind", snpdflist)
  plotobj <- snpdf %>%
    ggplot() +
    geom_rect(mapping=aes(xmin=start, xmax=end, ymin=sampleindex-0.49, ymax=sampleindex+0.49, fill=GTfact)) +
    scale_y_continuous(breaks = min(snpdf$sampleindex):max(snpdf$sampleindex),
                       labels = levels(factor(snpdf$sample))) +
    scale_fill_manual("Genotype Call", values=c("#AA0A3C", "#313695", "#fee090", "#cccccc")) +
    facet_grid(CHROM~.) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size=9, family = "Arial", angle = 45),
          axis.title.y = element_text(size=14, face="bold", family = "Arial"),
          axis.text.y = element_text(size=12, face="bold", family = "Arial")
          )

  retlist = list(snpdflist = snpdflist,
                 plot = plotobj)
  return(retlist)

}
