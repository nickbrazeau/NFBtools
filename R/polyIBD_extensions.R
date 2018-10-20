

#------------------------------------------------
#' @title Plot polyIBD sim as gtcov mat
#' @details this is temporary
#'
#' @export

polyIBDsim2GTcov <- function(polyIBDsim){

  haplotype1 <- data.frame(polyIBDsim$SampleHaplotypes$Sample1)
  colnames(haplotype1) <- paste0("Smpl1_", colnames(haplotype1))
  haplotype2 <- data.frame(polyIBDsim$SampleHaplotypes$Sample2)
  colnames(haplotype2) <- paste0("Smpl2_", colnames(haplotype2))

  IBDws_sample1 <- data.frame(polyIBDsim$IBDws$Sample1)
  IBDws_sample1[IBDws_sample1 == 0] <- "U"
  IBDws_sample1[IBDws_sample1 == 1] <- "I"

  IBDws_sample2 <- data.frame(polyIBDsim$IBDws$Sample2)
  IBDws_sample2[IBDws_sample2 == 0] <- "U"
  IBDws_sample2[IBDws_sample2 == 1] <- "I"

  IBDbs <- data.frame(polyIBDsim$IBDbs)
  IBDbs[IBDbs == 0] <- "U"
  IBDbs[IBDbs == 1] <- "I"

  snpdf <- cbind.data.frame(CHROM = polyIBDsim$CHROMPOS$CHROM,
                            POS = polyIBDsim$CHROMPOS$POS,
                            polyIBDsim$gtmatrix,
                            haplotype1,
                            haplotype2,
                            IBDws_sample1,
                            IBDws_sample2,
                            IBDbs)
  snpdf <- snpdf %>%
    tibble::as.tibble(.) %>%
    tidyr::gather(key="param", value="GT", 3:ncol(snpdf)) %>%
    dplyr::mutate(GTfact = factor(GT, levels=c(0,1,2, "U", "I"), labels=c("Ref", "Het", "Alt", "Not IBD", "IBD")))

  snpdflist <- split(snpdf, f=snpdf$CHROM)

  lagPOS <- function(df){
    df <- df %>%
      dplyr::mutate(start = lag(POS)) %>%
      dplyr::mutate(end = POS) %>%
      dplyr::mutate(param = factor(param)) %>%
      dplyr::mutate(paramindex = as.numeric(factor(param)))  # Will pull out levels
    return(df)
  }

  snpdflist <- lapply(snpdflist, lagPOS)
  snpdf <- snpdflist %>% bind_rows

  plotobj <- snpdf %>%
    ggplot() +
    geom_rect(mapping=aes(xmin=start, xmax=end, ymin=paramindex-0.49, ymax=paramindex+0.49, fill=GTfact)) +
    scale_y_continuous(breaks = min(snpdf$paramindex):max(snpdf$paramindex),
                       labels = levels(factor(snpdf$param))) +
    scale_fill_manual("Genotype Call", values=c("#AA0A3C", "#313695", "#fee090", "#8c510a", "#01665e", "#cccccc")) +
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



