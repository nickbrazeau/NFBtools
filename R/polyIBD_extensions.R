

#------------------------------------------------
#' @title tidy polyIBD sim as snpdf
#' @details this is temporary
#'
#' @export

polyIBDsim2snpdf <- function(polyIBDsim){

  haploid1 <- data.frame(polyIBDsim$haploid$haploid1)
  colnames(haploid1) <- c(paste0("Smpl1_Haplotype", seq(1, ncol(haploid1))))
  haploid2 <- data.frame(polyIBDsim$haploid$haploid2)
  colnames(haploid2) <- c(paste0("Smpl2_Haplotype", seq(1, ncol(haploid1))))
  IBD <- data.frame(polyIBDsim$IBD)
  IBD[IBD == 0] <- "U"
  IBD[IBD == 1] <- "I"
  colnames(IBD) <- c(paste0("HaploPair_IBD", seq(1, ncol(haploid1))))

  snpdf <- cbind.data.frame(CHROM = polyIBDsim$CHROMPOS$CHROM,
                            POS = polyIBDsim$CHROMPOS$POS,
                            polyIBDsim$gtmatrix,
                            haploid1,
                            haploid2,
                            IBD)
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
}



#------------------------------------------------
#' @title Plot polyIBD sim as gtcov mat
#' @details this is temporary
#'
#' @export


polyIBDsimsnpdf2GTcov <- function(snpdf){
plot <- snpdf %>%
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

  return(plot)

}



