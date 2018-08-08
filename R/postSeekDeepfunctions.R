#' @title check if seekdeep dat
#' @noRd

is.SeekDeepDat <- function(x){inherits(x, "SeekDeepDat")}


#' @title skdp_filter_simplifier
#'
#' @description filters
#'
#' @param skdpclustinfo_df from seekdeep
#'
#' @export

SeekDeepOutput2SeekDeepDat <- function(skdpclustinfo_df, readcountcutoff = 0){
  # Drops clusters without a lot of read support and recalculates haplotype/cluster fractions.
  #
  # Args:
  #   skdpclustinfo_df: From SeekDeep Process Clusters -- selectedClustersInfo.tab.txt
  #   readcountcutoff: Read Depth Cutoff
  #
  # Returns:
  #  simplier and filtered df



  # Filter
  skdpclustinfo_df_simp <- skdpclustinfo_df
  skdpclustinfo_df_simp <- skdpclustinfo_df_simp[skdpclustinfo_df_simp$c_ReadCnt > readcountcutoff, ] # filter at the cluster level (which is the within sample haplotype)

  # Error handling
  if(nrow(skdpclustinfo_df_simp) > 0){
    filtered_c_ReadCnt_sum <- aggregate(skdpclustinfo_df_simp$c_ReadCnt ~ skdpclustinfo_df_simp$s_Sample, function(x){sum(x)}, data = skdpclustinfo_df_simp) # get denominator
    colnames(filtered_c_ReadCnt_sum) <- c("s_Sample", "filtered_c_ReadCnt_denom")
    skdpclustinfo_df_simp <-  dplyr::left_join(skdpclustinfo_df_simp, filtered_c_ReadCnt_sum, by=c("s_Sample"))



    skdpclustinfo_df_simp$c_AveragedFrac_adj <- skdpclustinfo_df_simp$c_ReadCnt/skdpclustinfo_df_simp$filtered_c_ReadCnt_denom # adjusted average fraction by cluster


    skdpclustinfo_df_simp <- skdpclustinfo_df_simp[, colnames(skdpclustinfo_df_simp) %in% c("s_Sample", "c_AveragedFrac_adj", "h_popUID", "c_Consensus", "filtered_c_ReadCnt_denom", "c_ReadCnt")] # keep specific columns

  } else{
    stop("There was an error filtering the reads. Contact the developer")
  }

  class(skdpclustinfo_df_simp) <- append(  class(skdpclustinfo_df_simp) , "SeekDeepDat" )
  return(skdpclustinfo_df_simp)


}


#---------------------------------------------------------------------------------

#' @title SeekDeepDat2HapPlotter
#'
#' @description plots
#'
#' @param input from skdp_filter_simplifier
#'
#' @export

SeekDeepDat2HapPlotter <- function(input, target="Target"){

  # error handle
  if(!is.SeekDeepDat(input)){
    stop("Input must be of class SeekDeepDat See the SeekDeepOutput2SeekDeepDat function.")
  }


  skdpclustinfo_df_simp <- input # NFB fix this to input throughout...

  # Color setup
  # stackoverflow, # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  n <- length(unique(skdpclustinfo_df_simp$h_popUID))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] # pul out qualitative paletes
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col=sample(col_vector, n)


  # Drops clusters without a lot of read support and recalculates haplotype/cluster fractions.
  #
  # Args:
  #   skdpclustinfo_df_simp: Plots output from above
  #
  # Returns:
  #  stacked ggplot



  skdpclustinfo_df_simp %>%
    ggplot() +
    geom_bar(aes(y = c_AveragedFrac_adj, x = s_Sample, fill = h_popUID), stat="identity") +
    ggtitle(paste("Haplotype (Cluster) Frequencies by s_Sample", target)) +
    xlab("Sample") +
    ylab("Haplotype/Cluster Frequency") +
    guides(fill=guide_legend(title="Haplotyper (Cluster)")) +
    theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, family = "Helvetica", face = "bold", size=14)) +
    theme(axis.text.y=element_text(family = "Helvetica", face = "bold", size=8)) +
    theme(axis.title.x=element_text(family = "Helvetica", face = "bold", size=10)) +
    theme(axis.title.y=element_text(family = "Helvetica", face = "bold", size=10)) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size=18)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = col)

}


#---------------------------------------------------------------------------------

#' @title SeekDeepDat2ExonAnnotation
#'
#' @description annotates exons
#'
#' @param skdpclustinfo_df_simp from skdp_filter_simplifier
#'
#' @export

SeekDeepDat2ExonAnnotation <- function(input,
                                     gff, geneid,
                                     ampliconrefseqpath, forwardprimerpath, reverseprimerpath,
                                     ncbigeneticcode=1){
  #
  # IMPORTANT -- this function assumes that after you have trimmed your primers you are only in an exonic region
  # i.e. we are only reading genes from gff
  #
  # Args:
  #   skdpclustinfo_df: From SeekDeep Process Clusters -- selectedClustersInfo.tab.txt
  #   gffpath: path to reference gff
  #   ampliconrefseqpath: seekdeep will provide this after extraction under targetRefSeqs
  #   forwardprimerpath: seekdeep will provide this after extraction under targetRefSeqs
  #   reverseprimerpath: seekdeep will provide this after extraction under targetRefSeqs
  #
  # Returns:
  #  annotated df



  # error handle
  if(!is.SeekDeepDat(input)){
    stop("Input must be of class SeekDeepDat See the SeekDeepOutput2SeekDeepDat function.")
  }


  skdpclustinfo_df_simp <- input # NFB fix this to input throughout...


  #####  READ IN GFF GENE INFORMATION
  geneid.end <- grep('##FASTA',readLines(gff))
  geneid.gff <- read.delim(file = gff,
                           nrows = geneid.end-1, comment= "#", header=F)
  colnames(geneid.gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "info")


  geneid.gff <- subset(geneid.gff, geneid.gff$feature == "gene") # subset to only genes, don't want the other mRNA, etc data

  ############################################
  ##### EXTRACT GENEID from GFF INFO   #######
  ############################################
  geneid.gff$GeneID <- stringr::str_split_fixed(geneid.gff$info, ";", n=2)[,1] # give it to columns to parse on
  geneid.gff$GeneID <- gsub("ID=", "", geneid.gff$GeneID, fixed=T)
  # These are the gene identifiers that we care about

  ############################################
  #####    Subset Gene from GFF        #######
  ############################################
  geneid.gff <- data.frame(geneid.gff[geneid.gff$GeneID == geneid, ], stringsAsFactors = F)

  # read in gff fasta
  gffseq <- seqinr::read.fasta(file=gff, seqtype = c("DNA"), forceDNAtolower = F, strip.desc = T) # of note this is a hacky solution...review it and improve

  # subset
  gffseq <- gffseq[names(gffseq) %in% c(geneid)]



  ####################################
  #####     REFERNET FASTA     #######
  ####################################
  ampliconrefseq <- Biostrings::readDNAStringSet(filepath = ampliconrefseqpath, format="fasta")
  ampliconrefseq <- ampliconrefseq[grepl("Pf3d7", names(ampliconrefseq))]


  ############################################
  #### Now Trim Forward and Reverse Primer  ##
  ############################################
  Lprimer <- Biostrings::readDNAStringSet(filepath=forwardprimerpath, format = "fasta")
  Rprimer <- Biostrings::reverseComplement(Biostrings::readDNAStringSet(filepath=reverseprimerpath, format = "fasta"))

  Pf3D7haplotypeRef <- Biostrings::trimLRPatterns(Lpattern = Lprimer[[1]],
                                                  Rpattern = Rprimer[[1]],
                                                  subject = ampliconrefseq)  # trim off primers

  #### THIS IS NOW THE REFERENT HAPLOTYPE (primers trimmed that we can compare/align to for variants)

  ############################################
  ####    Make a Dataframe of Variants     ###
  ############################################

  skdpconsens_dnalist <- split(skdpclustinfo_df_simp, factor(skdpclustinfo_df_simp$h_popUID))

  compare_DNA <- function(skdpclustinfo_df_simp,dnastringobject){
    x <- as.integer(Biostrings::DNAString(skdpclustinfo_df_simp$c_Consensus[1]))
    y <- as.integer(Biostrings::DNAString(dnastringobject))

    SNPpos <- which(x != y) #https://www.biostars.org/p/16880/
    SNP <- seqinr::s2c(skdpclustinfo_df_simp$c_Consensus[1])[which(x != y)] # nucleotide bp that are associated with snps

    SNPdf <- data.frame(h_popUID=rep(skdpclustinfo_df_simp$h_popUID[1], length(SNPpos)),
                        AmpPos=SNPpos, SNP=SNP)
    return(SNPdf)   # a single consensus cluster may have several variants -- need to account for this --> will do left_join
  }

  # find snps
  skdpconsens_SNPlist <- lapply(X=skdpconsens_dnalist, FUN=compare_DNA, Pf3D7haplotypeRef[[1]])

  # call snps
  skdpconsens_SNP <- do.call("rbind", skdpconsens_SNPlist)
  ############################################
  ##### error handling if no snps in amplicon
  ############################################
  if(nrow(skdpconsens_SNP) !=0){


    ############################################
    ###    Identify First base in Amplicon   ###
    ###    from "global" Gene Fasta          ###
    ############################################

    RefSeqGenePos <- BioStrings::matchPattern(pattern = Pf3D7haplotypeRef[[1]],
                                  subject = DNAString(seqinr::c2s(gffseq[[1]]))) # find pos of amplicon in gene


    ### SANITY CHECK
    ampliconfull <- ampliconrefseq@ranges@width
    ampliconprimertrim <- ampliconfull - Rprimer@ranges@width - Lprimer@ranges@width
    if(ampliconprimertrim != Pf3D7haplotypeRef@ranges@width){
      stop("The amplicon and haplotype are of different lengths. There was an issue in primer triming. See primer trim in code")
    }


    ############################################
    ### Updated Amplicon Position to Gene Pos ##
    ############################################
    skdpconsens_SNP$GenePos <- skdpconsens_SNP$AmpPos + RefSeqGenePos@ranges@start - 1 # minus one here so the first base doesn't get counted twice in our amp count and our vcf count

    ############################################
    ####    Referent Amplicon w/in Gene   ######
    ############################################
    AArefseq <- seqinr::translate(gffseq[[1]], sens=F, numcode=ncbigeneticcode)

    ############################################
    ####    Make Mutant Amplicon Haplotype  ####
    ############################################
    haplist <- split(skdpconsens_SNP, f=factor(skdpconsens_SNP$h_popUID))

    mutatehap <- function(haplistobj){
      mutseq <- gffseq[[1]]
      mutseq[haplistobj$GenePos] <- as.character(haplistobj$SNP) # make mutation seq by putting in all SNPs
      AAmutseq <- seqinr::translate(mutseq, sens=F, numcode=ncbigeneticcode) # translate

      haplistobj$CODON <- ceiling(haplistobj$GenePos/3) # ceiling will so that 1/3, 2/3, 3/3 all got to #1 CODON
      haplistobj$AAREF <- AArefseq[haplistobj$CODON]
      haplistobj$AAALT <- AAmutseq[haplistobj$CODON]
      haplistobj$MUT_Type <- ifelse(haplistobj$AAREF == haplistobj$AAALT, "Syn", "Nonsyn")
      return(haplistobj)
    }

    skdpconsens_SNP <- do.call("rbind", lapply(haplist, mutatehap))

    #################################################################
    #######                        RETURN                       #####
    #################################################################
    skdpvcf <- dplyr::left_join(skdpclustinfo_df_simp, skdpconsens_SNP, by=c("h_popUID"))

  } else { skdpvcf <- NULL }

  return(skdpvcf)
}
