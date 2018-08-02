#' @title vcfR2vcffilter_CHROMPOS
#'
#'
#' @export

#-----------------------------------------------------
# Filter VCF based on given positions among the chromosomes
#------------------------------------------------------
vcffilter_CHROMPOS <- function(vcfRobject = NULL,
                               chromposdf = NULL # this is expected to be from GFF2VariantAnnotation_Short
){
  chromposdflist = split(chromposdf, f=factor(chromposdf$GeneID))
  chromposlong <- do.call("rbind", parallel::mclapply(chromposdflist, function(dat){

    i <- length(seq(from=dat$start, to=dat$end))

    temp <- tibble::tibble(CHROM = rep(dat$seqname, i),
                           POS = seq(from=dat$start, to=dat$end),
                           info = rep(dat$info, i),
                           GeneID = rep(dat$GeneID, i),
                           Description = rep(dat$Description),
                           remove = rep("Y", i)
    )
    return(temp)

  }))


  vcftidy <- vcfR2tidy(vcfRobject)

  passloci <- vcftidy$fix %>%
    dplyr::select(CHROM, POS) %>%
    dplyr::left_join(x=., y=chromposlong) %>%
    dplyr::mutate(keep = ifelse(is.na(remove), TRUE, FALSE)) %>%
    dplyr::select(keep) %>%
    as.matrix(.)

  meta <- append(vcfRobject@meta, paste("##Additional Filters for excluding hypervariable regions determined by user"))
  fix <- vcfRobject@fix[ passloci, ]
  gt <- vcfRobject@gt[ passloci, ]


  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  newvcfR

}


#' @title vcfR2vcffilter_INFO
#'
#'
#' @export

## note if there is total agreement in the alignments, i.e. 1/1 where the reads are 0,100 0,100 etc then the
## readposranksum and mqranksum won't be calculated since there is no p-value (there is no read diff)
## this can happen when the reference genome is different from samples and there are snps where
## samples are in total agreement but diff from ref

#-----------------------------------------------------
# Filter VCF based on INFO Fields
#------------------------------------------------------

vcffilter_info <- function(vcfRobject = NULL,
                           infoMQ=55,
                           infoQD=25,
                           infoSOR=2,
                           infoAF = 0.05,
                           infoDP = NULL, # this is a percentile cutoff
                           infoFS = NULL,
                           infoMQRankSum = NULL,
                           infoReadPosRankSum = NULL,
                           biallelic = TRUE,
                           SNPs = FALSE){

  require(vcfR)
  require(tidyverse)

  vcf <- vcfRobject # legacy
  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(biallelic==TRUE){
    vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  }
  if(SNPs==TRUE){
    vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  }

  # store loci objects on info fields
  infolist <- list()
  if(!is.null(infoMQ)){
    infolist <- append(infolist, "MQ")
  }
  if(!is.null(infoQD)){
    infolist <- append(infolist, "QD")
  }
  if(!is.null(infoSOR)){
    infolist <- append(infolist, "SOR")
  }
  if(!is.null(infoAF)){
    infolist <- append(infolist, "AF")
  }
  if(!is.null(infoDP)){
    infolist <- append(infolist, "DP")
  }
  if(!is.null(infoFS)){
    infolist <- append(infolist, "FS")
  }

  if(!is.null(infoMQRankSum)){
    infolist <- append(infolist, "MQRankSum")
  }
  if(!is.null(infoReadPosRankSum)){
    infolist <- append(infolist, "ReadPosRankSum")
  }

# extract format list information
tidyvcf <- vcfR::vcfR2tidy(vcf)
infodf <- tidyvcf$fix


#--------------------------------------------------------
# filter loci
#--------------------------------------------------------
if(!is.null(infoMQ)){
  infodf <- infodf %>%
    dplyr::mutate(MQ = ifelse(MQ < infoMQ, "DROP", MQ))
}
if(!is.null(infoQD)){
  infodf <- infodf %>%
    dplyr::mutate(QD = ifelse(QD < infoQD, "DROP", QD))
}
if(!is.null(infoSOR)){
  infodf <- infodf %>%
    dplyr::mutate(SOR = ifelse(SOR > infoSOR, "DROP", SOR))
}
if(!is.null(infoAF)){
  infodf <- infodf %>%
    dplyr::mutate(AF = ifelse(AF< infoAF, "DROP", AF))
}
if(!is.null(infoDP)){
  DPpercentile <- quantile(infodf$DP, c(infoDP, 1-infoDP))
  infodf <- infodf %>%
    dplyr::mutate(DP = ifelse(DP >= DPpercentile[1] & DP <= DPpercentile[2],
                              DP, "DROP"))
}
if(!is.null(infoFS)){
  infodf <- infodf %>%
    dplyr::mutate(FS = ifelse(FS > infoFS, "DROP", FS))
}
if(!is.null(infoMQRankSum)){
  infodf <- infodf %>%
    dplyr::mutate(MQRankSum = ifelse(MQRankSum < infoMQRankSum, "DROP", MQRankSum))
}
if(!is.null(infoReadPosRankSum)){
  infodf <- infodf %>%
    dplyr::mutate(ReadPosRankSum = ifelse(ReadPosRankSum < infoReadPosRankSum, "DROP", ReadPosRankSum))
}

#--------------------------------------------------------
# apply filter
#--------------------------------------------------------

infodf <- infodf[ , colnames(infodf) %in% c("CHROM", "POS", infolist) ]
passedloci <- infodf %>%
  dplyr::mutate(CHROMPOS = paste0(CHROM, POS)) %>%
  dplyr::mutate(excl = apply(., 1, function(x){all(x != "DROP")})) %>%
  dplyr::select(excl)

passedloci <- passedloci$excl

fix <- as.matrix(vcfR::getFIX(vcf, getINFO = T)[ passedloci ,])
gt <- as.matrix(vcf@gt[ passedloci , ])
meta <- append(vcf@meta, "##Additional Filters provided by NFB filter tools")

# Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

newvcfR

}







#' @title vcfR2vcffilter_FORMAT
#'
#'
#' @export

#-----------------------------------------------------
# Filter VCF based on Format Fields
#------------------------------------------------------
vcffilter_format <- function(vcfRobject = NULL,
                             formatGQ=NULL,
                             formatDP = NULL,
                             formatSP = NULL,
                             prop.loci.missing = 0.5, # this is given a loci, how many samples can have missing information before it is dropped
                             biallelic = TRUE,
                             SNPs = FALSE){

  require(vcfR)
  require(tidyverse)

  vcf <- vcfRobject # legacy
  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(biallelic==TRUE){
    vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  }
  if(SNPs==TRUE){
    vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  }
  # store loci objects on info fields
  formatlist <- list()
  if(!is.null(formatGQ)){
    formatlist <- append(formatlist, "gt_GQ")
  }
  if(!is.null(formatDP)){
    formatlist <- append(formatlist, "gt_DP")
  }
  if(!is.null(formatSP)){
    formatlist <- append(formatlist, "gt_SP")
  }

  # extract format list information
  tidyvcf <- vcfR::vcfR2tidy(vcf)
  formatdf <- tidyvcf$gt


  #--------------------------------------------------------
  # filter loci
  #--------------------------------------------------------
  if(!is.null(formatGQ)){
    formatdf <- formatdf %>%
      dplyr::mutate(gt_GQ = ifelse(gt_GQ < formatGQ, "DROP", gt_GQ))
  }
  if(!is.null(formatDP)){
    formatdf <- formatdf %>%
      dplyr::mutate(gt_DP = ifelse(gt_DP < formatDP, "DROP", gt_DP))
  }
  if(!is.null(formatSP)){
    formatdf <- formatdf %>%
      dplyr::mutate(gt_SP = ifelse(gt_SP > formatSP, "DROP", gt_SP)) # this is a fisher-score p-value for likelihood of strand bias. Higher worse
  }

  #--------------------------------------------------------
  # apply filter
  #--------------------------------------------------------

  formatdf <- formatdf[ , colnames(formatdf) %in% c("ChromKey", "POS", "Indiv", formatlist) ]
  formatdf <- formatdf %>%
    dplyr::mutate(excl = apply(., 1, function(x){any(x == "DROP")})) %>%
    dplyr::select(ChromKey, POS, Indiv, excl) %>%
    tidyr::spread(., key="Indiv", value="excl") %>%
    dplyr::select(-c(ChromKey, POS)) %>%
    cbind(FORMAT=rep(FALSE, nrow(.)), .)

  vcf@gt[as.matrix(formatdf)] <- NA

  #--------------------------------------------------------
  # Subset by loci missingness
  #--------------------------------------------------------
  if(!is.null(prop.loci.missing)){
    locimissingness <- apply(vcf@gt, 1, function(x){sum(is.na(x))})
    locimissingnessprop <- locimissingness/(ncol(vcf@gt)-1)


    # filter loci with too much missing information
    vcf@gt <- vcf@gt[ locimissingnessprop < prop.loci.missing , ]
  }

  fix <- as.matrix(vcfR::getFIX(vcf, getINFO = T)[ locimissingnessprop < prop.loci.missing ,])
  gt <- as.matrix(vcf@gt)
  meta <- append(vcf@meta, "##Additional Filters provided by NFB filter tools")

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  newvcfR

}


#' @title vcfR2vcffilter_SEGSITES
#'
#'
#' @export
#-----------------------------------------------------
# Read vcfR object to segregated sites
#------------------------------------------------------

vcfR2segsites <- function(vcfRobj = NULL, err = 0.025){
  vcf <- vcfRobj # legacy
  if(!identical(vcf, vcf[is.biallelic(vcf)])){
    stop("VCF must be biallelic for this to work.")
  }

  ad <- vcfR::extract.gt(vcf, element = "AD")
  refad <- masplit(ad, record=1, sort=0, decreasing = 0)
  altad <- masplit(ad, record=2, sort=0, decreasing = 0)

  NRAF <- altad/(altad + refad)
  NRAF[is.nan(NRAF)] <- NA # the 0,0 are returning Nans

  segsites <- apply(NRAF, 1, function(x){
      if(all(is.na(x))){
        return(FALSE)
      } else {
        return( (max(x, na.rm=T) - min(x, na.rm=T)) > err )

      }

  })

  vcfRobj@gt <- vcfRobj@gt[segsites,]

  fix <- as.matrix(vcfR::getFIX(vcfRobj, getINFO = T)[segsites,])
  gt <- as.matrix(vcfRobj@gt)
  meta <- append(vcfRobj@meta, paste("##Additional Filters for segregating sites, such that within-sample AF must vary more than:", err))

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  newvcfR

}


##--------------------------------------------------------------------
## LINKAGE DISEQUILIBRIUM
##--------------------------------------------------------------------

#' @title cormat
#' @description Calculates genetic autocorrelation among variants.
#' @detials Expects a matrix of allele frequencies with rows as loci and columns as samples
#'
#'@author Bob Verity & Nick Brazeau
#'
#' @export
#'

#--------------------------------------------------------
# Calculate pairwise correlation matrix between observations. Input is a matrix with n rows corresponding to n multivariate measurements, output is a n-by-n correlation matrix. NA values are imputed as the mean.
#--------------------------------------------------------
corMat <- function(m) {
  tol <- 1e-9 # tolerance for denominator if AF become the same
  n <- nrow(m) # number of loci
  c <- matrix(NA,n,n)
  for (i in 1:n) {
    x1 <- unlist(m[i,])
    mu1 <- mean(x1,na.rm=TRUE)
    for (j in i:n) {
      x2 <- unlist(m[j,])
      mu2 <- mean(x2,na.rm=TRUE)
      c[i,j] <- c[j,i] <- sum((x1-mu1)*(x2-mu2),na.rm=TRUE)/sqrt( sum((x1-mu1)^2,na.rm=TRUE)*sum((x2-mu2)^2,na.rm=TRUE) + tol)
    }
  }
  return(c)
}



#' @title genautocorr setup
#' @description Calculates genetic autocorrelation for later linkage disequilibrium filtering.
#' @detials From an object of class \code{vcfR}, calculate the genetic autocorrelation as the estimate of linkage disequilibrium by genetic distance.
#' @param vcffile A variant call file (vcf) path. This VCF will be converted to an object of class \code{vcfR}.
#'
#'@author Bob Verity & Nick Brazeau
#'
#' @export
#'

genautocorr <- function(vcffile = NULL, vcfR = NULL, biallelicsnps=TRUE){

  require(vcfR)
  require(tidyverse)

  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      stop("vcfR object must be of class vcfR")
    }
    vcf <- vcfR
  } else{
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=T) # read vcf
  }
  if(biallelicsnps == FALSE){
    stop("Must take in biallelic SNPs")
  }
  vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic

  # extract the within sample allele frequencies
  ad <- vcfR::extract.gt(vcf, element = "AD")
  refad <- masplit(ad, record=1, sort=0, decreasing = 0)
  altad <- masplit(ad, record=2, sort=0, decreasing = 0)

  NRAF <- altad/(altad + refad)
  NRAF[is.nan(NRAF)] <- NA # the 0,0 are returning Nans


  # extract distances via positions from vcfR object
  CHROM <- vcfR::getCHROM(vcf)
  POS <- vcfR::getPOS(vcf)
  vcfpos <- cbind.data.frame(CHROM, POS)
  vcfdf <- cbind.data.frame(vcfpos, NRAF)

  vcflist <- split(vcfdf, vcfdf$CHROM)

  # -----------------------------------------------------
  # calculations
  #------------------------------------------------------

  cormatgendistwrapper <- function(vcfdf_fromlist){

    # get correlation matrix. NA values are imputed as the mean
    df1 <- vcfdf_fromlist[, !colnames(vcfdf_fromlist) %in% c("CHROM", "POS")]
    c <- corMat(df1)
    # get distance between SNPs. This can be extracted from the row names of the vcf
    df2 <- vcfdf_fromlist[, colnames(vcfdf_fromlist) %in% c("POS")]
    gendist <- as.matrix(dist(df2))

    ret <- list(vcfAF = vcfdf_fromlist, corMat=c, gendist=gendist)
    return(ret)
  }


  retlist <- parallel::mclapply(vcflist, cormatgendistwrapper)

  return(retlist)

}



#' @title vcfR2LDfiltered
#' @description Filtering an object of class \code{vcfR} for linkage disequilibrium via genetic autocorrelation.
#' @param vcffile A variant call file (vcf) path. This VCF will be converted to an object of class \code{vcfR}.
#'
#'
#' @export
#'

vcfR2LDfiltered <- function(vcffile = NULL, vcfR = NULL, genautocorrresult=NULL, threshDist=1e3, biallelicsnps=TRUE){

  require(vcfR)
  require(tidyverse)

  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(genautocorrresult)){
    stop("Must specify a linkage disequilibrium threshold (see help and tutorial).")
  }
  if(is.null(genautocorrresult)){
    stop("Must specify a genetic autocorrelation results object using the genautocorr function")
  }

  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      stop("vcfR object must be of class vcfR")
    }
    vcf <- vcfR
  } else{
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=T) # read vcf
  }
  if(biallelicsnps == FALSE){
    stop("This must be used with biallelic SNPs")
  }

  vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic

  vcfdf <- cbind.data.frame(vcf@fix, vcf@gt)
  vcflist <- split(vcfdf, f=factor(vcfdf$CHROM))
  if(length(vcflist) != length(genautocorrresult)){
    stop("The number of chromosomes in the vcfR objec and the results from the genetic autocorrelation analysis differ.")
  }
  # -----------------------------------------------------
  # Filter based on distance
  #------------------------------------------------------

  filter_autocorr <- function(vcflist, genautocorrresult){
    gendist <- as.matrix(genautocorrresult$gendist)
    diag(gendist) <- Inf	# block self-comparison
    while (any(gendist<threshDist)) {
      w <- which(gendist<threshDist, arr.ind=TRUE)
      vcflist <- vcflist[-w[1,1],]
      gendist <- gendist[-w[1,1],-w[1,1]]

    }
    return(vcflist)
  }

  updatedvcflist <- parallel::mcmapply(filter_autocorr, vcflist, genautocorrresult, SIMPLIFY = F)
  updatedvcfdf <- do.call(rbind, updatedvcflist)

  # -----------------------------------------------------
  # Return to vcfR object
  #------------------------------------------------------
  fix <- as.matrix(updatedvcfdf[,1:8])
  gt <- as.matrix(updatedvcfdf[,9:ncol(updatedvcfdf)])
  meta <- append(vcf@meta, paste("##Filtered for genetic autocorrelation by polyIBD filter tools with a threshold distance of", threshDist))

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)

}






