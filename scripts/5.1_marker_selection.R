#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#    Functions for Marker Selection tab     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~#
#      MAF      #
#~~~~~~~~~~~~~~~#


#~~ Get minor allele frequencies from SNP data

# Still need to account for different types of SNP dataframes

GetMAFs <- function(df) {
  
 # if(snp_col == T){
  #  SNP_ID <- colnames(df)
  #  df <- as.matrix(df)
  #  df <- t(df)
  #} else{
  #  df <- as.matrix(df)
  #}
  
  SNP_ID <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  
  
  ## calc_n
  n0 <- apply(df==0,1,sum,na.rm=T)
  n1 <- apply(df==1,1,sum,na.rm=T)
  n2 <- apply(df==2,1,sum,na.rm=T)
  
  n <- n0 + n1 + n2
  
  ## calculate allele frequencies
  p <- ((2*n0)+n1)/(2*n)
  q <- 1 - p
  MAF <- pmin(p, q)
  
  out <- data.frame(SNP_ID) %>%
    mutate(MAF = MAF)
  

}


#~~ Plot minor allele frequency distribution from SNP data

PlotMAF <- function(df) {
  
  df_maf <- GetMAFs(df)
  
  # Plot

  # Theme
  ggthemr(palette = "pale", layout = "clean", 
          line_weight = 0.7, text_size = 20, type = "outer")
  
  # plot maf dist
  MAFplot <- ggplot(df_maf, aes(x=MAF)) +
    geom_histogram(aes(y=..density..), alpha=0.9, col = "grey33", fill = "grey33", binwidth = 0.06) +
    geom_density(alpha=0.6, col = "grey33", fill = "grey33") +
    labs(y = "Density", x = "MAF") +
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  
  return(MAFplot)
  
}



#~~ Filter SNP data based on max and min minor allele frequencies

FilterMaf <- function(df, min_maf = 0.2, max_maf = 0.4) {
  
  #if(snp_col == T){
  #  SNP_ID <- colnames(df)
  #  df <- as.matrix(df)
  #  df <- t(df)
  #} else{
  #  df <- as.matrix(df)
  #}
  
  SNP_ID <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  
  ## calc_n
  n0 <- apply(df==0,1,sum,na.rm=T)
  n1 <- apply(df==1,1,sum,na.rm=T)
  n2 <- apply(df==2,1,sum,na.rm=T)
  
  n <- n0 + n1 + n2
  
  ## calculate allele frequencies
  p <- ((2*n0)+n1)/(2*n)
  q <- 1 - p
  MAF <- pmin(p, q)
  
  

  out <- SNP_ID %>%
    cbind(data.frame(df)) %>%
    mutate(MAF = MAF) %>%
    filter(MAF >= min_maf & MAF <= max_maf) %>%
    select(-MAF) %>%
    remove_rownames %>% column_to_rownames(var=".")
  
  out <- as.data.frame(t(out))
  out <- out


}



#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Genotyping Rate     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~ Remove SNPs with a genotyping rate below threshold

Filter_SNPgeno <- function(df, geno_rate = 0.9) {
  
  SNP_ID <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  
  out <- data.frame(df) %>%
    mutate(SNP_MISS = rowSums(is.na(.)) / ncol(.),
           SNP_GENO = 1 - SNP_MISS) %>%
    cbind(data.frame(SNP_ID)) %>%
    filter(SNP_GENO >= geno_rate) %>%
    select(-c(SNP_MISS, SNP_GENO)) %>%
    remove_rownames %>% column_to_rownames(var="SNP_ID")
  
  out <- as.data.frame(t(out))
  out <- out

  
}



#~~ Remove individuals with a genotyping rate below threshold


Filter_INDgeno <- function(df, geno_rate = 0.9) {
  
  out <- data.frame(df) %>%
    mutate(IND_MISS = (rowSums(is.na(.)) / ncol(.)),
           IND_GENO = 1 - IND_MISS) %>%
    filter(IND_GENO >= geno_rate) %>%
    select(-c(IND_MISS, IND_GENO))

  out <- out
  
  
}



#~~~~~~~~~~~~~~~~#
#      sMLH      #
#~~~~~~~~~~~~~~~~#

# mouse_msats <- read.csv("data/mouse_msats.csv", row.names = T)

# mouse_geno <- convert_raw


# smlh / mlh
# distribution
# subsampling


#~~~~~~~~~~~~~~~~#
#      HWE       #
#~~~~~~~~~~~~~~~~#

#~~ HWE Chi-Square test

GetHWE_chisq <- function(df) {

  SNP_ID <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  
  n0 <- apply(df==0,1,sum,na.rm=T)
  n1 <- apply(df==1,1,sum,na.rm=T)
  n2 <- apply(df==2,1,sum,na.rm=T)
  n <- n0+n1+n2
  
  obs <- cbind(n0, n1, n2)
  p <- ((2*n0)+n1)/(2*n)
  q <- (1-p)
  
  expected <- cbind(p*p, 2*p*q, q*q)
  expected <- expected*n
  chisq <- (obs-expected)
  chisq <- (chisq*chisq) /expected
  chisq <- apply(chisq,1,sum)
  chisq.p <- 1-pchisq(chisq,df=1)
  

  out <- data.frame(SNP_ID) %>%
    mutate(Estimate = chisq,
           P_value = chisq.p)

  out <- out
  

}



#~~ Filter SNP data based on chisq p_values

FilterHWE_chisq <- function(df, pvalue = 0.0001) {
  
  SNP_ID <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  
  n0 <- apply(df==0,1,sum,na.rm=T)
  n1 <- apply(df==1,1,sum,na.rm=T)
  n2 <- apply(df==2,1,sum,na.rm=T)
  n <- n0+n1+n2
  
  obs <- cbind(n0, n1, n2)
  p <- ((2*n0)+n1)/(2*n)
  q <- (1-p)
  
  expected <- cbind(p*p, 2*p*q, q*q)
  expected <- expected*n
  chisq <- (obs-expected)
  chisq <- (chisq*chisq) /expected
  chisq <- apply(chisq,1,sum)
  chisq.p <- 1-pchisq(chisq,df=1)
  


  out <- SNP_ID %>%
    cbind(data.frame(df)) %>%
    mutate(P_value = chisq.p) %>%
    filter(P_value > pvalue) %>%
    select(-P_value) %>%
    remove_rownames %>% column_to_rownames(var=".")
  
    out <- out
    
}


#~~ HWE Fisher's Exact Test


GetHWE_fisher <- function(df) {
  
  SNP_ID <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  
  n0 <- apply(df==0,1,sum,na.rm=T)
  n1 <- apply(df==1,1,sum,na.rm=T)
  n2 <- apply(df==2,1,sum,na.rm=T)

  z <- cbind(n0, ceiling(n1/2), floor(n1/2), n2)
  z <- lapply( split( z, 1:nrow(z) ), matrix, ncol=2 )
  z <- lapply( z, fisher.test )
  
  
  out <- cbind(OddsRatio = as.numeric(unlist(lapply(z, "[[", "estimate"))),
               P_value = as.numeric(unlist(lapply(z, "[[", "p.value"))))
  
  out <- cbind(data.frame(SNP_ID), out)

}




#~~ Filter SNP data based on Fisher's exact test p_values

FilterHWE_fisher <- function(df, pvalue = 0.0001) {
  
  SNP_ID <- colnames(df)
  df <- as.matrix(df)
  df <- t(df)
  
  n0 <- apply(df==0,1,sum,na.rm=T)
  n1 <- apply(df==1,1,sum,na.rm=T)
  n2 <- apply(df==2,1,sum,na.rm=T)
  
  z <- cbind(n0, ceiling(n1/2), floor(n1/2), n2)
  z <- lapply( split( z, 1:nrow(z) ), matrix, ncol=2 )
  z <- lapply( z, fisher.test )
  
  
  out <- cbind(OddsRatio = as.numeric(unlist(lapply(z, "[[", "estimate"))),
               P_value = as.numeric(unlist(lapply(z, "[[", "p.value"))))
  
  out <- cbind(data.frame(SNP_ID), out) %>%
    cbind(data.frame(df)) %>%
    filter(P_value > pvalue) %>%
    select(-c(OddsRatio, P_value)) %>%
    remove_rownames %>% column_to_rownames(var="SNP_ID")
  
  out <- out
  
}



#~~~~~~~~~~~~~~~~~~~#
#      Ho / He      #
#~~~~~~~~~~~~~~~~~~~#




