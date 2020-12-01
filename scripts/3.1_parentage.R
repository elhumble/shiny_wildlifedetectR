# library(dplyr)
# library(inbreedR)
# library(data.table)
# library(plyr)
# library(tidyr)
# library(readxl)

#~~~~~~~~~~~~~~~~~~~~~~#
#     SNP Functions    #
#~~~~~~~~~~~~~~~~~~~~~~#


alleleFreqs_plink <- function(db){
  
  # Prepare database of allele freqs from plink raw
  
  # Remove first six plink cols
  db <- db %>%
    dplyr::select(-c(FID, IID, PAT, MAT, SEX, PHENOTYPE))
  
  # Missing data to NA
  db[(db!=0) & (db!=1) & (db!=2)] <- NA
  
  SNP_ID <- colnames(db)
  df <- as.matrix(db)
  df <- t(df)
  
  ## calc_n
  n0 <- apply(df==0,1,sum,na.rm=T)
  n1 <- apply(df==1,1,sum,na.rm=T)
  n2 <- apply(df==2,1,sum,na.rm=T)
  
  n <- n0 + n1 + n2
  
  ## calculate allele frequencies of db
  p <- ((2*n0)+n1)/(2*n)
  q <- 1 - p
  freqs <- data.frame(SNP_ID = SNP_ID, p = p, q = q) %>%
    mutate(SNP_ID = as.character(SNP_ID))

  return(freqs)
}



prep_snp_trio <- function(trio) {
  
    # Remove first six plink cols
  trio <- trio %>%
    dplyr::select(-c(FID, IID, PAT, MAT, SEX, PHENOTYPE))
  
  SNP_ID <- colnames(trio)
  
  # Missing data to NA
  trio[(trio!=0) & (trio!=1) & (trio!=2)] <- NA
  
  # Individuals by column and add column of locus names
  trio <- as.matrix(trio)
  trio <- as.data.frame(t(trio)) %>%
    mutate(SNP_ID = SNP_ID)
  
  # Remove loci with missing data (strict)
  trio <- trio %>%
    filter(!is.na(V3) &
             !is.na(V1) &
             !is.na(V2))
  
  return(trio)
}



prep_snp_duo <- function(duo) {
  
  # Remove first six plink cols
  duo <- duo %>%
    dplyr::select(-c(FID, IID, PAT, MAT, SEX, PHENOTYPE))
  
  SNP_ID <- colnames(duo)
  
  # Missing data to NA
  duo[(duo!=0) & (duo!=1) & (duo!=2)] <- NA
  
  # Individual by column and add column of locus names
  duo <- as.matrix(duo)
  duo <- as.data.frame(t(duo)) %>%
    mutate(SNP_ID = SNP_ID)
  
  # Remove loci with missing data (strict)
  duo <- duo %>%
    filter(!is.na(V1) &
             !is.na(V2))
  
  return(duo)
}


#~~~~~~~~~~~~~~~~~~~~~~#
#     Test 1 SNPs      #
#~~~~~~~~~~~~~~~~~~~~~~#

test1_snps <- function(df = df, 
                       theta = 0.1){
  
  # Likelihood analysis
  Th <- theta
  
  df <- df %>%
    mutate(prob = case_when(V1 == 0 & V2 == 0 & V3 == 0 ~ as.character(((4*Th)+((1-Th)*(p)))*((5*Th)+((1-Th)*(p)))/((1+(3*Th))*(1+(4*Th)))), # Eq1
                            V1 == 2 & V2 == 2 & V3 == 2 ~ as.character(((4*Th)+((1-Th)*(q)))*((5*Th)+((1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq1
                            V1 == 0 & V2 == 1 & V3 == 1 ~ as.character((4*((2*Th)+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th)))), # Eq2
                            V1 == 2 & V2 == 1 & V3 == 1 ~ as.character((4*((2*Th)+(1-Th)*(q))*((3*Th)+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq2
                            V1 == 1 & V2 == 0 & V3 == 1 ~ as.character(4*((3*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq3
                            V1 == 1 & V2 == 1 & V3 == 0 ~ as.character(4*((3*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq3
                            V1 == 1 & V2 == 2 & V3 == 1 ~ as.character(4*((3*Th)+((1-Th)*(q)))*(Th+((1-Th)*(p)))/((1+(3*Th))*(1+(4*Th)))), # Eq3 
                            V1 == 1 & V2 == 1 & V3 == 2 ~ as.character(4*((3*Th)+((1-Th)*(q)))*(Th+((1-Th)*(p)))/((1+(3*Th))*(1+(4*Th)))), # Eq3 
                            V1 == 1 & V2 == 1 & V3 == 1 ~ as.character(4*((2*Th)+((1-Th)*(p)))*((2*Th)+((1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq5
                            V1 == 0 & V2 == 0 & V3 == 1 ~ as.character(2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*(p)))/((1+(3*Th))*(1+(4*Th)))), # Eq6
                            V1 == 0 & V2 == 1 & V3 == 0 ~ as.character(2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*(p)))/((1+(3*Th))*(1+(4*Th)))), # Eq6
                            V1 == 2 & V2 == 2 & V3 == 1 ~ as.character(2*((3*Th)+((1-Th)*(q)))*((4*Th)+((1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq6
                            V1 == 2 & V2 == 1 & V3 == 2 ~ as.character(2*((3*Th)+((1-Th)*(q)))*((4*Th)+((1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq6
                            V1 == 1 & V2 == 0 & V3 == 2 ~ as.character(2*((2*Th)+((1-Th)*(p)))*((2*Th)+((1-Th)*(q)))/((1+(3*Th))*(1+(4*Th)))), # Eq7 
                            V1 == 1 & V2 == 2 & V3 == 0 ~ as.character(2*((2*Th)+((1-Th)*(q)))*((2*Th)+((1-Th)*(p)))/((1+(3*Th))*(1+(4*Th)))), # Eq7 
                            
                            # If we want to add in a mutation rate/error into these equations, I think we multiply the denominator by the rate (e.g. 0.001)
                            # I think you normally multiply the numerator, but these are upside down right? (because of the way Rob did it in his spreadsheet)
                            
                            # Exclusion scenarios. Assigning an 'E' to any loci that 'cause' an exlusion. 
                            
                            V1 == 0 & V2 == 2 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 0 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 1 & V2 == 0 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 1 & V2 == 2 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 2 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 0 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 0 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 2 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 0 & V3 == 1 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 1 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 1 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 2 & V3 == 1 ~ as.character('E'))) # Exclusion
  
  ### To exclude, or to calculate likelihood ratio ###
  
  # If there is a 'E' probability in the dataframe, print out an exclusion statement. If not, calculate the likelihood ratio.
  if('E' %in% df$prob) {
    # Print how many loci are 'causing' the exclusion, and list them.
    # Position of loci
    loci_pos <- which('E' == df$prob)
    # Name of loci
    loci_name <- SNP_ID[loci_pos]
    loci_name <- as.data.frame(loci_name)
    
    # Print statement, then number of loci 'causing' the exclusion. Note, I think there is better teminology we can use than 'causing' an exclusion!
    cat("Exclusion. Number of loci causing the exclusion - ", nrow(loci_name),"\n\nLoci that caused the exclusion:\n")
    # Then print list of loci.
    print(loci_name)
    
    #If we want to include a mutation rate/error threshold (e.g. 1%, or 0.01). 
   # cat("\nChance of seeing at least one mutation: \n", ((0.01)^(nrow(df)))*100, "%")
    # We can integrate a mutation/error rate much better; this is just an example. And I'm not 100% sure that it's even correct.
    
    # Incorporating a mutaion rate:
    ### Defaults:
    # microsat mutaiton rate = 10^-3
    # SNP mutaiton rate = 10^-8
    # E.g.
    mutation_rate <- (10^(-8))
    cat("\nChance of seeing ", nrow(loci_name), "mismatches based on the mutaiton rate: \n",
       # (1-((1-((mutation_rate)^(nrow(loci_name))))^(nrow(df))))*100, "%")
        (1-((1-((mutation_rate)^(nrow(df))))^(nrow(loci_name))))*100, "%")

    
    
    ## Calculate likelihood ratio
    
  } else {
    
    
    # Multiplying all loci probabilities together to get the cumulative probability
    cuml_prob <- prod(as.numeric(df$prob), na.rm = T)
    # Now taking the reciprocal for the cumulative likelihood
    cuml_likelihood <- 1/cuml_prob
    cat("Cumulative likelihood: ", cuml_likelihood)
    #cat("\nRMP: ", cuml_prob)
  }
  
}


#~~~~~~~~~~~~~~~~~~~~~~#
#     Test 2 SNPs      #
#~~~~~~~~~~~~~~~~~~~~~~#
  
test2_snps <- function(df = df, 
                       theta = 0.1) {
  
  Th <- theta
  
  df <- df %>%
    mutate(prob = case_when(V1 == 0 & V2 == 0 & V3 == 0 ~ as.character((1+(3*Th))/((4*Th)+((1-Th)*(p)))), #Eq9
                            V1 == 2 & V2 == 2 & V3 == 2 ~ as.character((1+(3*Th))/((4*Th)+((1-Th)*(q)))), #Eq9
                            V1 == 0 & V2 == 1 & V3 == 1 ~ as.character((1+(3*Th))/(2*((2*Th)+((1-Th)*(p))))), #Eq4
                            V1 == 2 & V2 == 1 & V3 == 1 ~ as.character((1+(3*Th))/(2*((2*Th)+((1-Th)*(q))))), #Eq4
                            V1 == 1 & V2 == 0 & V3 == 1 ~ as.character((1+(3*Th))/(2*(Th+((1-Th)*(q))))), #Eq2
                            V1 == 1 & V2 == 1 & V3 == 0 ~ as.character((1+(3*Th))/((4*Th)+((1-Th)*(p+q)))), #Eq14
                            V1 == 1 & V2 == 2 & V3 == 1 ~ as.character((1+(3*Th))/(2*(Th+((1-Th)*(p))))), #Eq2
                            V1 == 1 & V2 == 1 & V3 == 2 ~ as.character((1+(3*Th))/((4*Th)+((1-Th)*(p+q)))), #Eq14
                            V1 == 1 & V2 == 1 & V3 == 1 ~ as.character((1+(3*Th))/((4*Th)+((1-Th)*(p+q)))), #Eq15
                            V1 == 0 & V2 == 0 & V3 == 1 ~ as.character((1+(3*Th))/(2*((3*Th)+((1-Th)*(p))))), #Eq1
                            V1 == 0 & V2 == 1 & V3 == 0 ~ as.character((1+(3*Th))/((3*Th)+((1-Th)*(p)))), #Eq10
                            V1 == 2 & V2 == 2 & V3 == 1 ~ as.character((1+(3*Th))/(2*((3*Th)+((1-Th)*(q))))), #Eq1
                            V1 == 2 & V2 == 1 & V3 == 2 ~ as.character((1+(3*Th))/((3*Th)+((1-Th)*(q)))), #Eq10
                            V1 == 1 & V2 == 0 & V3 == 2 ~ as.character((1+(3*Th))/((2*Th)+((1-Th)*(q)))), #Eq11
                            V1 == 1 & V2 == 2 & V3 == 0 ~ as.character((1+(3*Th))/((2*Th)+((1-Th)*(p)))), #Eq11
                            
                            # If we want to add in a mutation rate/error into these equations, I think we multiply the numerator by the rate (e.g. 0.001)
                            
                            # Exclusion scenarios. Assigning an 'E' to any loci that 'cause' an exlusion. 
                            
                            V1 == 0 & V2 == 2 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 0 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 1 & V2 == 0 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 1 & V2 == 2 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 2 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 0 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 0 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 2 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 0 & V3 == 1 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 1 & V3 == 0 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 1 & V3 == 2 ~ as.character('E'), # Exclusion
                            V1 == 0 & V2 == 2 & V3 == 1 ~ as.character('E'))) # Exclusion
  
  
  ### To exclude, or to calculate likelihood ratio ###
  
  # If there is a 'E' probability in the dataframe, print out an exclusion statement. If not, calculate the likelihood ratio.
  if('E' %in% df$prob) {
    # Print how many loci are 'causing' the exclusion, and list them.
    # Position of loci
    loci_pos <- which('E' == df$prob)
    # Name of loci
    loci_name <- SNP_ID[loci_pos]
    loci_name <- as.data.frame(loci_name)
    
    # Print statement, then number of loci 'causing' the exclusion. Note, I think there is better teminology we can use than 'causing' an exclusion!
    cat("Exclusion. Number of loci causing the exclusion - ", nrow(loci_name),"\n\nLoci that caused the exclusion:\n")
    # Then print list of loci.
    print(loci_name)
    
    # Incorporating a mutaion rate:
    ### Defaults:
    # microsat mutaiton rate = 10^-3
    # SNP mutaiton rate = 10^-8
    # E.g.
    mutation_rate <- (10^(-8))
    cat("\nChance of seeing ", nrow(loci_name), "mismatches based on the mutaiton rate: \n",
        (1-((1-((mutation_rate)^(nrow(df))))^(nrow(loci_name))))*100, "%")
    
    ## Calculate likelihood ratio
    
  } else {
    
    # Calculation likelihood ratio
    # Multiplying all loci probabilities together to get the cumulative likelihood
    cuml_likelihood <- prod(as.numeric(df$prob), na.rm = T)
    cat("Cumulative likelihood: ", cuml_likelihood)

  }
  
}

  
#~~~~~~~~~~~~~~~~~~~~~~#
#     Test 3 SNPs      #
#~~~~~~~~~~~~~~~~~~~~~~#
  
test3_snps <- function(df = df,
                       theta = 0.1){
  
  Th <- theta
  df <- left_join(duo, freqs, by = "SNP_ID") %>%
    mutate(prob = case_when(V1 == 0 & V2 == 0 ~ as.character((1+(2*Th))/(((1-Th)*(p))+(3*Th))), #Eq5
                            V1 == 2 & V2 == 2 ~ as.character((1+(2*Th))/(((1-Th)*(q))+(3*Th))), #Eq5
                            V1 == 1 & V2 == 0 ~ as.character((1+(2*Th))/(2*(((1-Th)*(p))+(2*Th)))), #Eq3
                            V1 == 1 & V2 == 2 ~ as.character((1+(2*Th))/(2*(((1-Th)*(q))+(2*Th)))), #Eq3
                            V1 == 0 & V2 == 1 ~ as.character((1+(2*Th))/(2*(((1-Th)*(p))+(2*Th)))), #Eq4
                            V1 == 2 & V2 == 1 ~ as.character((1+(2*Th))/(2*(((1-Th)*(q))+(2*Th)))), #Eq4
                            V1 == 1 & V2 == 1 ~ as.character(((1+(2*Th))*((((p)+(q))*(1-Th))+(2*Th)))/(4*(((1-Th)*(p))+Th)*(((1-Th)*(q))+Th))), #Eq2
                            
                            # If we want to add in a mutation rate/error into these equations, I think we multiply the numerator by the rate (e.g. 0.001)
                            
                            # Exclusion scenarios. Assigning an 'E' to any loci that 'cause' an exlusion.
                            
                            V1 == 0 & V2 == 2 ~ as.character('E'), # Exclusion
                            V1 == 2 & V2 == 0 ~ as.character('E'))) # Exclusion
  
  ### To exclude, or to calcualte likelihood ratio ###
  
  # If there is a 'E' probability in the dataframe, print out an exclusion statement. If not, calculate the likelihood ratio.
  if('E' %in% df$prob) {
    # Print how many loci are 'causing' the exclusion, and list them.
    # Position of loci
    loci_pos <- which('E' == df$prob)
    # Name of loci
    loci_name <- SNP_ID[loci_pos]
    loci_name <- as.data.frame(loci_name)
    
    # Print statement, then number of loci 'causing' the exclusion. Note, I think there is better teminology we can use than 'causing' an exclusion!
    cat("Exclusion. Number of loci causing the exclusion - ", nrow(loci_name),"\n\nLoci that caused the exclusion:\n")
    # Then print list of loci.
    print(loci_name)
    
    # Incorporating a mutaion rate:
    ### Defaults:
    # microsat mutaiton rate = 10^-3
    # SNP mutaiton rate = 10^-8
    # E.g.
    mutation_rate <- (10^(-8))
    cat("\nChance of seeing ", nrow(loci_name), "mismatches based on the mutaiton rate: \n",
        (1-((1-((mutation_rate)^(nrow(df))))^(nrow(loci_name))))*100, "%")
    
    ## Calculation likelihood ratio
    
  } else {
    # Multiplying all loci probabilities together to get the cumulative likelihood
    cuml_likelihood <- prod(as.numeric(df$prob), na.rm = T)
    cat("Cumulative likelihood: ", cuml_likelihood)
  }
  
} 



#~~~~~~~~~~~~~~~~~~~~~~#
#      MS Functions    #
#~~~~~~~~~~~~~~~~~~~~~~#

#path <- "data/dog_db_case_1.xlsx"

# db <- path %>%
#  excel_sheets() %>%
#  set_names() %>%
#  purrr::map(read_excel, path = path)

#unknown <- read_excel("data/dog_trio_case_1.xlsx")

#~~ Frequency data


freqs2freqs <- function(db) {
  
  db_size <- db[[1]] %>%
    gather(Locus, Size, -Allele) %>%
    na.omit() %>%
    .[c(2, 1, 3)] %>%
    group_by(Locus)
  
  # gather allele freqs into one column
  
  db_freq <- db[[2]] %>%
    gather(Locus, Freq, -Allele) %>%
    na.omit() %>%
    .[c(2,1,3)] %>%
    group_by(Locus) 
  
  # join dataframes
  freqs <- left_join(db_freq, db_size) %>% ### numeric here...check this...
    mutate(Size = as.character(Size))
  
  return(freqs)
  
}


freqs2freqs_genalex <- function(db) {
  
  db <- db[[2]] %>%
    .[-c(1:10,12),-1]
  loci <- c("Size", db[1,c(2:ncol(db))])
  
  freqs <- db %>%
    .[-1,] %>%
    `colnames<-`(unlist(loci)) %>%
    gather(Locus, Freq, -Size) %>%
    drop_na() %>%
    mutate(Freq = as.numeric(Freq))
  
  return(freqs)

}


genos2freqs <- function(db){
  
  freqs <- db %>%
    select(-ID) %>%
    gather(Locus, Size) %>%
    separate(Locus, c("Locus", "Dip")) %>%
    group_by(Locus, Size) %>%
    dplyr::summarise(n = n()) %>%
    group_by(Locus) %>%
    mutate(Freq = n/sum(n)) %>%
    select(-n) %>%
    mutate(Size = as.character(Size))
  
  return(freqs)
  
}



#~~ Top 2 alleles

top2freq <- function(db) {
  
  db <- db %>%
    group_by(Locus) %>%
    top_n(2, Freq) %>%
    arrange(Locus, dplyr::desc(Freq)) %>%
    ungroup() %>%
    dplyr::mutate(Allele = rep(c("p", "q"), nrow(.)/2)) %>%
    select(Locus, Allele, Freq) %>%
    spread(Allele, Freq)

  return(db)
}


#~~ Allele frequencies / genos matching


extract_known_alleles <- function(unknown, freqs) {
  
  un <- unknown %>%
    .[1,] %>%
    gather(Locus, Size, -ID) %>%
    separate(Locus, c("Locus", "Dip")) %>%
    arrange(Locus, Size) %>% # arrange loci by size
    mutate(Dip = rep(c("p", "q"), nrow(.)/2)) %>%
    mutate(Size = as.character(Size)) %>%
    dplyr::left_join(freqs, by = c("Locus", "Size")) %>%
    select(Locus, Dip, Freq) %>%
    spread(Dip, Freq) %>%
    select(Locus, p, q)
  
  return(un)
  
}



#~~ Prepare trio genotypes

prep_ms_trio <- function(trio){
  
  #trio <- read_excel("data/badger_trio_case_1.xlsx")
  
  inds <- trio[1]
  trio <- select(trio, -ID)
  
  n_inds <- nrow(trio)
  n_loc <- (ncol(trio))/2
  loc_names <- unique(gsub("..$|..$|_.$|_.$|-.$|-.$", "", colnames(trio))) # def need a more general regex here
  
  n <- ncol(trio)
  df1 <- trio[seq(1,n,2)] %>%
    gather()
  df2 <- trio[seq(2,n,2)] %>%
    gather()
  
  # arrange alleles by size
  df <- cbind(df1, df2)
  
  df <- df %>%
    `colnames<-`(c("one", "size1", "two", "size2")) %>%
    mutate(p = case_when(size1 >= size2 ~ size2,
                         size1 < size2 ~ size1),
           q = case_when(size1 <= size2 ~ size2,
                         size1 > size2 ~ size1)) %>%
    unite(allele, c("p", "q"), sep = "_") %>%
    mutate(ind = rep(1:n_inds, n_loc)) %>%
    select(-c(size1, size2, two)) %>%
    separate(one, c("locus", "sup")) %>%
    select(-sup) %>%    
    mutate(locus = factor(locus,levels=unique(loc_names))) %>% 
    spread(ind, allele) %>%
    select(-locus)
  
  df <- t(df)
  colnames(df) <- loc_names
  
  df <- cbind(inds, df)
  return(df)
  
}


prep_ms_duo <- function(duo){
  
  duo <- duo %>%
    .[c(1:2),]
  
  inds <- duo[1]
  duo <- select(duo, -ID)
  
  n_inds <- nrow(duo)
  n_loc <- (ncol(duo))/2
  loc_names <- unique(gsub("..$|..$|_.$|_.$|-.$|-.$", "", colnames(duo))) # def need a more general regex here
  
  n <- ncol(duo)
  df1 <- duo[seq(1,n,2)] %>%
    gather()
  df2 <- duo[seq(2,n,2)] %>%
    gather()
  
  # arrange alleles by size
  df <- cbind(df1, df2)
  
  df <- df %>%
    `colnames<-`(c("one", "size1", "two", "size2")) %>%
    mutate(p = case_when(size1 >= size2 ~ size2,
                         size1 < size2 ~ size1),
           q = case_when(size1 <= size2 ~ size2,
                         size1 > size2 ~ size1)) %>%
    unite(allele, c("p", "q"), sep = "_") %>%
    mutate(ind = rep(1:n_inds, n_loc)) %>%
    select(-c(size1, size2, two)) %>%
    separate(one, c("locus", "sup")) %>%
    select(-sup) %>%    
    mutate(locus = factor(locus,levels=unique(loc_names))) %>% 
    spread(ind, allele) %>%
    select(-locus)
  
  df <- t(df)
  colnames(df) <- loc_names
  
  df <- cbind(inds, df)
  df
  
}


#~~ Recode ms alleles

numericSNPmat <- function(g) {
  
  # if(ploidy(g) != 2) stop("'g' must have diploid data")
  
  # recode pp, pq, qq, qp, qr etc. to 01234...
  
  nums <- c('1.1' = 0, '1.2' = 1, '2.1' = 1, '2.2' = 2, '3.3' = 3, '1.3' = 4, '3.1' = 4,
            '2.3' = 5, '3.2' = 5, '4.4' = 6, '1.4' = 7, '4.1' = 7, '2.4' = 8, '4.2' = 8,
            '3.4' = 9, '4.3' = 9)
  
  #g.mat <- as.matrix(g, sep = "_", strata = FALSE, one.col = TRUE)
  
  mat <- apply(g[,-1], 2, function(loc) { # apply over columns ()
    loc <- strsplit(loc, split = "_")
    lvls <- unique(unlist(loc))
    lvls <- lvls[!is.na(lvls)]
    #if(length(lvls) > 2) return(rep(NA, nrow(g.mat)))
    unname(sapply(loc, function(x) nums[paste(match(x, lvls), collapse = ".")]))
    
  })
  
  rownames(mat) <- mat[, 1]
  
  to.delete <- which(apply(mat, 2, function(x) all(is.na(x))))
  
  if(length(to.delete) > 0) {
    mat <- mat[, -to.delete]
    to.delete <- paste(colnames(mat)[to.delete], collapse = ", ")
    warning("The following loci have been removed: ", to.delete)
    
  }
  mat
}


#~~~~~~~~~~~~~~~~~~~~~~#
#      Test 1 MS       #
#~~~~~~~~~~~~~~~~~~~~~~#

# Likelihood analysis
# Need to triple check I've covered all possible scenarios

test1_ms <- function(df, Th = 0.1) {
  
  Th <- Th
  
  # incorporate an NA error. E.g. if p and q == NA stop.
  
  df <- df %>%
    mutate(prob = case_when(V1 == 0 & V2 == 0 & V3 == 0 ~((4*Th)+((1-Th)*(p)))*((5*Th)+((1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            
                            # V1 == 2 & V2 == 2 & V3 == 2 ~ ((1+(3*Th))*(1+(4*Th))/(4*Th)+((1-Th)*q))*((5*Th)+((1-Th)*q)), # This will never be the case as the child's alleles are always assigned first
                            
                            V1 == 0 & V2 == 0 & V3 == 1 ~ (2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 0 & V3 == 4 ~ (2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 0 & V3 == 7 ~ (2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 1 & V3 == 0 ~ (2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 4 & V3 == 0 ~ (2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 7 & V3 == 0 ~ (2*((3*Th)+((1-Th)*(p)))*((4*Th)+((1-Th)*p)))/((1+(3*Th))*(1+(4*Th))),
                            
                            V1 == 0 & V2 == 1 & V3 == 1 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 1 & V3 == 4 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 1 & V3 == 7 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 4 & V3 == 1 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 4 & V3 == 4 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 4 & V3 == 7 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 7 & V3 == 1 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 7 & V3 == 4 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 0 & V2 == 7 & V3 == 7 ~ (4*(2*Th+(1-Th)*(p))*((3*Th)+(1-Th)*(p)))/((1+(3*Th))*(1+(4*Th))),
                            
                            V1 == 1 & V2 == 0 & V3 == 1 ~ (4*((3*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 1 & V3 == 0 ~ (4*((3*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 1 & V3 == 2 ~ (4*((3*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 2 & V3 == 1 ~ (4*((3*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            
                            V1 == 1 & V2 == 0 & V3 == 2 ~ (2*((2*Th)+((1-Th)*(p)))*((2*Th)+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 2 & V3 == 0 ~ (2*((2*Th)+((1-Th)*(p)))*((2*Th)+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            
                            V1 == 1 & V2 == 0 & V3 == 5 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 0 & V3 == 8 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 4 & V3 == 2 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 7 & V3 == 2 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 2 & V3 == 4 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 2 & V3 == 7 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 5 & V3 == 0 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 8 & V3 == 0 ~ (4*((2*Th)+(1-Th)*(p))*(Th+(1-Th)*(q)))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            
                            V1 == 1 & V2 == 1 & V3 == 1 ~ (4*((2*Th)+((1-Th)*(p)))*((2*Th)+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            
                            V1 == 1 & V2 == 1 & V3 == 4 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 1 & V3 == 7 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 1 & V3 == 5 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 1 & V3 == 8 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 4 & V3 == 1 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 7 & V3 == 1 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 5 & V3 == 1 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            V1 == 1 & V2 == 8 & V3 == 1 ~ (8*((2*Th)+((1-Th)*(p)))*(Th+((1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))), # I think this is correct
                            
                            V1 == 1 & V2 == 4 & V3 == 5 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 4 & V3 == 8 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 7 & V3 == 5 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 7 & V3 == 8 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 5 & V3 == 4 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 5 & V3 == 7 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 8 & V3 == 4 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th))),
                            V1 == 1 & V2 == 8 & V3 == 7 ~ (8*((Th+(1-Th)*(p)))*((Th+(1-Th)*(q))))/((1+(3*Th))*(1+(4*Th)))),
           # Exclusion scenarios
           exclude = case_when(V1 == 0 & V2 == 0 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 0 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 1 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 1 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 1 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 4 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 4 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 4 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 7 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 7 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 7 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 0 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 1 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 4 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 7 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 0 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 1 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 4 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 7 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 0 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 1 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 4 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 7 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 8 ~ 'Y',
                               V1 == 1 & V2 == 0 & V3 == 0 ~ 'Y',
                               V1 == 1 & V2 == 0 & V3 == 4 ~ 'Y',
                               V1 == 1 & V2 == 0 & V3 == 7 ~ 'Y',
                               V1 == 1 & V2 == 4 & V3 == 0 ~ 'Y',
                               V1 == 1 & V2 == 4 & V3 == 4 ~ 'Y',
                               V1 == 1 & V2 == 4 & V3 == 7 ~ 'Y',
                               V1 == 1 & V2 == 7 & V3 == 0 ~ 'Y',
                               V1 == 1 & V2 == 7 & V3 == 4 ~ 'Y',
                               V1 == 1 & V2 == 7 & V3 == 7 ~ 'Y',
                               V1 == 1 & V2 == 2 & V3 == 2 ~ 'Y',
                               V1 == 1 & V2 == 2 & V3 == 5 ~ 'Y',
                               V1 == 1 & V2 == 2 & V3 == 8 ~ 'Y',
                               V1 == 1 & V2 == 5 & V3 == 2 ~ 'Y',
                               V1 == 1 & V2 == 5 & V3 == 5 ~ 'Y',
                               V1 == 1 & V2 == 5 & V3 == 8 ~ 'Y',
                               V1 == 1 & V2 == 8 & V3 == 2 ~ 'Y',
                               V1 == 1 & V2 == 8 & V3 == 5 ~ 'Y',
                               V1 == 1 & V2 == 8 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 1 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 2 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 5 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 1 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 1 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 1 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 1 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 2 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 5 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 1 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 2 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 5 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 2 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 2 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 2 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 5 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 5 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 5 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 8 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 8 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 8 & V3 == 7 ~ 'Y',
                               TRUE ~ "N"))
                            
  
  ### To exclude, or to calculate likelihood ratio ###
  
  # If there is a 'E' probability in the dataframe, print out an exclusion statement. If not, calculate the likelihood ratio.
  # if('Y' %in% df$exclude) {
  #   # Print how many loci are 'causing' the exclusion, and list them.
  #   # Position of loci
  #   loci_pos <- which('Y' == df$exclude)
  #   # Name of loci
  #   loci_name <- SNP_ID[loci_pos]
  #   loci_name <- as.data.frame(loci_name)
  #   
  #   # Print statement, then number of loci 'causing' the exclusion. Note, I think there is better teminology we can use than 'causing' an exclusion!
  #   cat("Exclusion. Number of loci causing the exclusion - ", nrow(loci_name),"\n\nLoci that caused the exclusion:\n")
  #   # Then print list of loci.
  #   print(loci_name)
  #   
  #   # Incorporating a mutaion rate:
  #   ### Defaults:
  #   # microsat mutaiton rate = 10^-3
  #   # SNP mutation rate = 10^-8
  #   # E.g.
  #   mutation_rate <- (10^(-3))
  #   cat("\nChance of seeing ", nrow(loci_name), "mismatches based on the mutation rate: \n",
  #       (1-((1-((mutation_rate)^(nrow(loci_name))))^(nrow(df))))*100, "%")
  # 
  #   ### Calculate likelihood ratio ###
  #   
  # } else {
  #   # Now multiplying all loci probabilities together to get the cumulative probability
  #   cuml_prob <- prod(df$prob, na.rm = T)
  #   
  #   # Now taking the reciprocal for the cumulative likelihood
  #   cuml_likelihood <- 1/cuml_prob
  #   #cuml_likelihood
  #   return(cuml_prob)
  #   
  # }
  
  return(df)
  
}


#~~~~~~~~~~~~~~~~~~~~~~#
#      Test 2 MS       #
#~~~~~~~~~~~~~~~~~~~~~~#


test2_ms <- function(df, Th = 0.1) {
  
  Th <- Th
  
  df <- df %>%
    mutate(prob = case_when(V1 == 0 & V2 == 0 & V3 == 0 ~ (1+(3*Th))/((4*Th)+(p*(1-Th))),
                            
                            V1 == 0 & V2 == 0 & V3 == 1 ~ (1+(3*Th))/(2*((3*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 0 & V3 == 4 ~ (1+(3*Th))/(2*((3*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 0 & V3 == 7 ~ (1+(3*Th))/(2*((3*Th)+(p*(1-Th)))),
                            
                            V1 == 0 & V2 == 1 & V3 == 0 ~ (1+(3*Th))/((3*Th)+(p*(1-Th))),
                            V1 == 0 & V2 == 4 & V3 == 0 ~ (1+(3*Th))/((3*Th)+(p*(1-Th))),
                            V1 == 0 & V2 == 7 & V3 == 0 ~ (1+(3*Th))/((3*Th)+(p*(1-Th))),
                            
                            V1 == 0 & V2 == 1 & V3 == 1 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 1 & V3 == 4 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 1 & V3 == 7 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 4 & V3 == 1 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 4 & V3 == 4 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 4 & V3 == 7 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 7 & V3 == 1 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 7 & V3 == 4 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            V1 == 0 & V2 == 7 & V3 == 7 ~ (1+(3*Th))/(2*((2*Th)+(p*(1-Th)))),
                            
                            V1 == 1 & V2 == 0 & V3 == 1 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 0 & V3 == 5 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 0 & V3 == 8 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 2 & V3 == 1 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            V1 == 1 & V2 == 2 & V3 == 4 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            V1 == 1 & V2 == 2 & V3 == 7 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            
                            V1 == 1 & V2 == 0 & V3 == 2 ~ (1+(3*Th))/((2*Th)+(q*(1-Th))),
                            V1 == 1 & V2 == 4 & V3 == 2 ~ (1+(3*Th))/((2*Th)+(q*(1-Th))),
                            V1 == 1 & V2 == 7 & V3 == 2 ~ (1+(3*Th))/((2*Th)+(q*(1-Th))),
                            V1 == 1 & V2 == 2 & V3 == 0 ~ (1+(3*Th))/((2*Th)+(p*(1-Th))),
                            V1 == 1 & V2 == 5 & V3 == 0 ~ (1+(3*Th))/((2*Th)+(p*(1-Th))),
                            V1 == 1 & V2 == 8 & V3 == 0 ~ (1+(3*Th))/((2*Th)+(p*(1-Th))),
                            
                            V1 == 1 & V2 == 1 & V3 == 0 ~ (1+(3*Th))/((4*Th)+((1-Th)*(p+q))),
                            V1 == 1 & V2 == 1 & V3 == 1 ~ (1+(3*Th))/((4*Th)+((1-Th)*(p+q))),
                            V1 == 1 & V2 == 1 & V3 == 2 ~ (1+(3*Th))/((4*Th)+((1-Th)*(p+q))),
                            
                            V1 == 1 & V2 == 1 & V3 == 4 ~ (1+(3*Th))/(2*((3*Th)+((1-Th)*(p+q)))),
                            V1 == 1 & V2 == 1 & V3 == 7 ~ (1+(3*Th))/(2*((3*Th)+((1-Th)*(p+q)))),
                            V1 == 1 & V2 == 1 & V3 == 5 ~ (1+(3*Th))/(2*((3*Th)+((1-Th)*(p+q)))),
                            V1 == 1 & V2 == 1 & V3 == 8 ~ (1+(3*Th))/(2*((3*Th)+((1-Th)*(p+q)))),
                            
                            V1 == 1 & V2 == 4 & V3 == 1 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 4 & V3 == 5 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 4 & V3 == 8 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 7 & V3 == 1 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 7 & V3 == 5 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 7 & V3 == 8 ~ (1+(3*Th))/(2*(Th+(q*(1-Th)))),
                            V1 == 1 & V2 == 5 & V3 == 1 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            V1 == 1 & V2 == 5 & V3 == 4 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            V1 == 1 & V2 == 5 & V3 == 7 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            V1 == 1 & V2 == 8 & V3 == 1 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            V1 == 1 & V2 == 8 & V3 == 4 ~ (1+(3*Th))/(2*(Th+(p*(1-Th)))),
                            V1 == 1 & V2 == 8 & V3 == 7 ~ (1+(3*Th))/(2*(Th+(p*(1-Th))))),
           
           # Exclusion scenarios. Assigning an 'E' to any loci that 'cause' an exlusion.
           
           exclude = case_when(V1 == 0 & V2 == 0 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 0 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 0 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 1 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 1 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 1 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 4 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 4 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 4 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 7 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 7 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 7 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 0 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 1 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 4 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 7 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 2 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 0 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 1 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 4 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 7 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 5 & V3 == 8 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 0 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 1 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 4 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 7 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 2 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 5 ~ 'Y',
                               V1 == 0 & V2 == 8 & V3 == 8 ~ 'Y',
                               V1 == 1 & V2 == 0 & V3 == 0 ~ 'Y',
                               V1 == 1 & V2 == 0 & V3 == 4 ~ 'Y',
                               V1 == 1 & V2 == 0 & V3 == 7 ~ 'Y',
                               V1 == 1 & V2 == 4 & V3 == 0 ~ 'Y',
                               V1 == 1 & V2 == 4 & V3 == 4 ~ 'Y',
                               V1 == 1 & V2 == 4 & V3 == 7 ~ 'Y',
                               V1 == 1 & V2 == 7 & V3 == 0 ~ 'Y',
                               V1 == 1 & V2 == 7 & V3 == 4 ~ 'Y',
                               V1 == 1 & V2 == 7 & V3 == 7 ~ 'Y',
                               V1 == 1 & V2 == 2 & V3 == 2 ~ 'Y',
                               V1 == 1 & V2 == 2 & V3 == 5 ~ 'Y',
                               V1 == 1 & V2 == 2 & V3 == 8 ~ 'Y',
                               V1 == 1 & V2 == 5 & V3 == 2 ~ 'Y',
                               V1 == 1 & V2 == 5 & V3 == 5 ~ 'Y',
                               V1 == 1 & V2 == 5 & V3 == 8 ~ 'Y',
                               V1 == 1 & V2 == 8 & V3 == 2 ~ 'Y',
                               V1 == 1 & V2 == 8 & V3 == 5 ~ 'Y',
                               V1 == 1 & V2 == 8 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 1 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 2 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 5 ~ 'Y',
                               V1 == 2 & V2 == 0 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 1 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 1 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 1 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 1 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 2 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 5 ~ 'Y',
                               V1 == 2 & V2 == 4 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 1 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 2 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 5 ~ 'Y',
                               V1 == 2 & V2 == 7 & V3 == 8 ~ 'Y',
                               V1 == 2 & V2 == 2 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 2 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 2 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 5 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 5 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 5 & V3 == 7 ~ 'Y',
                               V1 == 2 & V2 == 8 & V3 == 0 ~ 'Y',
                               V1 == 2 & V2 == 8 & V3 == 4 ~ 'Y',
                               V1 == 2 & V2 == 8 & V3 == 7 ~ 'Y',
                               TRUE ~ "N"))
  
  ### To exclude, or to calculate likelihood ratio ###
  
  # If there is a 'E' probability in the dataframe, print out an exclusion statement. If not, calculate the likelihood ratio.
  # if('Y' %in% df$exclude) {
  #   # Print how many loci are 'causing' the exclusion, and list them.
  #   # Position of loci
  #   loci_pos <- which('Y' == df$exclude)
  #   # Name of loci
  #   loci_name <- SNP_ID[loci_pos]
  #   loci_name <- as.data.frame(loci_name)
  #   
  #   # Print statement, then number of loci 'causing' the exclusion. Note, I think there is better teminology we can use than 'causing' an exclusion!
  #   cat("Exclusion. Number of loci causing the exclusion - ", nrow(loci_name),"\n\nLoci that caused the exclusion:\n")
  #   # Then print list of loci.
  #   print(loci_name)
  #   
  #   # Incorporating a mutaion rate:
  #   ### Defaults:
  #   # microsat mutaiton rate = 10^-3
  #   # SNP mutaiton rate = 10^-8
  #   # E.g.
  #   mutation_rate <- (10^(-3))
  #   cat("\nChance of seeing ", nrow(loci_name), "mismatches based on the mutaiton rate: \n",
  #       (1-((1-((mutation_rate)^(nrow(loci_name))))^(nrow(df))))*100, "%")
  #   
  #   
  #   ### Calculate likelihood ratio ###
  #   
  # } else {
  #   
  #   # Now multiplying all loci probabilities together to get the cumulative probability
  #   df <- df %>%
  #     mutate(round = round(prob, 2))
  #   
  #   cuml_prob <- prod(df$round, na.rm = T)
  #   
  #   # Now taking the reciprocal for the cumulative likelihood
  #   cuml_likelihood <- 1/cuml_prob
  #   return(cuml_prob)
  #   #cuml_likelihood
  # }
  return(df)
}



#~~~~~~~~~~~~~~~~~~~~~~#
#      Test 3 MS       #
#~~~~~~~~~~~~~~~~~~~~~~#

test3_ms <- function(df, Th = 0.1){
  
  Th <- Th
  
  df <- df %>%
    mutate(prob = case_when(V1 == 0 & V2 == 0 ~ (1+(2*Th))/(((1-Th)*(p)+(3*Th))),
                            
                            V1 == 0 & V2 == 1 ~ (1+(2*Th))/(2*(((1-Th)*(p)+(2*Th)))),
                            V1 == 0 & V2 == 4 ~ (1+(2*Th))/(2*(((1-Th)*(p)+(2*Th)))),
                            V1 == 0 & V2 == 7 ~ (1+(2*Th))/(2*(((1-Th)*(p)+(2*Th)))),
                            V1 == 1 & V2 == 0 ~ (1+(2*Th))/(2*(((1-Th)*(p)+(2*Th)))),
                            V1 == 1 & V2 == 2 ~ (1+(2*Th))/(2*(((1-Th)*(q)+(2*Th)))),
                            
                            V1 == 1 & V2 == 1 ~ ((1+(2*Th))*((((p)+(q))*(1-Th))+(2*Th)))/(4*(((1-Th)*(p))+Th)*(((1-Th)*(q))+Th)),
                            
                            V1 == 1 & V2 == 4 ~ (1+(2*Th))/(4*((1-Th)*(p)+Th)),
                            V1 == 1 & V2 == 7 ~ (1+(2*Th))/(4*((1-Th)*(p)+Th)),
                            V1 == 1 & V2 == 5 ~ (1+(2*Th))/(4*((1-Th)*(q)+Th)),
                            V1 == 1 & V2 == 8 ~ (1+(2*Th))/(4*((1-Th)*(q)+Th))),
           
           # Exclusion scenarios. Assigning an 'E' to any loci that 'cause' an exlusion.
           
           exclude = case_when(V1 == 0 & V2 == 2 ~ 'Y',
                               V1 == 0 & V2 == 3 ~ 'Y',
                               V1 == 0 & V2 == 5 ~ 'Y',
                               V1 == 0 & V2 == 6 ~ 'Y',
                               V1 == 0 & V2 == 8 ~ 'Y',
                               V1 == 0 & V2 == 9 ~ 'Y',
                               TRUE ~ "N"))
  
  ### To exclude, or to calculate likelihood ratio ###
  
  # If there is a 'E' probability in the dataframe, print out an exclusion statement. If not, calculate the likelihood ratio.
  # if('Y' %in% df$exclude) {
  #   # Print how many loci are 'causing' the exclusion, and list them.
  #   # Position of loci
  #   loci_pos <- which('Y' == df$exclude)
  #   # Name of loci
  #   loci_name <- SNP_ID[loci_pos]
  #   loci_name <- as.data.frame(loci_name)
  #   
  #   # Print statement, then number of loci 'causing' the exclusion. Note, I think there is better teminology we can use than 'causing' an exclusion!
  #   cat("Exclusion. Number of loci causing the exclusion - ", nrow(loci_name),"\n\nLoci that caused the exclusion:\n")
  #   # Then print list of loci.
  #   print(loci_name)
  #   
  #   # Incorporating a mutaion rate:
  #   ### Defaults:
  #   # microsat mutaiton rate = 10^-3
  #   # SNP mutaiton rate = 10^-8
  #   # E.g.
  #   mutation_rate <- (10^(-3))
  #   cat("\nChance of seeing ", nrow(loci_name), "mismatches based on the mutaiton rate: \n",
  #       (1-((1-((mutation_rate)^(nrow(loci_name))))^(nrow(df))))*100, "%")
  #   
  #   
  #   ### Calculate likelihood ratio ###
  #   
  # } else {
  #   
  #   df <- df %>%
  #     mutate(prob_round = round(prob,2))
  #   
  #   # Now multiplying all loci probabilities together to get the cumulative probability
  #   cuml_prob <- prod(df$prob, na.rm = T)
  #   
  #   # Now taking the reciprocal for the cumulative likelihood
  #   cuml_likelihood <- 1/cuml_prob
  #   return(cuml_prob)
  #   cuml_likelihood
  #   
  #   }
  return(df)
}
  








