# Microsatellite individualisation

# options(digits=10)
# library(dplyr)
# library(data.table)
# library(readxl)
# library(purrr)


#~~~~~~~~~~~~~~~~~~~#
#     Function      #
#~~~~~~~~~~~~~~~~~~~#

#path <- "data/badger_case_1.xlsx"
#path <- "data/badger_case_1_genalex.xlsx"

# db <- path %>%
#   excel_sheets() %>%
#   set_names() %>%
#   purrr::map(read_excel, path = path)

#unknown <- read_excel("data/badger_unknown_case_1.xlsx")

individualisation <- function(db = db,
                              unknown = unknown, 
                              Th = 0, 
                              allele_freq = T,
                              known_alleles = T) {
  
  if(allele_freq == T){
    
    #~~ Prepare database of allele frequencies
    
    # gather allele sizes into one column
    
    # db_size <- db[[1]] %>%
    #   gather(locus, size, -Allele) %>%
    #   na.omit() %>%
    #   .[c(2,1,3)] %>%
    #   group_by(locus)
    
    # gather allele freqs into one column
    
    # db_freq <- db[[2]] %>%
    #   gather(locus, freq, -Allele) %>%
    #   na.omit() %>%
    #   .[c(2,1,3)] %>%
    #   group_by(locus)
    
     # join dataframes
    # db <- left_join(db_freq, db_size) %>%
    #   mutate(size = as.character(size))
    
    # Genalex:
    
    db <- db[[2]] %>%
      .[-c(1:10,12),-1]
    loci <- c("size", db[1,c(2:ncol(db))])
    
    db <- db %>%
      .[-1,] %>%
      `colnames<-`(unlist(loci)) %>%
      gather(locus, freq, -size) %>%
      drop_na() %>%
      mutate(freq = as.numeric(freq),
             size = as.numeric(size))
    
  } 
  
  if(allele_freq == F){
    
    # Read in genotypes and calculate allele freqs
    
    db <- db %>%
        select(-ID) %>%
        gather(locus, size) %>%
        separate(locus, c("locus", "dip")) %>%
        group_by(locus, size) %>%
        dplyr::summarise(n = n()) %>%
        group_by(locus) %>%
        mutate(freq = n/sum(n))
    
    }
  
  #~~~~~~~~~~~~~~~~~~~#
  #      Unknown      #
  #~~~~~~~~~~~~~~~~~~~#
  
  # If unknown individual's alleles are known, read in unknown
  # Extract allele freqs from db
  
  if(known_alleles == TRUE){
    
    un <- unknown %>%
      gather(locus, size, -ID) %>%
      mutate(size = as.numeric(size)) %>%
      separate(locus, c("locus", "dip")) %>%
      mutate(dip = rep(c("p", "q"), nrow(.)/2)) %>%
      left_join(db, by = c("locus", "size")) %>%
      select(locus, dip, freq) %>%
      spread(dip, freq) %>%
      select(locus, p, q)
    
    #~~~~~~~~~~~~~~~~~~~#
    #        RMP        #
    #~~~~~~~~~~~~~~~~~~~#
      
    # Rob's method (same result at Pmatch(mat))
    
    rmp <- un %>%
      mutate(p_het = (2*(Th+((1-Th)*p))*(Th+((1-Th)*q)))/((1+Th)*(1+(2*Th))),
             p_hom = (((2*Th)+((1-Th)*p))*((3*Th)+((1-Th)*p)))/((1+Th)*(1+(2*Th))),
             prob = case_when(p == q ~ p_hom,
                              TRUE ~ p_het))
    
    cuml_prob <- prod(rmp$prob, na.rm = T)
    
    # Take the reciprocal for the cumulative likelihood
    cuml_likelihood <- 1/cuml_prob
      
    cuml_prob
    return(cuml_prob)
    
    
    }
  
  # If unknown individual's alleles are not known, take two most common as conservative approach

  if(known_alleles == F){
    
    top_2 <- db %>%
      group_by(locus) %>%
      top_n(2, freq) %>%
      arrange(locus, desc(freq)) %>%
      mutate(allele = rep(c("p", "q"))) %>%
      select(locus, allele, freq) %>%
      spread(allele, freq)
    
    #~~~~~~~~~~~~~~~~~~~#
    #        RMP        #
    #~~~~~~~~~~~~~~~~~~~#
    
    rmp <- top_2 %>%
      mutate(p_het = (2*(Th+((1-Th)*p))*(Th+((1-Th)*q)))/((1+Th)*(1+(2*Th))),
             p_hom = (((2*Th)+((1-Th)*p))*((3*Th)+((1-Th)*p)))/((1+Th)*(1+(2*Th))),
             prob = case_when(p == q ~ p_hom,
                              TRUE ~ p_het))
    cuml_prob <- prod(rmp$prob, na.rm = T)
      
    # Take the reciprocal for the cumulative likelihood
    cuml_likelihood <- 1/cuml_prob
    
    cuml_prob
      
    return(cuml_prob)
    
  }
  
}
  

# if(marker == "SNP)

