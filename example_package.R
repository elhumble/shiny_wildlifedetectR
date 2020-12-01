# Example usage of functions

source("scripts/1.0_LoadAlignment.R")
source("scripts/1.1_seqsimilarityID.R")
source("scripts/1.2_phyloID.R")
source("scripts/2.1_individualisation.R")
source("scripts/3.1_parentage.R")
source("scripts/4.1_haplotype_network.R")
#source("scripts/5.1_marker_selection.R")
source("scripts/theme_emily.R")


#~~~~~~~~~~~~~~~~~~~#
#     Species ID    #
#~~~~~~~~~~~~~~~~~~~#

#~~ Rhino example

# Read in files

seq1 <- readDNAStringSet("data/rhino_unknown.fa", format="fasta")
seq2 <- readDNAStringSet("data/rhino_reference.fa", format="fasta")
seq <- DNAStringSet(c(seq1, seq2))

# Alignment, align and trim

x <- LoadAlignment(seq, align = T)
ape::image.DNAbin(x, col = c("#CB2314", "#FAD510", "#4daf4a", "#5BBCD6", "grey", "black"))

# Phylogeny

tree <- phyloID(DNAbin = x, 
                model = "auto", 
                outgroup = T, 
                gamma = FALSE, 
                bootstraps = 10)

plotPhylo(tree)

# Seq Similarity

seqsimilarityID(DNAbin = x, 
                model = "K80", 
                hits = 10)

# Macaw example

seq1 <- readDNAStringSet("data/macaw_unknown.fa", format="fasta")
seq2 <- readDNAStringSet("data/macaw_reference.fa", format="fasta")
seq <- DNAStringSet(c(seq1, seq2))


# Alignment, align and trim

x <- LoadAlignment(seq, align = F)
ape::image.DNAbin(x, col = c("#CB2314", "#FAD510", "#4daf4a", "#5BBCD6", "grey", "black"))

# Phylogeny

tree <- phyloID(DNAbin = x, 
                model = "auto", 
                outgroup = T, 
                gamma = FALSE, 
                bootstraps = 10)

plotPhylo(tree)

# Seq Similarity

seqsimilarityID(DNAbin = x, 
                model = "K80", 
                hits = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Individualisation    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Allele frequencies

# Badger example

u <- read_excel("data/badger_trace.xlsx")
path <- "data/badger_reference_db.xlsx"

db <- path %>%
  excel_sheets() %>%
  set_names() %>%
  purrr::map(read_excel, path = path)

# RMP

individualisation(db = db,
                  unknown = u,
                  Th = 0.12,
                  allele_freq = T,
                  known_alleles = T)

#~~ Genotype table example

u <- read_excel("data/individualisation_unknown.xls")
db <- read_excel("data/individualisation_genotype_db.xls")

# Likelihood

1/individualisation(db = db,
                  unknown = u,
                  Th = 0,
                  allele_freq = F,
                  known_alleles = T)



#~~~~~~~~~~~~~~~~~~~#
#     Parentage     #
#~~~~~~~~~~~~~~~~~~~#

#~~ SNPs

# Testing with fur seal genotypes
# Important that locus names (columns) are exactly the same for both the unknown and db files

#~~ Test 1

db <- fread("data/fur_seal_SNP_genotypes.raw", header = T)
trio <- fread("data/fur_seal_trio.raw")

# Get allele frequencies from plink db file

freqs <- alleleFreqs_plink(db)

# Prepare trio file

trio <- prep_snp_trio(trio)

# Join data

df <- dplyr::left_join(trio, freqs, by = "SNP_ID")

# Run test

test1_snps(df)

#~~ Test 2

test2_snps(df)

#~~ Test 3

duo <- fread("data/fur_seal_duo.raw")

# Prepare duo file

duo <- prep_snp_duo(duo)
df <- left_join(duo, freqs, by = "SNP_ID")

test3_snps(df)

#~~ Microsatellites

#~~ Test 1: Joint parentage

# Missing puppy (unknown allele sizes)

unknown <- read_excel("data/dog_trio.xlsx")

# Genalex

path <- "data/dog_reference_db.xlsx"

db <- path %>%
   excel_sheets() %>%
   set_names() %>%
   purrr::map(read_excel, path = path)

freqs <- freqs2freqs_genalex(db)


# Get top two most frequent alleles

freqs <- top2freq(freqs)

# Recode family genotypes

unknown <- prep_ms_trio(unknown)
recode <- numericSNPmat(unknown)

# Prepare df for LRT

df_ms <- data.frame(t(recode)) %>%
  `colnames<-`(c("V1", "V2", "V3")) %>%
  mutate(Locus = rownames(.)) %>%
  left_join(freqs, by = "Locus") %>%
  mutate(p = as.numeric(p),
         q = as.numeric(q))

# Run LRT

case1 <- test1_ms(df = df_ms, Th = 0.11)

1/(prod(case1$prob, na.rm = T))

#~~ Test 2: Paternity

# African grey parrot

unknown <- read_excel("data/AGP_trio.xlsx")

path <- "data/AGP_reference_db.xlsx"

db <- path %>%
  excel_sheets() %>%
  set_names() %>%
  purrr::map(read_excel, path = path)

# Genalex
freqs <- freqs2freqs_genalex(db)

# Extract allele frequencies for unknown

freqs <- extract_known_alleles(unknown, freqs)

# Recode family genotypes

unknown <- prep_ms_trio(unknown)
recode <- numericSNPmat(unknown)

# Prepare df for LRT

df_ms <- data.frame(t(recode)) %>%
  `colnames<-`(c("V1", "V2", "V3")) %>%
  mutate(Locus = rownames(.)) %>%
  left_join(freqs, by = "Locus") %>%
  mutate(p = as.numeric(p),
         q = as.numeric(q))

# Run LRT

case2 <- test2_ms(df = df_ms, Th = 0.01)

prod(case2$prob, na.rm = T)

#~~ Test 3: Motherless paternity

# Motherless paternity

unknown <- read_excel("data/tiger_duo.xlsx")

path <- "data/tiger_reference_db.xlsx"

db <- path %>%
  excel_sheets() %>%
  set_names() %>%
  purrr::map(read_excel, path = path)

# Genalex
freqs <- freqs2freqs_genalex(db)

# Extract allele frequencies for unknown

freqs <- extract_known_alleles(unknown, freqs)

# Recode family genotypes

unknown <- prep_ms_duo(unknown)
recode <- numericSNPmat(unknown)

# Prepare df for LRT

df_ms <- data.frame(t(recode)) %>%
  `colnames<-`(c("V1", "V2")) %>%
  mutate(Locus = rownames(.)) %>%
  left_join(freqs, by = "Locus") %>%
  mutate(p = as.numeric(p),
         q = as.numeric(q))

# Run LRT

case3 <- test3_ms(df = df_ms, Th = 0.01)
prod(case3$prob, na.rm = T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Haplotype Network      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

seqs <- readDNAStringSet("data/rtb_cockatoo.fa", format = "fasta")

haplotype_network(seqs, pop_header = T, 
                  scale_ratio = 0.3,
                  x_coord = -20,
                  y_coord = 4)

