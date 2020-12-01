# library(pegas)
# library(ape)
# library(msa)
# library(RColorBrewer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#     Haplotype Network Analysis      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#database <- readDNAStringSet("data/rhino_hap.fa", format = "fasta")

haplotype_network <- function(database = database,
                              pop_header = F,
                              pop_file = pop_file,
                              scale_ratio = scale_ratio,
                              x_coord = x_coord,
                              y_coord = y_coord){
  
  
  #~~ Read in fasta files: I've just coded this so the user has to upload two separate files
  
  # d <- readDNAStringSet(database, format="fasta")
  # unknown <- readDNAStringSet(unknown, format = "fasta")
  #seqs <- DNAStringSet(c(database, unknown))
  #seqs <- DNAStringSet(database)
  seqs <- ape::as.DNAbin(database)
  
  #~~ Replace sample names with pop IDs
  
  if(pop_header == TRUE) {
    
    names(seqs) <- sub('.*(_|-)', '', names(seqs)) # Probably need something different here
    
  }
  
  if(pop_header == FALSE) {
    
    # Read in csv files with pop IDs in the order of the associated individual in the fasta file, then replace.
    #pops <- read.csv(pop_file, header = FALSE)
    pops <- as.matrix(pop_file)
    names(seqs) <- pops
    
  }
  
  
  #~~ Compute haplotype network
  
  # Do we need this????
  
  # Compute distance matrix. By default this uses the "K80" evolutioanry model, and gamma is false. Can easily change this and/or add a model testing step.
  #e <- dist.dna(d)
  
  # Extract haplotypes from the list of sequences.
  h <- pegas::haplotype(seqs)
  
  # Compute haplotype network. Can include Templeton's probabilities if needed.
  net <- pegas::haploNet(h)
  
  #~~ Plot haplotype network
  
  # Label haplotype groups (have just included the first here as it is clear from the colours which group the unknown is in)
  # Although do we even need this: size of circles is proportional
  
  lab <- as.character(attr(net, "freq"))
  attr(net, "labels") <- lab
  
  ind.hap<-with(
    utils::stack(setNames(attr(h, "index"), rownames(h))),
    table(hap=ind, pop=names(seqs)[values])
  )
  
  
  # Code to standardise values by total number of seqs -- hope this works.
  
  ind.hap <- apply(ind.hap, 1, function(x) x/sum(colSums(ind.hap)))
  ind.hap <- t(ind.hap)
  
  # Plotting network
  
  #pal <- brewer.pal(n = 9, name = "Set1")
  #pal <- c(colorRampPalette(brewer.pal(10,"Paired"))(50),colorRampPalette(brewer.pal(10,"Set3"))(50))
  #set.seed(150)
  #pal <- sample(pal) # get a good sample
  #pal <- brewer.pal.info[brewer.pal.info$category == "Set3"]
  #pals_vector<-unlist(mapply(brewer.pal, pal$maxcolors, rownames(pal)))
  #pal <- sample(pal)
  
  
  # 18 colours
  pal <- c("cadetblue3", "chartreuse1", "chocolate1","darkmagenta","darkgreen", 
           "firebrick","darkgoldenrod1","bisque4","blue2", "deeppink", "saddlebrown", 
           "yellow","gray75","red", "seagreen", "cyan", "plum2","gray15")
  
  #pie(rep(1, length(pal)), col = pal , main="")
  
  
  plot(net, size=attr(net, "freq"), bg=pal, scale.ratio=scale_ratio, pie=ind.hap)
  
  # Adding legend
  
  legend(x_coord, y_coord, colnames(ind.hap), col=pal, pch=19, ncol=2)
  
  
}

