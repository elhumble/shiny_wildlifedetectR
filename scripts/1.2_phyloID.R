#~~ Generate and plot phylogeny

# Specify if their is an outgroup species (at the end of the sequences).
# NB: I probably don't even need 'DNAbin = FinalAlignment' - 
# I could probably just use the 'FinalAlignment' object throughout because it's already in the environment from the LoadAlignment.R function. 

phyloID <- function(DNAbin = FinalAlignment, model = "auto", outgroup = TRUE, gamma = FALSE, bootstraps = bootstraps) {
  
  if(model == "auto"){
    
    # Converting the DNAbin object (from the LoadAlignment.R script) into phyDat
    phy <- phangorn::phyDat(DNAbin, type = "DNA")
    
    # Run the model testing
    # Should we perform this model test with the outgorup included???
    # NB: this modelTest function may not work within a GUI interface. We may need to use the "capture.output" funciton to stop it printing out.
    # NB: I = FALSE, because the distance function to construct the phylogenetic tree does not include invariable sites modeling.
    allmodels <- phangorn::modelTest(phy, model = "all", G = TRUE, I = FALSE, tree = NULL, multicore = FALSE)
    
    # Using only the available models when making a neighbour-joining tree as not all models that were tested are available to make a tree.
    available <- c("JC", "JC+G","K80", "K80+G", "F81", "F81+G","K81", "K81+G")
    availablemodels <- allmodels[allmodels$Model %in% available,]
    #NB: is JC69 the same as JC??? Check the other models in '?dist.dna' to see if they are called other names and hence still applicable.
    
    # Selecting the model with the lowest BIC
    optimalmodel <- availablemodels[which.min(availablemodels$BIC),1]
    
    # Selecting the nucleotide substitution portion
    model <- strsplit(optimalmodel, "[+]")[[1]][1]
    
    # Indicating whether we are using gamma distributed rates
    gamma <- strsplit(optimalmodel, "[+]")[[1]][2] == "G"
  }
  
  # Now to make a neighbour-joining tree
  
  # Making a distance matrix
  d <- ape::dist.dna(DNAbin, model = model, gamma = gamma)
  
  # Now to make the tree
  tree <- ape::bionj(d)
  
  # If they use an outgroup
  if(outgroup == TRUE){
    
    # Selecting the outgroup sequence. NB: make sure their outgroup sequence is last.
    seqnames <- names(as.list(DNAbin))
    # Number of sequences
    n <- length(seqnames)
    # Select the last sequence (i.e. the outgroup)
    outgroup <- seqnames[n]
    
    # Now root the tree with the outgroup
    tree <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)
    
    # Bootstrapping
    boot <- ape::boot.phylo(phy = tree, DNAbin, B = bootstraps,  function(d) root(bionj(dist.dna(d, model = model, gamma = gamma)), outgroup = outgroup))
    # NB: again, this prints to the console so may need to use the "capture.output" funciton
    
    # Now plot the tree with bootstrap support as a percentage
    #plot(tree, cex=0.5, main = "Rooted Tree"); nodelabels((boot/bootstraps), cex = 0.5)
    
    tree$node.labels <- boot/bootstraps
    return(tree)
    # ggtree
    
    # node_df <- data.frame(1:nrow(tree$edge))
    # colour <- c(rep("cornflowerblue", 1), rep("gray55", nrow(node_df)-1))
    # node_df$colour <- colour
    # colnames(node_df) <- c("node", "colour")
    # 
    # # boot <- c(boot/bootstraps)
    # 
    # p <- ggtree(tree, aes(color = I(colour))) %<+% node_df +
    #   geom_tiplab(size=4) + 
    #   theme(legend.position="right") 
    # 
    # d <- p$data
    # d <- d[!d$isTip,]
    # d$label <- as.numeric(d$label)
    # #d <- d[d$label > 80,]
    # 
    # p <- p + geom_text(data = d, aes(x = branch, label=label), size = 3, vjust=-.5, color='grey20')
    # 
    # print(p)
    
    
  }
  
  # NEED TO UPDATE PLOTTING BIT ON THIS!!!!!!!
  
  
  # If they don't use an outgroup
  if(outgroup == FALSE){
    
    # Bootstrapping
    boot <- ape::boot.phylo(phy = tree, DNAbin, B = bootstraps,  FUN = function(d) bionj(dist.dna(d, model = model, gamma = gamma)))
    # NB: again, this prints to the console so may need to use the "capture.output" funciton
    
    # Now plot the tree with bootstrap support as a percentage
    #plot(tree, cex=0.5, main = "Unrooted Tree"); nodelabels((boot/bootstraps), cex = 0.5)
    
    tree$node.labels <- boot/bootstraps
    
    return(tree)
    
    # ggtree
    
    # node_df <- data.frame(1:nrow(tree$edge))
    # colour <- c(rep("cornflowerblue", 1), rep("gray55", nrow(node_df)-1))
    # node_df$colour <- colour
    # colnames(node_df) <- c("node", "colour")
    # 
    # # boot <- c(boot/bootstraps)
    # 
    # p <- ggtree(tree, aes(color = I(colour))) %<+% node_df +
    #   geom_tiplab(size=4) + 
    #   theme(legend.position="right") 
    # 
    # d <- p$data
    # d <- d[!d$isTip,]
    # d$label <- as.numeric(d$label)
    # #d <- d[d$label > 80,]
    # 
    # p + geom_text(data = d, aes(x = branch, label=label), size = 3, vjust=-.5, color='grey20')

  }
  
}


plotPhylo <- function(tree){
  
  cols <- c("#5BBCD6", rep("grey20", length(tree$tip.label)-1))
  plot(tree, lwd = 1, edge.color="grey20", tip.color = cols, use.edge.length = FALSE)
  i <- c(1)
  j <-c(2:length(tree$tip.label)-1)
  #tiplabels(tree$tip.label[i], i, adj = 0, frame = "r", bg = "cornflowerblue")
  
}


# plotPhylo <- function(tree){
#   
#   node_df <- data.frame(1:nrow(tree$edge))
#   colour <- c(rep("cornflowerblue", 1), rep("gray55", nrow(node_df)-1))
#   node_df$colour <- colour
#   colnames(node_df) <- c("node", "colour")
#   
#   # boot <- c(boot/bootstraps)
#   
#   p <- ggtree(tree, aes(color = I(colour))) %<+% node_df +
#     geom_tiplab(size=4) + 
#     theme(legend.position="right") 
#   
#   d <- p$data
#   d <- d[!d$isTip,]
#   d$label <- as.numeric(d$label)
#   #d <- d[d$label > 80,]
#   
#   p + geom_text(data = d, aes(x = branch, label=label), 
#                 size = 3, vjust=-.5, color='grey20') +
#     ggplot2::xlim(0, 0.08)
#   
# }



# Could make the tree look alot nicer
# E.g. colour the unknown sequence and the outgroup
# It may be better to ladderize the tree before plotting with 'ladderize(tree)'
          