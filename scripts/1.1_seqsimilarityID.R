#~~ Get seq similarity table


#DNAbin <- x

seqsimilarityID <- function(DNAbin = FinalAlignment, model = model, hits = hits) {
  
  # Constructing a distance matrix using the alignment from the LoadAlignment.R script. This calculates a proportion of sites that differ between each pair.
  dist.matrix <- ape::dist.dna(DNAbin, model = model, as.matrix = TRUE)
  
  # Convert it to a dataframe
  dist.df <- as.data.frame(dist.matrix)
  
  # Now we can order the data frame by the distance from the unknown sequence (in ascending order). NB: the unknown sequence is the first one.
  ordered.df <- dist.df[order(dist.df[1]), ]
  # This is the column of ordered distances
  res.distances <- ordered.df[1]
  
  ### Now making a table with sequences names, sequence lengths, sequence similarites and number of differences
  ## Sequences names
  Name <- as.data.frame(rownames(ordered.df))
  
  ## Sequence lenghts - they should all be the same length, but it might be good to include anyway.
  # Number of sequences
  n.sequences <- length(as.list(DNAbin))
  # Sequence length
  seq.length <- dim(as.matrix(DNAbin))[2]
  # List of lengths
  Length <- as.data.frame(rep.int(seq.length, n.sequences))
  
  ## Sequence similarities - as a %
  Similarity <- (1-(res.distances))*100
  
  ## Number of differences
  Differences <- res.distances*seq.length
  
  ### Combine these 
  #Table <- cbind(Name, Length, Differences, Similarity)
  Table <- cbind(Name, Length, Similarity)
  
  # Rename columns and rows
  #colnames(Table) <- c("Names","Length","# of differences", "Similarity (%)")
  colnames(Table) <- c("Names","Length", "Similarity (%)")
  
  rownames(Table) <- NULL
  # Get rid of differences and similarity of the unknown sequence to itself.
  #Table$`# of differences`[1] <- "N/A"
  Table$`Similarity (%)`[1] <- "N/A"
  
  ## How many sequences to report specified by 'hits' (relevant if it's a massive database)
  if(hits=="all"){
    return(Table)
  } else{
    return(Table[1:(as.numeric(hits)+1),])
  }

}
