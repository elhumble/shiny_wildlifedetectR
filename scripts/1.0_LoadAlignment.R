#~~ Check alignment

# Function, where db is the reference database, unkown is the unknown sequence you want, and align is whether you want R to perfrom a sequence alignment

# unknown <- "data/rhino_unknown2.fa"
# db <- "data/rhino_reference2.fa"

# seq1 <- readDNAStringSet(unknown, format="fasta")
# seq2 <- readDNAStringSet(db, format = "fasta")
# seq <- DNAStringSet(c(seq1, seq2))


LoadAlignment <- function(seq, align = FALSE) {
	
	# Reading in the reference and unknown sequence files, then combining them

	# If aligning the data
	if(align == TRUE) {
		
		# Performing the alignment. 
	  # I prefer the 'Muscle' method (not 'ClustalW'), but for some reason this method did not work on a few of my files when testing (it worked for some)
		alignment <- msa::msa(inputSeqs = seq, method = "ClustalW", order = "input")
		
		# Converting the sequences to read into the ape package
		x <- msa::msaConvert(alignment, "phangorn::phyDat")
		x <- ape::as.DNAbin(x)
		
		# Trim the ends of sequences - make both ends blunt (no overhangs)
		x <- ips::trimEnds(as.matrix(x), min.n.seq = length(as.list(x)))
		
		# Then convert
		x <- ape::as.alignment(x)
		x <- ape::as.DNAbin.alignment(x)
		}
	
  # If the data is already aligned
	if(align == FALSE) {
	  
	  # Converting the sequences to read into the ape package
	  x <- ape::as.DNAbin(seq)
	  
	  # Trim the ends of sequences - make both ends blunt (no overhangs)
	  x <- ips::trimEnds(as.matrix(x), min.n.seq = length(as.list(x)))
	  
	  # Then convert
	  x <- ape::as.alignment(x)
	  x <- ape::as.DNAbin.alignment(x)
	}
	
  # Saving the alignment into the global environment for the next species ID functions. 
	# Possibly change this (i.e. assign an ealrier variable e.g. DNAbin) if that's going to be easier to work with. 
  #FinalAlignment <<- x
  # Then do the alignment check
  #par(mfrow=c(1,2))
  #checkAlignment(x, check.gaps = FALSE, plot = TRUE, what = 1)
  #checkAlignment(x, check.gaps = FALSE, plot = TRUE, what = 3)

  return(x)
  
	
}

# NB: Could also look at the alignment with ape:
# par(mar=c(1,1))
# ape::image.DNAbin(x[,ape::seg.sites(x)])
