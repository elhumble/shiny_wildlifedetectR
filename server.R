
#~~ Server functions

shinyServer <- 
  function(input, output) {
    
    #~~~~~~~~~~~~~~~~~~~~~~~~#
    #       Species ID       #
    #~~~~~~~~~~~~~~~~~~~~~~~~#    
    
    #~~ Phylogeny
    
    observeEvent(input$SpID_action1,{
      
      # Load fasta files
      req(input$SpID_unknown)
      seq1 <- readDNAStringSet(input$SpID_unknown$datapath, format="fasta")
      
      req(input$SpID_ref)
      seq2 <- readDNAStringSet(input$SpID_ref$datapath, format="fasta")
      
      seq <- DNAStringSet(c(seq1, seq2))
      
      # Alignment
      
      if(input$SpID_align == T){
        x <- LoadAlignment(seq, align = T)
      }
      
      if(input$SpID_align == F){
        x <- LoadAlignment(seq, align = F)
      }
      
      # Plot Alignment
      
      output$SpIDAlign <- renderPlot({
        par(mar=c(5,20,4,2)+0.1) # could do something here with counting characters in the IDs
        ape::image.DNAbin(x, col = c("#CB2314", "#FAD510", "#4daf4a", "#5BBCD6", "grey", "black"))
        #ape::image.DNAbin(x)
        #checkAlignment(x, check.gaps = FALSE, plot = TRUE, what = 1) 
      })
      
      # Phylogeny
      
      output$SpIDPhylo <- renderPlot({
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Running bootstraps", value = 0.3)
        
        tree <- phyloID(DNAbin = x, 
                        model = "auto", 
                        outgroup = input$SpID_PhyloOut, 
                        gamma = FALSE, 
                        bootstraps = as.numeric(input$SpID_PhyloBoot))
        
        progress$inc(1, detail = paste("Bootstrapping complete"))
        
        plotPhylo(tree)
        
      })
      
    })
    
    #~~ Sequence similarity
    
    observeEvent(input$SpID_action2,{
      
      # Load data
      
      req(input$SpID_unknown)
      seq1 <- readDNAStringSet(input$SpID_unknown$datapath, format="fasta")
      
      req(input$SpID_ref)
      seq2 <- readDNAStringSet(input$SpID_ref$datapath, format="fasta")
      
      seq <- DNAStringSet(c(seq1, seq2))
      
      #Alignment
      
      if(input$SpID_align == T){
        x <- LoadAlignment(seq, align = T)
      }
      
      if(input$SpID_align == F){
        x <- LoadAlignment(seq, align = F)
      }
      
      # Plot Alignment
      
      output$SpIDAlign <- renderPlot({
        par(mar=c(5,20,4,2)+0.1) # could do something here with counting characters in the IDs
        ape::image.DNAbin(x, col = c("#CB2314", "#FAD510", "#4daf4a", "#5BBCD6", "grey", "black"))
        #checkAlignment(x, check.gaps = FALSE, plot = TRUE, what = 1) 
        
      })
      
      # Sequence Similarity Matrix
      
      output$SpIDSeqSim <- renderTable({
        seqsimilarityID(DNAbin = x, 
                        model = input$SpID_SeqSimModel, 
                        hits = input$SpID_SeqSimHits)
      })
      
    })
    
    #~~ Output report phylo
    
    output$SpID_report_phylo <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "sp_ID_report_phylo.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "sp_ID_report_phylo.Rmd")
        file.copy("sp_ID_report_phylo.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        
        req(input$SpID_unknown)
        seq1 <- readDNAStringSet(input$SpID_unknown$datapath, format="fasta")
        req(input$SpID_ref)
        seq2 <- readDNAStringSet(input$SpID_ref$datapath, format="fasta")
        
        params <- list(seq1 = seq1,
                       seq2 = seq2,
                       unk = input$SpID_unknown$name,
                       ref = input$SpID_ref$name,
                       align = input$SpID_align,
                       outgroup = input$SpID_PhyloOut,
                       boot = input$SpID_PhyloBoot)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Generating report", value = 0.7)
        
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
        progress$inc(1, detail = paste("Report generated"))
        
      }
    )
    
    
    #~~ Output report seq sim
    
    output$SpID_report_seqsim <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "sp_ID_report_seqsim.html", #SpID_report_phylo.pdf
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "sp_ID_report_seqsim.Rmd")
        file.copy("sp_ID_report_seqsim.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        
        req(input$SpID_unknown)
        seq1 <- readDNAStringSet(input$SpID_unknown$datapath, format="fasta")
        req(input$SpID_ref)
        seq2 <- readDNAStringSet(input$SpID_ref$datapath, format="fasta")
        
        params <- list(seq1 = seq1,
                       seq2 = seq2,
                       unk = input$SpID_unknown$name,
                       ref = input$SpID_ref$name,
                       align = input$SpID_align,
                       model = input$SpID_SeqSimModel,
                       hits = input$SpID_SeqSimHits)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Generating report", value = 0.7)
        
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
        progress$inc(1, detail = paste("Report generated"))
        
      }
    )
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~#
    #    Individualization   #
    #~~~~~~~~~~~~~~~~~~~~~~~~#    
    
    observeEvent(input$action2, {
      
      if(input$method2 == F){
        
        output$Individualisation <- renderUI({
          
          # Read in input files ready for function
          
          req(input$UnknownProfile)
          u <- read_excel(input$UnknownProfile$datapath)
          
          req(input$Ref_db)
          db <- read_excel(input$Ref_db$datapath)
          
          if(input$method10 == "Random match probability"){
            
            str0 <- paste0("Match: Random match probability is ", "<b>",
            individualisation(db = db,
                              unknown = u,
                              Th = input$IndividualizationTheta,
                              allele_freq = input$method2,
                              known_alleles = input$method11),"</b>")
            
            str1 <- paste0("The probability of a match if this genotype 
                           derived from an individual other than the suspect is ", "<b>",
                           individualisation(db = db,
                                             unknown = u,
                                             Th = input$IndividualizationTheta,
                                             allele_freq = input$method2,
                                             known_alleles = input$method11),"</b>")
            
            #HTML(paste(str0, str1, sep = '<br/>'))
          }
          
          if(input$method10 == "Likelihood ratio approach"){
            
            str0 <- paste0("Match: Cumulative likelihood is ", "<b>",
                           1/individualisation(db = db,
                                             unknown = u,
                                             Th = input$IndividualizationTheta,
                                             allele_freq = input$method2,
                                             known_alleles = input$method11),"</b>")
            
            str1 <- paste0("It is ", "<b>",1/individualisation(db = db,
                                                               unknown = u,
                                                               Th = input$IndividualizationTheta,
                                                               allele_freq = input$method2,
                                                               known_alleles = input$method11),"</b>", 
                           " times more likely to obtain this DNA profile if it originated from the suspect than 
                 from an unknown unrelated individual in the population.")
            
            # HTML(paste(str0, str1, sep = '<br/>'))
            
          }
          
          HTML(paste(str0, str1, sep = '<br/>'))
          
          
        })
      }
      
      if(input$method2 == T){
        
        output$Individualisation <- renderUI({
          
          # Read in input files ready for function
          
          req(input$UnknownProfile)
          u <- read_excel(input$UnknownProfile$datapath)
          
          req(input$Ref_db)
          path <- input$Ref_db$datapath
          
          db <- path %>%
            excel_sheets() %>%
            set_names() %>%
            purrr::map(read_excel, path = path)
          
          
          if(input$method10 == "Random match probability"){
            
            str0 <- paste0("Match: Random match probability is ", "<b>",
                           individualisation(db = db,
                                             unknown = u,
                                             Th = input$IndividualizationTheta,
                                             allele_freq = input$method2,
                                             known_alleles = input$method11),"</b>")
            
            str1 <- paste0("The probability of a match if this genotype 
                           derived from an individual other than the suspect is ", "<b>",
                           individualisation(db = db,
                                             unknown = u,
                                             Th = input$IndividualizationTheta,
                                             allele_freq = input$method2,
                                             known_alleles = input$method11),"</b>")
            #HTML(paste(str0, str1, sep = '<br/>'))
          }
          
          if(input$method10 == "Likelihood ratio approach"){
            
            str0 <- paste0("Match: Cumulative likelihood is ", "<b>",
                           individualisation(db = db,
                                             unknown = u,
                                             Th = input$IndividualizationTheta,
                                             allele_freq = input$method2,
                                             known_alleles = input$method11),"</b>")
            
            str1 <- paste0("It is ", "<b>",1/individualisation(db = db,
                                                               unknown = u,
                                                               Th = input$IndividualizationTheta,
                                                               allele_freq = input$method2,
                                                               known_alleles = input$method11),"</b>", 
                           " times more likely to obtain this DNA profile if it originated from the suspect than 
                 from an unknown unrelated individual in the population.")
            
            # HTML(paste(str0, str1, sep = '<br/>'))
            
          }
          
          HTML(paste(str0, str1, sep = '<br/>'))
          
          
        })
      }
    })
    
    #~~ Output report individualisation
    
    output$ind_report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "ind_report.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "ind_report.Rmd")
        file.copy("ind_report.Rmd", tempReport, overwrite = TRUE)
        
        
        # Prep df for test
        if(input$method2 == T){
          
          req(input$UnknownProfile)
          u <- read_excel(input$UnknownProfile$datapath)
          
          req(input$Ref_db)
          path <- input$Ref_db$datapath
          
          db <- path %>%
            excel_sheets() %>%
            set_names() %>%
            purrr::map(read_excel, path = path)
          
        }
        
        if(input$method2 == F){
          
          req(input$UnknownProfile)
          u <- read_excel(input$UnknownProfile$datapath)
          
          req(input$Ref_db)
          db <- read_excel(input$Ref_db$datapath)
          
        }
  
        # Set up parameters to pass to Rmd document
        
        params <- list(df = db,
                       trace = u,
                       db = input$Ref_db$name,
                       unk = input$UnknownProfile$name,
                       method = input$method10,
                       Th = input$IndividualizationTheta,
                       allele_freq = input$method2,
                       known_alleles = input$method11)

        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Generating report", value = 0.7)
        
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
        progress$inc(1, detail = paste("Report generated"))
        
      }
    )
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~#
    #       Parentage        #
    #~~~~~~~~~~~~~~~~~~~~~~~~#   
    
    observeEvent(input$action5, {
      
      if (input$ParentageMarker == "Microsatellites") {
        
        req(input$ParentageProfile)
        fam <- read_excel(input$ParentageProfile$datapath)
        
        if (input$ParentageData == "Allele frequency table (Genalex format)") {
          # Load freq table
          req(input$ParentageDB)
          path <- input$ParentageDB$datapath
          
          db <- path %>%
            excel_sheets() %>%
            set_names() %>%
            purrr::map(read_excel, path = path)
          
          # Wrangle allele frequencies
          freqs <- freqs2freqs_genalex(db)
        }
        
        else {
          # Load genos
          req(input$ParentageDB)
          db <- read_excel(input$ParentageDB$datapath)
          
          # Freqs from genos
          freqs <- genos2freqs(db)
        }
        
        # Extract allele frequencies
        
        if (input$ParentageAlleles == T) {
          freqs <- extract_known_alleles(fam, freqs)
        }
        else {
          freqs <- top2freq(freqs)
        }
        
        
        # Recode family genotypes and prepare df for LRT
        
        if (input$ParentageTest == "Probability of motherless paternity") {
          fam <- prep_ms_duo(fam)
          
          recode <- numericSNPmat(fam)
          
          df_ms <- data.frame(t(recode)) %>%
            `colnames<-`(c("V1", "V2")) %>%
            mutate(Locus = rownames(.)) %>%
            left_join(freqs, by = "Locus") %>%
            mutate(p = as.numeric(p),
                   q = as.numeric(q))
        }
        else {
          
          fam <- prep_ms_trio(fam)
          
          recode <- numericSNPmat(fam)
          
          df_ms <- data.frame(t(recode)) %>%
            `colnames<-`(c("V1", "V2", "V3")) %>%
            mutate(Locus = rownames(.)) %>%
            left_join(freqs, by = "Locus") %>%
            mutate(p = as.numeric(p),
                   q = as.numeric(q))
        }
        
        #~~ Run LRT
        
        if (input$ParentageTest == "Probability of joint parentage"){
          out <- test1_ms(df = df_ms, Th = input$ParentageTheta)
        }
        
        if (input$ParentageTest == "Probability of paternity (single parentage)"){
          out <- test2_ms(df = df_ms, Th = input$ParentageTheta)
        }
        
        if (input$ParentageTest == "Probability of motherless paternity"){
          out <- test3_ms(df = df_ms, Th = input$ParentageTheta)
        }
        
        #out <- prod(out$prob, na.rm = T)
        
        
        #~~ Output statements
        
        if (input$ParentageTest == "Probability of joint parentage"){
          
          if('Y' %in% out$exclude)  {
            
            output$Parentage_statement <- renderUI({
              
              # Print how many loci are 'causing' the exclusion, and list them.
              # Position of loci
              loci_pos <- which('Y' == out$exclude)
              # Name of loci
              loci_name <- out[loci_pos,]
              loci_name <- loci_name$Locus
              
              # Print statement, then number of loci 'causing' the exclusion. 
              #Note, I think there is better teminology we can use than 'causing' an exclusion!
              str0 <- paste0(" ")
              str1 <- paste0("Exclusion. Number of loci causing the exclusion: ", length(loci_name),'<br/>')
              str2 <- paste0("Loci that caused the exclusion:")
              str3 <- paste(loci_name, collapse=', ')
              
              # Incorporating a mutation rate:
              mutation_rate <- (10^(-3))
              str4 <- paste0('<br/>',"Chance of seeing ", length(loci_name), " mismatches based on the mutation rate: ",
                             (1-((1-((mutation_rate)^(length(loci_name))))^(length(out))))*100, "%")

              HTML(paste(str0, str1, str2, str3, str4, sep = '<br/>'))
            })
            
          } else {
            
            output$Parentage_statement <- renderUI({
              str0 <- paste0(" ")
              str1 <- "Match:"
              str2 <- paste0("It is ", "<b>",1/(prod(out$prob, na.rm = T)),"</b>", " times more likely to see an offspring with this genotype if the suspects are
                     the biological parents, rather than if two other random individuals are the parents.")
              HTML(paste(str0, str1, str2, sep = '<br/>'))
              
            })
          }
        }
        
        if (input$ParentageTest == "Probability of paternity (single parentage)"){
          if('Y' %in% out$exclude)  {
            
            output$Parentage_statement <- renderUI({
              
              # Print how many loci are 'causing' the exclusion, and list them.
              # Position of loci
              loci_pos <- which('Y' == out$exclude)
              # Name of loci
              loci_name <- out[loci_pos,]
              loci_name <- loci_name$Locus
              
              # Print statement, then number of loci 'causing' the exclusion. 
              #Note, I think there is better teminology we can use than 'causing' an exclusion!
              str0 <- paste0(" ")
              str1 <- paste0("Exclusion. Number of loci causing the exclusion: ", length(loci_name),'<br/>')
              str2 <- paste0("Loci that caused the exclusion:")
              str3 <- paste(loci_name, collapse=', ')
              
              # Incorporating a mutation rate:
              mutation_rate <- (10^(-3))
              str4 <- paste0('<br/>',"Chance of seeing ", length(loci_name), " mismatches based on the mutation rate: ",
                              (1-((1-((mutation_rate)^(length(loci_name))))^(length(out))))*100, "%")
            
              HTML(paste(str0, str1, str2, str3, str4, sep = '<br/>'))
            })
            
          } else {
            
            output$Parentage_statement <- renderUI({
              str0 <- paste0(" ")
              str1 <- "Match:"
              str2 <- paste0("It is ", "<b>",(prod(out$prob, na.rm = T)),"</b>", " times more likely to see an offspring with this genotype if the 
                   suspect is the biological parent, rather than if a random individual is the parent.")
              HTML(paste(str0, str1, str2, sep = '<br/>'))
              
            })
          }
          
        }
        
        if (input$ParentageTest == "Probability of motherless paternity"){
          if('Y' %in% out$exclude)  {
            
            output$Parentage_statement <- renderUI({
              
              # Print how many loci are 'causing' the exclusion, and list them.
              # Position of loci
              loci_pos <- which('Y' == out$exclude)
              # Name of loci
              loci_name <- out[loci_pos,]
              loci_name <- loci_name$Locus
              
              # Print statement, then number of loci 'causing' the exclusion. 
              #Note, I think there is better teminology we can use than 'causing' an exclusion!
              str0 <- paste0(" ")
              str1 <- paste0("Exclusion. Number of loci causing the exclusion: ", length(loci_name),'<br/>')
              str2 <- paste0("Loci that caused the exclusion:")
              str3 <- paste(loci_name, collapse=', ')
              
              # Incorporating a mutation rate:
              mutation_rate <- (10^(-3))
              str4 <- paste0('<br/>',"Chance of seeing ", length(loci_name), " mismatches based on the mutation rate: ",
                             (1-((1-((mutation_rate)^(length(loci_name))))^(length(out))))*100, "%")
              
              HTML(paste(str0, str1, str2, str3, str4, sep = '<br/>'))
            })
            
          } else {
            
            output$Parentage_statement <- renderUI({
              str0 <- paste0(" ")
              str1 <- "Match:"
              str2 <- paste0("It is ", "<b>",(prod(out$prob, na.rm = T)),"</b>", " times more likely to see an offspring with this genotype if the suspect 
                   is the biological parent, rather than if a random individual is the parent.")
              HTML(paste(str0, str1, str2, sep = '<br/>'))
              
            })
          }
        }
        
      }
      
      #~~ SNPs
      
      if(input$ParentageMarker == "SNPs"){
      }
      
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~#
    #    Haplotype Network   #
    #~~~~~~~~~~~~~~~~~~~~~~~~#
    
    observeEvent(input$action4, {
      
      output$hapNet <- renderPlot({
        
        if(is.null(input$file_input))
        {
          #req(input$UnknownPopSeqs)
          #u_seq <- readDNAStringSet(input$UnknownPopSeqs$datapath, format = "fasta")
          
          req(input$ReferencePopSeqs)
          ref_seqs <- readDNAStringSet(input$ReferencePopSeqs$datapath, format = "fasta")
          
          haplotype_network(database = ref_seqs,
                            #unknown = u_seq,
                            pop_header = T,
                            scale_ratio = input$HapScaleRatio,
                            x_coord = input$HapX,
                            y_coord = input$HapY)
        }
        else
        {
          # req(input$UnknownPopSeqs)
          #u_seq <- readDNAStringSet(input$UnknownPopSeqs$datapath, format = "fasta")
          
          req(input$ReferencePopSeqs)
          ref_seqs <- readDNAStringSet(input$ReferencePopSeqs$datapath, format = "fasta")
          
          req(input$PopList)
          pop_file <- read.csv(input$PopList$datapath, header = FALSE)
          
          haplotype_network(database = ref_seqs,
                            # unknown = u_seq,
                            pop_header = F,
                            scale_ratio = input$HapScaleRatio,
                            x_coord = input$HapX,
                            y_coord = input$HapY)
        }
        
        
      })
      
    })
    

    #~~~~~~~~~~~~~~~~~~~~~~~#
    #    Marker Selection   #
    #~~~~~~~~~~~~~~~~~~~~~~~#
    
    
    observeEvent(input$action6, {
      #This function is repsonsible for loading in the selected file
      infile <- input$SNPdata
      if (is.null(infile)) {
        # User has not uploaded a file yet
        return(NULL)
      }
      
      input_data <- read.PLINK(input$SNPdata$datapath,
                               parallel = F) # parallel = T not supported on win
      #input_data <- read.PLINK("data/ORYX_500K_ld.raw")
      
      x <- as.matrix(input_data)
      df <- as.data.frame(t(x))
      
      df <- df[c(1:100),]
      SNP_ID <- rownames(df)
      #df[df=="2"]<-NA
      
      maf <- df %>% mutate(n0 = rowSums(. == 0, na.rm = T),
                           n1 = rowSums(. == 1, na.rm = T), # calc_n
                           n2 = rowSums(. == 2, na.rm = T),
                           n = n0 + n1 + n2,
                           p =  ((2*n0)+n1)/(2*n), # calc allele f
                           q = (1 - p),
                           MAF = pmin(p, q),
                           PIC = 1 - ((p^2) + (q^2)),
                           PId = (2*((p^2) + (q^2))^2)-((p^4) + (q^4))) %>%
        mutate(SNP_ID = SNP_ID) %>%
        select(SNP_ID, MAF, PIC, PId)
      
      snp_geno <- df %>%
        mutate(ind_per_snp = ncol(.) - rowSums(is.na(.)),
               geno_rate = (ind_per_snp / ncol(.))*100) %>%
        mutate(SNP_ID = SNP_ID)
      
      maf_geno <- left_join(maf, snp_geno, by = "SNP_ID")
      
      
      df <- reactive({
        
        df <- maf_geno %>%
          select(SNP_ID, MAF, geno_rate, PIC, PId) %>%
          # filter(geno_rate > 0 & MAF > 0.2)
          filter(geno_rate >= input$percent_geno & 
                   MAF >= input$maf_thresh)
        
      })
      
      
      
      genind <- reactive({ # needs to be the same as df
        
        SNP_ID <- df()$SNP_ID
        
        genind <- df() %>%
          select(-c(MAF, geno_rate, PIC, PId)) %>%
          left_join(maf_geno, by = "SNP_ID") %>%
          select(-c(SNP_ID, MAF, SNP_ID, ind_per_snp, geno_rate, PIC, PId)) %>%
          t(.)
        
        genind[genind == 0] <- "1/1" # homozygote reference
        genind[genind == 1] <- "1/2" # heterozygote
        genind[genind == 2] <- "2/2" # homozygote alternate
        
        colnames(genind) <- SNP_ID
        inds <- rownames(genind)
        genind <- df2genind(genind, ploidy=2, sep = "/", NA.char = "NA")
        #indNames(genind)
        #pop(genind) <- pop
        genind <- genind
        
      })
      
      
      genos <- reactive({
        
        genos <- maf_geno %>%
          #filter(geno_rate > 0 & MAF > 0.2) %>%
          filter(geno_rate >= input$percent_geno & 
                   MAF >= input$maf_thresh) %>%
          select(-c(SNP_ID, MAF, ind_per_snp, geno_rate, PIC, PId))
        
        
        genos[genos=="2"]<-0
        genos <- t(genos)
        
      })
      
      # output$MLH <- renderPlot({
      #   
      #   mlh <- data.frame(MLH = MLH(genos())) %>%
      #     mutate(ID = rownames(.))
      #   
      #   ggplot(mlh, aes(MLH)) +
      #     geom_histogram() +
      #     theme_emily()
      #   
      # })
      
      
      #~~ Summary table
      output$sum_stats <- function() {
        
        data.frame(x = c("Number of individuals", "Number of loci",
                         "Missing data (%)", "Expected heterozygosity", "Observed heterozygosity", "Probability of identity"),
                   b = c(ncol(maf_geno)-4, nrow(df()),
                         summary(genind())$NA.perc,
                         mean(summary(genind())$Hexp),
                         mean(summary(genind())$Hobs),
                         #mean(df()$PIC),
                         mean(df()$PId))) %>%
          # mutate(b = case_when(grepl("Number", b) ~ round(b),
          #                      !grepl("Number", b) ~ signif(b, 5))) %>%
          `colnames<-`(c("", "")) %>%
          knitr::kable("html") %>%
          kable_styling("striped", full_width = F, position = "left") %>%
          add_header_above(c("Summary Statistics" = 2), line = F) 
        
      }
      
      
      
      output$hwe <- DT::renderDataTable({
        hwe <- hw.test(genind(), B=0)
        
      })
      
      
      output$mytable = DT::renderDataTable({
        df()
      })
      
      
      output$geno <- renderPlot({
        ggplot(df(), aes(geno_rate)) +
          geom_histogram() +
          xlim(c(0,110)) +
          theme_emily() +
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)) +
          xlab("Genotyping rate") +
          ylab("Count")
        #   xlim(c(0,nind+1))
        
      })
      
      
      output$maf <- renderPlot({
        ggplot(df(), aes(MAF)) +
          geom_histogram() +
          xlim(c(0,0.5)) +
          theme_emily() +
          theme(axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)) +
          xlab("Minor Allele Frequency") +
          ylab("Count")
      })
      
      
    })
    
    
    
    
  }

