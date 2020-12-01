

#~~~~~~~~~~~~~~~~~~~~~#
#       Sidebar       #
#~~~~~~~~~~~~~~~~~~~~~#


shinyUI(pageWithSidebar(
  
  titlePanel(
    title=div(img(src="image.jpg", height = 150, width = 150), align = "center"), windowTitle = "wildlifedetectR"
    # div(img(src="image.jpg", height = 150, width = 150), align = "center", windowTitle = "wildlifedetectR")
    
    #title = h1("wildlife_detectR", align = "center"))
    #titlePanel((title=div("wildlife detectR",(img(src="heading.jpg", height = 150, width = 150)), align = "center"))
  ),
  
  sidebarPanel(
    
    ## conditionalPanel() functions for selected tab
    
    #~~ About
    
    conditionalPanel(condition="input.tabselected==1", 
                     h2("wildlife detectR"), 
                     h4("This app was built using R and shiny"),
                     h5(em("Please refer to the manual for instructions and analysis methods."))),
    
    #~~ Species Identification
    
    conditionalPanel(condition="input.tabselected==2", h3("Species Identification"),
                     fileInput("SpID_unknown", "Unknown sequence (fasta)"),
                     fileInput("SpID_ref", "Reference sequences (fasta)"),
                     radioButtons("SpID_align","Align", choices=c("Yes" = T, "No" = F)),
                     radioButtons("method1","Choose a method", choices=c("Phylogeny" = 1, "Sequence similarity" = 2))
    ),
    
    conditionalPanel(
      condition = "input.tabselected==2 & input.method1 == 1",
      radioButtons("SpID_PhyloOut","Outgroup", choices=c("Yes" = T, "No" = F)),
      selectInput("SpID_PhyloBoot", "Bootstraps",
                  choices = c(10,100,1000),
                  selected = 10),
      actionButton("SpID_action1", "Submit!"),
      downloadButton("SpID_report_phylo", "Generate report")
    ),
    
    conditionalPanel(
      condition = "input.tabselected==2 & input.method1 == 2",
      selectInput("SpID_SeqSimModel", "Substitution model",
                  choices = c("raw", "N", "TS", "TV", "JC69", "K80", 
                              "F81", "K81", "F84", "BH87", "T92", 
                              "TN93", "GG95", "logdet", "paralin", 
                              "indel", "indelblock"),
                  selected = "K80"),
      selectInput("SpID_SeqSimHits", "Number of hits to report",
                  choices = c("all", "1", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10")),
      actionButton("SpID_action2", "Submit!"),
      downloadButton("SpID_report_seqsim", "Generate report")
    ),
    
    
    
    #~~ Individualization
    
    conditionalPanel(
      condition = "input.tabselected==3",
      h3("Individualization"),
      fileInput("UnknownProfile", "Trace profile"),
      fileInput("Ref_db", "Reference population data"),
      radioButtons(
        "method2",
        "Reference population data format",
        choices = c(
          "Allele frequency table (Genalex format)" = T,
          "Genotypes" = F
        )
      ),
      radioButtons(
        "method10",
        "Choose a method",
        choices = c("Random match probability", "Likelihood ratio approach")
      ),
      radioButtons("method11", "Known alleles", choices =
                     c("Yes" = T, "No" = F)),
      numericInput(
        "IndividualizationTheta",
        "Theta",
        value = 0,
        min = 0,
        max = 1,
        step = 0.01
      ),
      #selectInput("ProfileInput", "File type", choices=c("Genetix", "Genalex", "Genepop", "DArT"), selected = "Genepop"),
      actionButton("action2", "Submit!"),
      downloadButton("ind_report", "Generate report")
    ), 
    
    #~~ Parentage Testing
    
    conditionalPanel(condition="input.tabselected==6", h3("Parentage"),
                     selectInput("ParentageTest", "Parentage test type", choices=c("Probability of joint parentage", "Probability of paternity (single parentage)", "Probability of motherless paternity"), selected = "Test 1"),
                     fileInput("ParentageProfile", "Parent offspring file"),
                     fileInput("ParentageDB", "Reference population data"),
                     selectInput("ParentageMarker", "Marker type", choices=c("Microsatellites", "SNPs"), selected = "Microsatellites"),
                     radioButtons("ParentageData","Data type", choices=c("Allele frequency table (Genalex format)", "Genotypes")),
                     radioButtons("ParentageAlleles","Known alleles", choices=c("Yes" = T, "No" = F)),
                     numericInput("ParentageTheta", "Theta", value = 0, min = 0, max = 1, step = 0.01),
                     #textInput("ParentageMutRate", "Mutation rate"),
                     actionButton("action5", "Submit!")
    ),
    
    #~~ Haplotype Network
    
    conditionalPanel(condition="input.tabselected==5", h3("Haplotype Network"), 
                     #fileInput("UnknownPopSeqs", "Unknown sequence file (fasta)"),
                     fileInput("ReferencePopSeqs", "Sequences file (fasta)"),
                     fileInput("PopList", "List of populations"),
                     sliderInput("HapScaleRatio", "Plotting scale ratio", min = 0, max = 3, value = 0.3, step = 0.1),
                     sliderInput("HapX", "Legend X coordinate", min = -50, max = 0, value = -8, step = 1),
                     sliderInput("HapY", "Legend Y coordinate", min = 0, max = 50, value = 14, step = 1),
                     actionButton("action4", "Submit!")
    ),
    

    #~~ Marker Selection
    conditionalPanel(condition="input.tabselected==7", h3("Marker Selection"), 
                     fileInput("SNPdata", "SNP dataset"),
                     sliderInput("percent_geno", "Select genotyping rate", min = 0, max = 100, value = 0, step = 1),
                     sliderInput("maf_thresh", "Minor allele frequency threshold", min = 0, max = 0.5, value = 0, step = 0.01),
                     # radioButtons("hwe_filter","Hardy-Weinberg Filter", choices=c("Yes" = T, "No" = F)),
                     actionButton("action6", "Submit!"))
  ),
  
  
  #~~~~~~~~~~~~~~~~~~~~~#
  #       Output        #
  #~~~~~~~~~~~~~~~~~~~~~#
  
  mainPanel(
    # recommend review the syntax for tabsetPanel() & tabPanel() for better understanding
    # id argument is important in the tabsetPanel()
    # value argument is important in the tabPanel()
    tabsetPanel(
      tabPanel("About", value=1,
               
               p(strong("wildlife detectR"), "is an interactive user-friendly tool for undertaking key wildlife 
                 forensic and population genetic analyses.", style = "font-size:20px"),
               p("The purpose of wildlife detectR is to streamline and standardize the analytical methods typically
                 used in DNA based wildlife forensic science. The program can generate statistics and reports for individualisation, 
                 parentage and species identification tests. Both traditional microsatellite and genomic-based single nucleotide 
                 polymorphism (SNP) markers can be input. It can also be used to generate population genetic statistics and 
                 undertake marker filtering to aid in the development of suitable datasets for wildlife forensic analysis.", style = "font-size:16px"),
               helpText(em("Equations and algorithms used for each test can be referred to in the manual which will be made available in the near future.")),
               helpText(em("Please send any queries to Emily Humble or Kyle Ewart."))
               # p("p creates a paragraph of text."),
               # p("A new p() command starts a new paragraph. Supply a style attribute to change the format of the entire paragraph.", style = "font-family: 'times'; font-si16pt"),
               # strong("strong() makes bold text."),
               # em("em() creates italicized (i.e, emphasized) text."),
               # br(),
               # code("code displays your text similar to computer code"),
               # div("div creates segments of text with a similar style. This division of text is all blue because I passed the argument 'style = color:blue' to div", style = "color:blue"),
               # br(),
               # p("span does the same thing as div, but it works with",
               #   span("groups of words", style = "color:blue"),
               #   "that appear inside a paragraph."),
               # helpText("Spiel about package, including links to the manual and associated paper.")
      ),
      
      tabPanel("Species Identification", value=2,
               helpText(strong("Unknown sequence:"),("Fasta file containing sequence for unknown individual")),
               helpText(strong("Reference sequences:"), ("Fasta file containing reference sequences")),
               helpText(strong("Align:"), ("Do you need to align the unknown sequence to the reference sequences?")),
               helpText(strong("Choose a method:"), ("Do you want to perform a phylogenetic or sequence similarity comparison?")),
               helpText(strong("Outgroup (phylogeny):"), ("Do your reference sequences contain an outgroup?")),
               helpText(strong("Bootstraps (phylogeny):"), ("Number of bootstrap replicates")),
               helpText(strong("Substitution model (sequence similarity):"), ("Evolutionary model for distance matrix")),
               helpText(strong("Number of hits to report (sequence similarity):"), ("Number of hits to the reference database to report (useful if your reference dataset is very large)")),
               hr(),
               plotOutput("SpIDAlign", width = "100%"),
               conditionalPanel(condition = "input.SpID_action1", plotOutput("SpIDPhylo", width = "100%")),
               conditionalPanel(condition = "input.SpID_action2", tableOutput("SpIDSeqSim"))),
      
      tabPanel("Individualization", value=3,
               helpText(strong("Trace profile:"),("File containing genotypes for case individual")),
               helpText(strong("Reference population data:"), ("File containing genotypes or allele frequencies for reference population")),
               helpText(strong("Reference population data format:"), ("Does your reference population data contain allele frequencies or genotypes?")),
               helpText(strong("Choose a method:"), ("Do you want to report random match probability or a likelihood ratio?")),
               helpText(strong("Known alleles:"), ("Are your alleles for your unknown equivalent to those in the potential matches file?")),
               helpText(strong("Theta:"), ("Value for theta")),
               hr(),
               htmlOutput("Individualisation")),
      
      tabPanel("Parentage", value = 6,
               helpText(strong("Parentage test type:"), ("Parentage test to perform")),
               helpText(strong("Parent offspring file:"), ("File containing genotypes for offspring and putative parent(s)")),
               helpText(strong("Reference population data:"), ("File containing genotypes or allele frequencies for reference population")),
               helpText(strong("Reference population data format:"), ("Does your reference population data contain allele frequencies or genotypes?")),
               helpText(strong("Known alleles:"), ("Does your reference population data contain the same alleles as your parent offspring file")),
               helpText(strong("Theta:"), ("Value for theta")),
               hr(),
               textOutput("Parentage_result"),
               htmlOutput("Parentage_statement")),
      
      tabPanel("Haplotype Network", value =5,
               helpText(strong("Sequences file:"),("Fasta file containing sequences of individuals from known locations and the sequence of individual from unknown location.")),
               helpText(strong("List of populations:"),("If population IDs are not contained within sequence headers of the database sequences fasta, please upload a csv file containing the population IDs.")),
               hr(),
               plotOutput("hapNet")),
      
      tabPanel("Marker Selection", value=7,
               helpText(strong("SNP dataset:"),("Upload a SNP dataset in PLINK raw format")),
               hr(),
               tableOutput("sum_stats"),
               #plotOutput("maf", width = "40%"),
               #plotOutput("geno", width = "40%"),
               splitLayout(
                 cellWidths = "50%",
                 cellArgs = list(style = "padding: 6px"),
                 plotOutput("maf"),
                 plotOutput("geno")
               )),
      id = "tabselected"
    )
  )
))