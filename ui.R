# ui.R
shinyUI(fluidPage(
titlePanel("Phy5"),

sidebarLayout(
    sidebarPanel(
        helpText("Distance Matrix"),
        
        selectInput("distance", 
                    label = "Choose a distance matrix",
                    choices = c("euclidean","manhattan", "maximum",
                                "canberra", "binary"),
                    selected = "manhattan"),
        
        selectInput("method",
                    label = "Choose an agglomeration method",
                    choices = c("ward.D2", "average", "single",
                                "complete", "mcquitty", "median",
                                "centroid"),
                    selected = "ward.D2"),
        
        fileInput("fnafile", "Upload fasta files",
                  accept = c("text/fasta",".fna"),
                  multiple = TRUE)
    ),
    
    mainPanel(
        textOutput("text1"),
        
        textOutput("text2"),
        
        downloadButton('pentamertable','Download: raw frequencies'),
        
        downloadButton('proptable','Download: proportions'),
        
        uiOutput("ui_plot"),
        
        downloadButton("image", "Download: Phylogenetic Tree")
    )
)
))
