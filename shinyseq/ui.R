library("BiocManager")
options(repos = BiocManager::repositories())
getOption("repos")
library(shiny)
library(shinythemes)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(shinyjs)
library(sortable)
library(tibble)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)
library(cluster)

options(shiny.maxRequestSize = 10000*1024^2)



# Define UI
ui <- shinyUI(navbarPage(title = "shinySeq",
                         theme=shinytheme('flatly'),
                         
                         tabPanel(title = "Upload",
                                  sidebarLayout(
                                    sidebarPanel(
                                      fileInput("file", label = "Count Matrix (tab separated .txt file)"),
                                      fileInput("file2", label = "MetaData target file (tab separated .txt file)"),
                                      actionButton(inputId="run","RUN"),
                                      hr(),
                                      actionButton(inputId="example","RUN with Example Data Set "),
                                      hr(),
                                      hr(),
                                      fileInput("file3", label = "Previously Calculated DESeq2 object (.rds)"),
                                      actionButton(inputId = 'load', "Load data"),
                                      textOutput("log"),
                                      
                                      # plotOutput("hist"), 
                                      width = 3), 
                                    
                                    mainPanel(
                                      
                                      DT::dataTableOutput("mytable2"),
                                      hr(),
                                      DT::dataTableOutput("mytable1"), 
                                      width = 9
                                      
                                    )
                                  )
                         ),
                         
                         tabPanel(title = "Sample Clusters", plotOutput(outputId = "clust", width = 800, height = 600)),
                         
                         tabPanel("3D-PCA", 
                                  sidebarLayout(
                                    sidebarPanel(
                                      uiOutput("groups"),width = 3, height =3), 
                                    
                                    mainPanel(
                                      plotlyOutput("plot"))
                                  )
                         ),
                         
                         tabPanel("2D-PCA", 
                                  sidebarLayout(
                                    sidebarPanel(
                                      uiOutput("groups2"),width = 4,
                                      uiOutput("xaxis"),
                                      uiOutput("yaxis"),
                                      plotlyOutput(outputId = "scree"),
                                      hr(),
                                      downloadButton("downloadData", "Download Eigenvals"),
                                      hr(),
                                      verbatimTextOutput("prop")), 
                                    
                                    mainPanel(
                                      plotlyOutput("plot2"), width = 8)
                                  )
                         ),
                         
                         tabPanel(title = "DESeq2",
                                  
                                  sidebarLayout(fluid = F,
                                                sidebarPanel(
                                                  radioButtons('multifactor', 
                                                               'Select number of factors in precomputed DESeq',
                                                               c('Multifactored','Monofactored'),
                                                               selected='Multifactored'),
                                                  uiOutput("numerator"),
                                                  uiOutput("denominator"),
                                                  uiOutput('contrast_v'),
                                                  hr(),
                                                  plotOutput(outputId = "MAPlot"),
                                                  # downloadButton("normcounts", "Download Normalized Counts"),
                                                  downloadButton("results", "Download DE results"),
                                                  width = 4
                                                ),
                                                
                                                mainPanel(
                                                  #plotOutput("group4")
                                                  DT::dataTableOutput("contrast"),
                                                  hr()
                                                )
                                  )
                                  
                                  
                                  
                                  
                                  
                                  
                         ),
                         
                         tabPanel(title = "Gene Cluster",
                                  
                                  sidebarLayout(fluid = F,
                                                sidebarPanel(
                                                  
                                                  
                                                  uiOutput("within_group"),
                                                  uiOutput("within_i"),
                                                  uiOutput("between_group"),
                                                  uiOutput('gene_filter'),
                                                  uiOutput('topn'),
                                                  uiOutput('order'),
                                                  actionButton(inputId="run_heat","Generate Heatmap"),
                                                  hr(),
                                                  # downloadButton("normcounts", "Download Normalized Counts"),
                                                  #downloadButton("results", "Download Cluster Results"),
                                                  width = 4
                                                ),
                                                
                                                mainPanel(
                                                  plotOutput(outputId = 'hmap_plot')
                                                  #plotOutput("group4")
                                                  #DT::dataTableOutput("contrast"),
                                                  #hr(),
                                                  
                                                  #verbatimTextOutput("spit"), width = 7
                                                )
                                  )
                                  
                                  
                                  
                                  
                                  
                                  
                         ),
                         
                         shinyjs::useShinyjs(),
                         div(class = "footer",
                             includeHTML("footer.html")
                         )
                         
)
)

