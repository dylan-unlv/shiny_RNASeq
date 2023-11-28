# Define server logic
server <- shinyServer(function(input, output, session) {
  
  
  observeEvent(input$run,{
    
    
    req(input$file)
    
    
    counts <- read.table(input$file$datapath,
                         header = T,
                         sep = "\t",
                         row.names = 1)
    counts = counts[ , order(names(counts))]
    
    req(input$file2)
    
    
    renderUI(actionButton("rmv", "x"),)
    
    target <- read.table(input$file2$datapath,
                         header = T,
                         sep = "\t")
    
    rownames(target) <- target$label
    
    target = target[order(rownames(target)), ]
    
    
    
    # load the file into new environment and get it from there
    e = new.env()
    
    progress <- Progress$new(session, min=1, max=20)
    on.exit(progress$close())
    
    progress$set(message = 'Execution in progress',
                 detail = 'This may take a while...')
    for (i in 1:4) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    
    # Plot the data
    #target2 <- target %>% select(label, group)
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(target, list(pageLength = 5, scrollX=T),rownames = F)
    })
    output$mytable2 <- DT::renderDataTable({
      DT::datatable(counts, list(pageLength = 5, scrollX=T))
    })
    
    output$prop <- renderPrint(summary(pca))
    
    
    
    # ################################
    # ################################
    # # DESEQ2
    #
    library("DESeq2")
    library("dplyr")
    library("tidyverse")
    library("plotly")
    #
    
    ## -------------------------------------------------------------------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = target,
                                  design = ~ group)
    
    #
    #
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    for (i in 4:12) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    
    withCallingHandlers({
      shinyjs::html("log", "DESeq2 Log: ")
      dds <- DESeq(dds)
    },
    message = function(m) {
      shinyjs::html(id = "log", html = m$message, add = TRUE)
      shinyjs::html(id = "log", html = "<br>", add = TRUE)
    }, collapse="<br>")
    
    
    
    #resultsNames(dds)
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    vsd <- varianceStabilizingTransformation(dds, blind=T)
    
    for (i in 13:17) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    # calculate the variance for each gene
    rv <- rowVars(assay(vsd))
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(vsd)[select,]))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    pVar.df <- as.data.frame(percentVar)
    pVar.df$x = as.factor(paste0("PC",rownames(pVar.df)))
    
    pVar.df = pVar.df[ , order(names(pVar.df))]
    pVar.df$percentVar = pVar.df$percentVar * 100
    pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)
    
    group = target$group
    intgroup.df <- as.data.frame(colData(vsd)[, "group", drop=FALSE])
    
    # assembly the data for the plot
    
    
    d <- data.frame(pca$x, name=rownames(pca$x))
    d2 <- left_join(target, d, by= c("label"="name"))
    
    ### -------------------------------------------------------------------------------------------------------------------
    
    hc2 <- hclust(dist(t(assay(vsd))), method="ward.D")
    
    output$clust <- renderPlot(plot(hc2, hang=-1, ylab="Height", las=2,
                                    xlab="Method: Euclidean distance - Ward criterion",
                                    main="Cluster Dendrogram"))
    
    
    output$groups <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice", label = "Select Groupings",
                    choices = colnames(d2)[colnames(d2) %in% colnames(target)], selected = "group")
      )
      
      
      
    })
    
    output$groups2 <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice2", label = "Select Groupings",
                    choices = colnames(d2)[colnames(d2) %in% colnames(target)], selected = "group")
      )
      
      
      
    })
    
    
    output$xaxis <- renderUI({
      
      
      tagList(
        selectInput(inputId = "pcx", label = "Select x-axis PC",
                    choices = colnames(d2)[grep(pattern = "PC", colnames(d2))], selected = "PC1")
      )
      
      
      
    })
    
    
    output$yaxis <- renderUI({
      
      
      tagList(
        selectInput(inputId = "pcy", label = "Select y-axis PC",
                    choices = colnames(d2)[grep(pattern = "PC", colnames(d2))], selected = "PC2")
      )
      
      
      
    })
    
    
    
    m <- list(
      l = 50,
      r = 50,
      b = 100,
      t = 100,
      pad = 4
    )
    output$plot <- renderPlotly({
      
      if(is.null(input$choice)){return()}
      plot_ly(data = d2, x = ~ PC1, y= ~ PC2, z = ~ PC3, 
              width = 800, height = 800,
              color = ~ get(input$choice),
              #colors = color2,
              marker = list(size = 8,
                            line = list(color = ~ label , width = 1))) %>%
        add_markers() %>%
        layout(autosize = F, margin =m, title = "First 3 Principal Components",
               scene = list(xaxis = list(title = 'PC1'),
                            yaxis = list(title = 'PC2'),
                            zaxis = list(title = 'PC3')
               ), margin = m
        ) %>% layout(legend = list(orientation = 'h', xanchor = "center",
                                   x = 0.5, y = -0.1))
    })
    
    
    
    x = reactive({input$pcx})
    y = reactive({input$pcy})
    
    
    plotdata <- reactive({
      
      dplyr::select(d2, !contains("PC"))
      
    })
    
    
    
    
    output$plot2 <- renderPlotly({
      
      if(is.null(input$choice2)){return()}
      
      d3   <- plotdata()
      d3$x <- d2[[input$pcx]]
      d3$y <- d2[[input$pcy]]
      
      plot_ly(data = d3 , x = ~ x, y =  ~ y, 
              width = 850, height = 800,
              color = ~ get(input$choice2),
              marker = list(size = 12,
                            line = list(color = ~ label, width = 1))) %>%
        add_markers() %>%
        layout(autosize = F, title = paste0(input$pcy, " vs ", input$pcx),
               xaxis = list(title = input$pcx),
               yaxis = list(title = input$pcy), margin = m) %>% 
        layout(legend = list(orientation = 'h', xanchor = "center",
                             x = 0.5, y = -0.2))
      
    })
    
    output$scree <- renderPlotly(
      
      plot_ly(data = pVar.df, x = ~x, y = ~ percentVar, type = "scatter", mode = "markers",
              width = 400, height = 400) %>%
        layout(autosize = F, margin = m,
               xaxis = list(categoryorder = "array",title = "nth PC", categoryarray = ~x),
               yaxis = list(title = "Percent Var", ticksuffix = "%"),
               title = "elbow plot")
    )
    
    
    
    
    for (i in 18:20) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    
    
    # save eigenvals to a new object 
    datasetInput <- pca$x
    
    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("eigenvalues.csv", sep = "")
      },
      content = function(file) {
        write.csv(datasetInput, file, row.names = T)
      }
    )
    
    
    
    
    #----------------------------------------------------------
    #----------------------------------------------------------
    #----------------------------------------------------------
    #-- Differential Expression Analysis
    
    # Make Normalized Counts Downloadable
    
    normalized_counts <- counts(dds, normalized=TRUE)
    colnames(normalized_counts) = paste0("norm.", colnames(normalized_counts))
    normalized_counts <- round(normalized_counts, digits = 0)
    # Downloadable csv of normalized counts dataset ----
    
    output$normcounts <- downloadHandler(
      filename = function() {
        paste("normalizedCounts.csv", sep = "")
      },
      content = function(file) {
        write.csv(normalized_counts, file, row.names = T, quote = F, col.names = NA)
      }
    )
    
    output$numerator <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice4", label = "Select Group of interest for DE",
                    choices = dds$group)
      )
    })
    
    output$denominator <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice3", label = "Select Base-level for DE",
                    choices = dds$group)
      )
    })
    
    
    
    dds.results <- reactive({
      
      req(input$choice3)
      req(input$choice4)
      
      results(dds, contrast=c("group", as.character(input$choice4), 
                              as.character(input$choice3)), alpha = 0.05)
      
    })        
    
    de.data <- reactive({
      
      req(input$choice3)
      req(input$choice4)
      round(as.data.frame(dds.results()),3)
      
      
    })
    
    nc <- as.data.frame(normalized_counts)
    nc$features <- rownames(nc)
    
    de.merge <- reactive({
      
      dds.nc <- de.data()
      dds.nc$features <- rownames(dds.nc)
      nc2 <- left_join(nc, dds.nc, by = "features")
      nc2 <- nc2 %>% 
        select("features", everything())
    })
    
    
    output$contrast <- DT::renderDataTable({
      
      DT::datatable(de.merge(), list(pageLength = 10, scrollX=T), rownames = F)
      
    })
    
    ma.data <- reactive({
      
      req(input$choice3)
      req(input$choice4)
      plotMA(dds.results(), main = paste0("MAPlot ", input$choice4, "_vs_", input$choice3)) 
      
    })
    
    output$MAPlot <-
      renderPlot(ma.data(), width = 400, height = 400)
    
    
    
    output$results <- downloadHandler(
      filename = function() {
        paste0(input$choice4, "_vs_", input$choice3, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(as.data.frame(de.merge()), file, row.names = T, quote = F, col.names = NA)
      }
    )
    
    
    
    
    output$spit <-
      
      renderPrint(summary(dds.results()))
    
  })
  
  
  
  
  #------ Example Data Set ------------
  
  observeEvent(input$example,{
    
    
    #req(input$file)
    
    counts <- read.table("example/countMatrix.txt",
                         header = T,
                         sep = "\t",
                         row.names = 1)
    counts = counts[ , order(names(counts))]
    
    #req(input$file2)
    
    renderUI(actionButton("rmv", "x"),)
    
    target <- read.table("example/targetFile.txt",
                         header = T,
                         sep = "\t")
    
    rownames(target) <- target$label
    
    target = target[order(rownames(target)), ]
    
    
    
    # load the file into new environment and get it from there
    e = new.env()
    
    progress <- Progress$new(session, min=1, max=20)
    on.exit(progress$close())
    
    progress$set(message = 'Execution in progress',
                 detail = 'This may take a while...')
    for (i in 1:4) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    
    # Plot the data
    #target2 <- target %>% select(label, group)
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(target, list(pageLength = 5, scrollX=T), rownames = F)
    })
    output$mytable2 <- DT::renderDataTable({
      DT::datatable(counts, list(pageLength = 5, scrollX=T))
    })
    
    output$prop <- renderPrint(summary(pca))
    
    
    
    # ################################
    # ################################
    # # DESEQ2
    #
    library("DESeq2")
    library("dplyr")
    library("tidyverse")
    library("plotly")
    #
    
    ## -------------------------------------------------------------------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = target,
                                  design = ~ group)
    
    #
    #
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    for (i in 4:12) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    
    withCallingHandlers({
      shinyjs::html("log", "DESeq2 Log: ")
      dds <- DESeq(dds)
    },
    message = function(m) {
      shinyjs::html(id = "log", html = m$message, add = TRUE)
      shinyjs::html(id = "log", html = "<br>", add = TRUE)
    }, collapse="<br>")
    
    
    
    #resultsNames(dds)
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    vsd <- varianceStabilizingTransformation(dds, blind=T)
    
    for (i in 13:17) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    #
    # ## -------------------------------------------------------------------------------------------------------------------
    # calculate the variance for each gene
    rv <- rowVars(assay(vsd))
    
    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(vsd)[select,]))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    pVar.df <- as.data.frame(percentVar)
    pVar.df$x = as.factor(paste0("PC",rownames(pVar.df)))
    
    pVar.df = pVar.df[ , order(names(pVar.df))]
    pVar.df$percentVar = pVar.df$percentVar * 100
    pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)
    group = target$group
    intgroup.df <- as.data.frame(colData(vsd)[, "group", drop=FALSE])
    
    # assembly the data for the plot
    
    #d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], name=colnames(vsd))
    
    d <- data.frame(pca$x, name=rownames(pca$x))
    
    d2 <- left_join(target, d, by= c("label"="name"))
    
    ### -------------------------------------------------------------------------------------------------------------------
    
    hc2 <- hclust(dist(t(assay(vsd))), method="ward.D")
    
    output$clust <- renderPlot(plot(hc2, hang=-1, ylab="Height", las=2,
                                    xlab="Method: Euclidean distance - Ward criterion",
                                    main="Cluster Dendrogram"))
    
    output$groups <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice", label = "Select Groupings",
                    choices = colnames(d2)[colnames(d2) %in% colnames(target)], selected = "group")
      )
      
      
      
    })
    
    output$groups2 <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice2", label = "Select Groupings",
                    choices = colnames(d2)[colnames(d2) %in% colnames(target)], selected = "group")
      )
      
      
      
    })
    
    
    output$xaxis <- renderUI({
      
      
      tagList(
        selectInput(inputId = "pcx", label = "Select x-axis PC",
                    choices = colnames(d2)[grep(pattern = "PC", colnames(d2))], selected = "PC1")
      )
      
      
      
    })
    
    
    output$yaxis <- renderUI({
      
      
      tagList(
        selectInput(inputId = "pcy", label = "Select y-axis PC",
                    choices = colnames(d2)[grep(pattern = "PC", colnames(d2))], selected = "PC2")
      )
      
      
      
    })
    
    
    m <- list(
      l = 50,
      r = 50,
      b = 100,
      t = 100,
      pad = 4
    )
    
    output$plot <- renderPlotly({
      
      if(is.null(input$choice)){return()}
      plot_ly(data = d2, x = ~ PC1, y= ~ PC2, z = ~ PC3, 
              width = 800, height = 800,
              color = ~ get(input$choice),
              #colors = color2,
              marker = list(size = 8,
                            line = list(color = ~ label , width = 1))) %>%
        add_markers() %>%
        layout(autosize = F, margin =m, title = "First 3 PC Dimensions",
               scene = list(xaxis = list(title = 'PC1'),
                            yaxis = list(title = 'PC2'),
                            zaxis = list(title = 'PC3')
               )
        ) %>% layout(legend = list(orientation = 'h', xanchor = "center",
                                   x = 0.5, y = -0.2))
    })
    
    
    x = reactive({input$pcx})
    y = reactive({input$pcy})
    
    
    plotdata <- reactive({
      
      dplyr::select(d2, !contains("PC"))
      
    })
    
    
    
    output$plot2 <- renderPlotly({
      
      if(is.null(input$choice2)){return()}
      
      d3   <- plotdata()
      d3$x <- d2[[input$pcx]]
      d3$y <- d2[[input$pcy]]
      
      plot_ly(data = d3 , x = ~ x, y =  ~ y, 
              width = 850, height = 800,
              color = ~ get(input$choice2),
              marker = list(size = 10,
                            line = list(color = ~ label, width = 1))) %>%
        add_markers() %>%
        layout(autosize = F, margin = m2, title = paste0(input$pcy, " vs ", input$pcx),
               xaxis = list(title = input$pcx),
               yaxis = list(title = input$pcy)) %>% layout(legend = list(orientation = 'h', xanchor = "center",
                                                                         x = 0.5, y = -0.2))
      
    })
    
    
    m2<- list(
      l = 70,
      r = 30,
      b = 100,
      t = 100,
      pad = 2
    )
    
    
    output$scree <- renderPlotly(
      
      plot_ly(data = pVar.df, x = ~x, y = ~ percentVar, type = "scatter", mode = "markers",
              width = 400, height = 400) %>%
        layout(autosize = F, margin = m2,
               xaxis = list(categoryorder = "array",title = "nth PC", categoryarray = ~x),
               yaxis = list(title = "Percent Var", ticksuffix = "%"),
               title = "elbow plot")
    )
    
    
    for (i in 18:20) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    
    
    # save eigenvals to a new object 
    datasetInput <- pca$x
    
    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("eigenvalues.csv", sep = "")
      },
      content = function(file) {
        write.csv(datasetInput, file, row.names = T)
      }
    )
    
    
    
    
    #----------------------------------------------------------
    #----------------------------------------------------------
    #----------------------------------------------------------
    #-- Differential Expression Analysis
    
    # Make Normalized Counts Downloadable
    
    normalized_counts <- counts(dds, normalized=TRUE)
    colnames(normalized_counts) = paste0("norm.", colnames(normalized_counts))
    normalized_counts <- round(normalized_counts, digits = 0)
    
    # Downloadable csv of normalized counts dataset ----
    
    # output$normcounts <- downloadHandler(
    #     filename = function() {
    #         paste("normalizedCounts.csv", sep = "")
    #     },
    #     content = function(file) {
    #         write.csv(normalized_counts, file, row.names = T, quote = F, col.names = NA)
    #     }
    # )
    
    output$numerator <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice4", label = "Select Group of interest for DE",
                    choices = dds$group)
      )
    })
    
    output$denominator <- renderUI({
      
      
      tagList(
        selectInput(inputId = "choice3", label = "Select Base-level for DE",
                    choices = dds$group)
      )
    })
    
    
    
    dds.results <- reactive({
      
      req(input$choice3)
      req(input$choice4)
      
      results(dds, contrast=c("group", as.character(input$choice4), 
                              as.character(input$choice3)), alpha = 0.05)
      
    })        
    
    de.data <- reactive({
      
      req(input$choice3)
      req(input$choice4)
      round(as.data.frame(dds.results()),3) 
      
      
      
    })
    
    nc <- as.data.frame(normalized_counts)
    nc$features <- rownames(nc)
    
    de.merge <- reactive({
      
      dds.nc <- de.data()
      dds.nc$features <- rownames(dds.nc)
      nc2 <- left_join(nc, dds.nc, by = "features")
      nc2 <- nc2 %>% 
        select("features", everything())
    })
    
    
    output$contrast <- DT::renderDataTable({
      
      DT::datatable(de.merge(), list(pageLength = 10, scrollX=T), rownames = F)
      
    })
    
    
    ma.data <- reactive({
      
      req(input$choice3)
      req(input$choice4)
      plotMA(dds.results(), main = paste0("MAPlot ", input$choice4, "_vs_", input$choice3)) 
      
    })
    
    output$MAPlot <-
      renderPlot(ma.data(), width = 400, height = 400)
    
    
    
    output$results <- downloadHandler(
      filename = function() {
        paste0(input$choice4, "_vs_", input$choice3, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(as.data.frame(de.merge()), file, row.names = T, quote = F, col.names = NA)
      }
    )
    
    
    
    
    output$spit <-
      
      renderPrint(summary(dds.results()))
    
    
    
    
    
  })
  
  
  
  
  
  
  #---- Load existing DESeq data ----
  
  observeEvent(input$load, {
    
    
    ###import deseq object
    req(input$file3)
    
    withCallingHandlers({
      shinyjs::html("log", "DESeq2 Log: ")
      dds <- readRDS(input$file3$datapath)
    },
    message = function(m) {
      shinyjs::html(id = "log", html = m$message, add = TRUE)
      shinyjs::html(id = "log", html = "<br>", add = TRUE)
    }, collapse="<br>")

    #output counts
    dds_count <- data.frame(counts(dds)) %>% rownames_to_column(var='gene') %>% relocate(gene)
    output$mytable2 <- DT::renderDataTable({
      DT::datatable(dds_count, list(pageLength = 10, scrollX=T),rownames = F)
    })
    
    #output metadata
    dds_meta <- data.frame(colData(dds))
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(dds_meta, list(pageLength = 5, scrollX=T),rownames = F)
    })
    target <- data.frame(colData(dds)) 
    
    #populate sample cluster
    hc2 <- hclust(dist(t(assay(dds))), method="ward.D")
    output$clust <- renderPlot(plot(hc2, hang=-1, ylab="Height", las=2,
                                    xlab="Method: Euclidean distance - Ward criterion",
                                    main="Cluster Dendrogram"))
    
    #populate PCA
    vsd <- varianceStabilizingTransformation(dds, blind=T)
    
      # calculate the variance for each gene
      rv <- rowVars(assay(vsd))
    
      # select the ntop genes by variance
      select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
    
      # perform a PCA on the data in assay(x) for the selected genes
      pca <- prcomp(t(assay(vsd)[select,]))
    
       # the contribution to the total variance for each component
      percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
      pVar.df <- as.data.frame(percentVar)
      pVar.df$x = as.factor(paste0("PC",rownames(pVar.df)))
      
      pVar.df = pVar.df[ , order(names(pVar.df))]
      pVar.df$percentVar = pVar.df$percentVar * 100
      pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)
      group = target$sample
      intgroup.df <- as.data.frame(colData(vsd)[, "sample", drop=FALSE])
    
      # assemble the data for the plot
      d <- data.frame(pca$x, name=rownames(pca$x))
    
      d2 <- left_join(target, d, by= c("sample"="name"))
      
      #contrasts (1)
      output$groups <- renderUI({
        
        
        tagList(
          selectInput(inputId = "choice", label = "Select Groupings",
                      choices = colnames(d2)[colnames(d2) %in% colnames(target)], selected = "sample")
        )
        
        
        
      })
      
      #contrasts (2)
      output$groups2 <- renderUI({
        
        
        tagList(
          selectInput(inputId = "choice2", label = "Select Groupings",
                      choices = colnames(d2)[colnames(d2) %in% colnames(target)], selected = "sample")
        )
        
        
        
      })
      
      
      output$xaxis <- renderUI({
        
        
        tagList(
          selectInput(inputId = "pcx", label = "Select x-axis PC",
                      choices = colnames(d2)[grep(pattern = "PC", colnames(d2))], selected = "PC1")
        )
        
        
        
      })
      
      
      output$yaxis <- renderUI({
        
        
        tagList(
          selectInput(inputId = "pcy", label = "Select y-axis PC",
                      choices = colnames(d2)[grep(pattern = "PC", colnames(d2))], selected = "PC2")
        )
        
        
        
      })
      

    
      
      
      #create group variable, should include all levels of DESeq design matrix
      

      
      m <- list(
        l = 50,
        r = 50,
        b = 100,
        t = 100,
        pad = 4
      )
      output$plot <- renderPlotly({
        
        if(is.null(input$choice)){return()}
        plot_ly(data = d2, x = ~ PC1, y= ~ PC2, z = ~ PC3, 
                width = 800, height = 800,
                color = ~ get(input$choice),
                #colors = color2,
                marker = list(size = 8,
                              line = list(color = ~ sample , width = 1))) %>%
          add_markers() %>%
          layout(autosize = F, margin =m, title = "First 3 Principal Components",
                 scene = list(xaxis = list(title = 'PC1'),
                              yaxis = list(title = 'PC2'),
                              zaxis = list(title = 'PC3')
                 ), margin = m
          ) %>% layout(legend = list(orientation = 'h', xanchor = "center",
                                     x = 0.5, y = -0.1))
      })
      
      
      
      x = reactive({input$pcx})
      y = reactive({input$pcy})
      
      
      plotdata <- reactive({
        
        dplyr::select(d2, !contains("PC"))
        
      })
      
      
      
      
      output$plot2 <- renderPlotly({
        
        if(is.null(input$choice2)){return()}
        
        d3   <- plotdata()
        d3$x <- d2[[input$pcx]]
        d3$y <- d2[[input$pcy]]
        
        plot_ly(data = d3 , x = ~ x, y =  ~ y, 
                width = 850, height = 800,
                color = ~ get(input$choice2),
                marker = list(size = 12,
                              line = list(color = ~ sample, width = 1))) %>%
          add_markers() %>%
          layout(autosize = F, title = paste0(input$pcy, " vs ", input$pcx),
                 xaxis = list(title = input$pcx),
                 yaxis = list(title = input$pcy), margin = m) %>% 
          layout(legend = list(orientation = 'h', xanchor = "center",
                               x = 0.5, y = -0.2))
        
      })
      
      output$scree <- renderPlotly(
        
        plot_ly(data = pVar.df, x = ~x, y = ~ percentVar, type = "scatter", mode = "markers",
                width = 400, height = 400) %>%
          layout(autosize = F, margin = m,
                 xaxis = list(categoryorder = "array",title = "nth PC", categoryarray = ~x),
                 yaxis = list(title = "Percent Var", ticksuffix = "%"),
                 title = "elbow plot")
      )
      
      
      ###DESeq 
      observeEvent(input$multifactor, {
        if (input$multifactor=='Multifactored'){
          if ('group' %in% names(colData(dds))){
            print('multi')
          } else { print('ERR: group variable not found in deseq object')}
        } else {
          
          #create group variable, should include all levels of DESeq design matrix
          cnames <- colnames(data.frame(colData(dds)) %>% select(-sample))
          gdat <- data.frame(colData(dds)) %>% unite(col='group', all_of(cnames), sep='_', remove=F) %>% pull(group)
          dds$group <- as.factor(gdat)
        }
        
        
      })
      
      
      
      
      
      normalized_counts <- counts(dds, normalized=TRUE)
      colnames(normalized_counts) = paste0("norm.", colnames(normalized_counts))
      normalized_counts <- round(normalized_counts, digits = 0)
      output$numerator <- renderUI({
        
        
        tagList(
          selectInput(inputId = "choice4", label = "Select Group of interest for DE",
                      choices = dds$group)
        )
      })
      
      output$denominator <- renderUI({
        
        
        tagList(
          selectInput(inputId = "choice3", label = "Select Base-level for DE",
                      choices = dds$group)
        )
      })
      
      output$contrast_v <- renderUI({
        
        
        tagList(
          selectInput(inputId = "choicecontrast", label = "Select Variable to Contrast",
                      choices = colnames(data.frame(colData(dds)) %>% select(-sample)))
        )
      })
      
      
      dds.results <- reactive({
        
        req(input$choice3)
        req(input$choice4)
        
        results(dds, contrast=c(as.character(input$choicecontrast), as.character(input$choice4), 
                                as.character(input$choice3)), alpha = 0.05)
        
      })        
      
      de.data <- reactive({
        
        req(input$choice3)
        req(input$choice4)
        round(as.data.frame(dds.results()),3)
        
        
      })
      
      nc <- as.data.frame(normalized_counts)
      nc$features <- rownames(nc)
      
      de.merge <- reactive({
        
        dds.nc <- de.data()
        dds.nc$features <- rownames(dds.nc)
        nc2 <- left_join(nc, dds.nc, by = "features")
        nc2 <- nc2 %>% 
          select("features",'baseMean', 
                 'log2FoldChange', 'lfcSE',
                 'stat', 'pvalue',
                 'padj', everything())
      })
      
      
      output$contrast <- DT::renderDataTable({
        
        DT::datatable(de.merge(), list(pageLength = 10, scrollX=T), rownames = F)
        
      })
      
      ma.data <- reactive({
        
        req(input$choice3)
        req(input$choice4)
        plotMA(dds.results(), main = paste0("MAPlot ", input$choice4, "_vs_", input$choice3)) 
        
      })
      
      output$MAPlot <-
        renderPlot(ma.data(), width = 400, height = 400)
      
      
      
      output$results <- downloadHandler(
        filename = function() {
          paste0(input$choice4, "_vs_", input$choice3, ".csv", sep = "")
          
        },
        content = function(file) {
          write.csv(as.data.frame(de.merge()), file, row.names = T, quote = F, col.names = NA)
        }
      )
      
      
      
      
      output$spit <-
        
        renderPrint(summary(dds.results()))
      
      
      
      
      ###Gene clusters
      
      groups <- names(data.frame(colData(dds)) %>% select(!c('sample','group')))
      output$within_group <- renderUI({
        
        
        tagList(
          selectInput(inputId = "within_choice", label = "Select within-group (only multi-factor experimental design)",
                      choices = groups)
        )
      })
      
      output$within_i <- renderUI({
        
        
        tagList(
          selectInput(inputId = "within_choice_i", label = "Select individual group",
                      choices = data.frame(colData(dds)) %>% pull(input$within_choice) %>% unique() )
        )
      })
      
      output$between_group <- renderUI({
        
        
        tagList(
          selectInput(inputId = "between_choice", label = "Select between-group",
                      choices = groups)
        )
      })
      
      output$gene_filter <- renderUI({
        
        
        tagList(
          selectInput(inputId = "gfilter_choice", label = "Select filtering method for cluster analysis",
                      choices = c('Top DE Genes', 'Top Var Genes'), selected = 'Top Var Genes')
        )
      })
      
      output$topn <- renderUI({
        
        tagList(
          selectInput(inputId='topnchoice', label='Select top N genes to cluster',
                      choices = c(100,250,500), selected=250)
        )
        
      })
      
      output$order <- renderUI({
        rank_list(text='Select order of between-groups',
                  labels=data.frame(colData(dds)) %>% pull(input$between_choice) %>% unique(),
                  input_id = 'btwnorder')
      })
      
      
      ###Print Heatmap
      observeEvent(input$run_heat, {
        req(input$between_choice)
        print(input$btwnorder)
        print(round(as.numeric(input$topnchoice) / 2, 0))
        #prepare for cluster
        meta <- data.frame(colData(dds))
        vst <- assay(vst(dds))
        samps <- meta %>% filter(get(input$within_choice)==input$within_choice_i) %>% pull(sample)
        rvars <- rowVars(vst[,which(colnames(vst) %in% samps)], useNames = T) 
        ngenes <- data.frame('gene'=names(rvars), 'vars'=rvars) %>% 
          mutate(idx=row_number()) %>% 
          arrange(-vars) %>% head(n=as.numeric(input$topnchoice)) %>% 
          pull(gene)
        vst_i <- vst[ngenes,which(colnames(vst) %in% samps)]
        vst_i <- vst_i[!is.na(rowVars(vst_i, useNames = T)),]
        vst_i <- vst_i[rowVars(vst_i, useNames = T)>0,]
        print(dim(vst_i))
        
        ###Cluster
        pamClusters <- cluster::pam(vst_i, k = 3) # pre-select k = 9 centers
        pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
        pamClusters$clustering <- factor(pamClusters$clustering,
                                         levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 
                                                    'Cluster 5', 'Cluster 6', 'Cluster 7', 'Cluster 8', 'Cluster 9'))
        
        ###gen heatmap
        mcolors<- list(temp=c('SummerActive'="#A50F15",
                              '4C'="#FEE5D9",
                              '12C'="#FCBBA1",
                              '20C'="#FC9272",
                              '25C'="#FB6A4A",
                              '30C'="#DE2D26"))
        
        colAnn <- HeatmapAnnotation(df= meta %>%  
                                      filter(sample %in% colnames(vst_i)) %>% 
                                      select(input$between_choice) %>% 
                                      arrange(factor(input$between_choice, levels=input$btwnorder)) %>% 
                                               as.data.frame(),
                                    which='col',
                                    na_col='white',
                                    col=mcolors, simple_anno_size = unit(3,'cm'),
                                    annotation_name_gp= gpar(fontsize = 35)) 
        boxplotCol <- HeatmapAnnotation(
          boxplot = anno_boxplot(
            vst_i,
            border = FALSE,
            gp = gpar(fill = '#CCCCCC'),
            pch = '.',
            size = unit(2, 'mm'),
            axis = TRUE,
            axis_param = list(
              gp = gpar(fontsize = 12),
              side = 'left')),
          annotation_width = unit(c(2.0), 'cm'),
          which = 'col')
        
        
        boxplotRow <- HeatmapAnnotation(
          boxplot = row_anno_boxplot(
            vst_i,
            border = FALSE,
            gp = gpar(fill = '#CCCCCC'),
            pch = '.',
            size = unit(2, 'mm'),
            axis = TRUE,
            axis_param = list(
              gp = gpar(fontsize = 12),
              side = 'top')),
          annotation_width = unit(c(2.0), 'cm'),
          which = 'row')
        
        nsamp <- round(as.numeric(input$topnchoice) / 100, 0)
        genelabels <- rowAnnotation(
          Genes = anno_mark(
            at = seq(1, nrow(vst_i), nsamp),
            labels = rownames(vst_i)[seq(1, nrow(vst_i), nsamp)],
            labels_gp = gpar(fontsize = 6, fontface = 'bold'),
            padding = 0.75),
          width = unit(2.0, 'cm') +
            
            max_text_width(
              rownames(vst_i)[seq(1, nrow(vst_i), nsamp)],
              gp = gpar(fontsize = 6,  fontface = 'bold')))
        
        myCol <- colorRampPalette(c('#ea6bdd', 'black','#50eea7'))(100)
        myBreaks <- seq(-2.5, 2.5, length.out = 100)
        
        #order samples by between group
        mmeta <- meta
        rownames(mmeta) <- c()
        print(input$between_choice)
        sorder <- mmeta %>% arrange(factor(get(input$between_choice), levels=input$btwnorder)) %>% pull(sample)
        
        hmap <- Heatmap(t(scale(t(vst_i))),
                        
                        # split the genes / rows according to the PAM clusters
                        split = pamClusters$clustering,
                        cluster_row_slices = FALSE,
                        
                        name = 'Gene\nZ-\nscore',
                        
                        col = colorRamp2(myBreaks, myCol),
                        
                        # parameters for the colour-bar that represents gradient of expression
                        heatmap_legend_param = list(
                          color_bar = 'continuous',
                          legend_direction = 'vertical',
                          legend_width = unit(7, 'cm'),
                          legend_height = unit(5.0, 'cm'),
                          title_position = 'topcenter',
                          title_gp=gpar(fontsize = 20, fontface = 'bold'),
                          labels_gp=gpar(fontsize = 20, fontface = 'bold')),
                        
                        # row (gene) parameters
                        cluster_rows = TRUE,
                        show_row_dend = TRUE,
                        #row_title = 'Statistically significant genes',
                        row_title_side = 'left',
                        row_title_gp = gpar(fontsize = 18,  fontface = 'bold'),
                        row_title_rot = 90,
                        show_row_names = FALSE,
                        row_names_gp = gpar(fontsize = 14, fontface = 'bold'),
                        row_names_side = 'left',
                        row_dend_width = unit(55,'mm'),
                        
                        # column (sample) parameters
                        cluster_columns = FALSE,
                        column_order = sorder[sorder %in% colnames(vst_i)],
                        show_column_dend = F,
                        column_title = '',
                        column_title_side = 'bottom',
                        column_title_gp = gpar(fontsize = 22, fontface = 'bold'),
                        column_title_rot = 0,
                        show_column_names = T,
                        column_names_gp = gpar(fontsize = 20, fontface = 'bold'),
                        column_names_max_height = unit(20, 'cm'),
                        column_dend_height = unit(25,'mm'),
                        
                        # cluster methods for rows and columns
                        clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                        clustering_method_columns = 'ward.D2',
                        clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                        clustering_method_rows = 'ward.D2',
                        
                        # specify top and bottom annotations
                        top_annotation = colAnn,
                        bottom_annotation = boxplotCol)
        
        output$hmap_plot <- renderPlot(
          draw(hmap + genelabels,
             heatmap_legend_side = 'left',
             annotation_legend_side = 'right',
             row_sub_title_side = 'left'), 
          height = 1200, 
          width = 1200
        )
      })
      
      
      
  })
  
    
  
  
  
  
})