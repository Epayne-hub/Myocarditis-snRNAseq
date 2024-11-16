#' Interactive visualization of Seurat (v3) results
#'
#' @param seurat Seurat object
#' @param appTitle App title
#' @param ... Additional arguments (currently not used)
#'
#' @author Charlotte Soneson/Mechthild Lütge/Roman Stadler
#'
#' @import shiny
#' @import ggplot2
#' @import plotly
#' @import shinycssloaders
#' @importFrom dplyr mutate group_by top_n summarise filter
#' @importFrom readr read_tsv
#' @importFrom magrittr %>%
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#'   dashboardBody tabBox
#' @importFrom pheatmap pheatmap
#' @importFrom Cairo Cairo
#' @importFrom colourpicker colourInput
#' @export
#'
SingleCellBrowserSeurat3 <- function(seurat, appTitle = "SingleCellBrowser", ...) {
  require(plotly)
  require(shinycssloaders)
  
  if (packageVersion("Seurat") < "3.0") {
    stop("Seurat Version < 3.0 is loaded. Function needs package Seurat Version >= 3.0")
  }
  
  gg_fill_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  pLayout <- function() {
    shinydashboard::dashboardPage(
      skin = "purple",
      
      shinydashboard::dashboardHeader(title = appTitle, titleWidth = nchar(appTitle) * 20),
      
      shinydashboard::dashboardSidebar(),
      
      shinydashboard::dashboardBody(shiny::fluidRow(do.call(
        shinydashboard::tabBox,
        c(
          width = 12,
          
          ## ======================================================== ##
          ## dimension reduction: clusters
          ## ======================================================== ##
          list(
            shiny::tabPanel(
              "Dimensional reduction: color by cluster",
              shiny::fluidRow(
                shiny::column(
                  3,
                  shiny::selectInput(
                    inputId = "tsne.cluster",
                    choices = c(colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)])),
                    label = "Color by",
                    selected = "dataset",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  2,
                  offset = 1,
                  shiny::selectInput(
                    inputId = "dr.method",
                    choices = c("tsne", "umap"),
                    label = "DR technique",
                    selected = "umap",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  2,
                  offset = 1,
                  shiny::selectInput(
                    inputId = "split.grp",
                    choices = c(colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)])),
                    label = "split plot by",
                    selected = NULL,
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                # shiny::column(
                #   2, offset = 1, shiny::selectInput(inputId = "legend",
                #                                     choices = c("right","none", "left", "top", "bottom"),
                #                                     label = "legend position",
                #                                     selected = "right",
                #                                     selectize = TRUE, multiple = FALSE)
                # ),
                shiny::column(
                  2,
                  offset = 1,
                  shiny::selectInput(
                    inputId = "pointsize",
                    choices = c("0.001", "0.01", "0.1", "0.5", "1", "2", "3"),
                    label = "Point Size",
                    selected = "0.1",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                # shiny::column(
                #   2, offset = 3, shiny::downloadButton(outputId = "export_TSNE",
                #                                        label = "export as PDF")
                # )
              ),
              # shiny::fluidRow(
              #   shiny::column(
              #     6, shiny::selectInput(inputId = "colGroup",
              #                           choices = seq(from=0, to=30, by=1),
              #                           label = "choose group",
              #                           selected =levels(seurat$dataset),
              #                           selectize = TRUE, multiple = TRUE)
              #   )
              # ),
              shiny::uiOutput("tsne.by.cluster.plot.ui")
            )
          ),
          
          ## ======================================================== ##
          ## dimension reduction: genes
          ## ======================================================== ##
          
          list(
            shiny::tabPanel(
              "Dimensional reduction: color by gene",
              shiny::fluidRow(
                shiny::column(
                  3,
                  shiny::selectInput(
                    inputId = "tsne.gene",
                    choices = rownames(seurat),
                    label = "Selected gene",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  2,
                  offset = 1,
                  shiny::selectInput(
                    inputId = "dr.method.feature",
                    choices = c("tsne", "umap"),
                    label = "DR technique",
                    selected = "umap",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  2,
                  offset = 1,
                  shiny::selectInput(
                    inputId = "split.grp2",
                    choices = c(colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)])),
                    label = "split plot by",
                    selected = NULL,
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  2,
                  offset = 1,
                  shiny::selectInput(
                    inputId = "pointsize",
                    choices = c("0.001", "0.01", "0.1", "0.5", "1", "2", "3"),
                    label = "Point Size",
                    selected = "0.1",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                )
                # shiny::column(
                #   2, offset = 2, shiny::downloadButton(outputId = "export_TSNEgene",
                #                                        label = "export as PDF")
                # )
              ),
              shiny::uiOutput("tsne.by.gene.plot.ui")
            )
          ),
          
          ## ======================================================== ##
          ## Tab with column data distribution
          ## ======================================================== ##
          list(
            shiny::tabPanel(
              "Cell annotations",
              shiny::fluidRow(
                shiny::column(
                  6,
                  shiny::selectInput(
                    inputId = "col.choice",
                    choices = colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)]),
                    label = "Cell annotation",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  6,
                  shiny::selectInput(
                    inputId = "color.choice",
                    choices = colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)]),
                    selected = "dataset",
                    label = "Cell coloring",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                # shiny::column(
                #   2, shiny::downloadButton(outputId = "export_cellAnno",
                #                            label = "export as PDF")
                # )
                
              ),
              
              # shiny::fluidRow(
              #   shiny::column(
              #     6, shiny::selectInput(inputId = "colGroupAnno",
              #                           choices = seq(from=0, to=30, by=1),
              #                           label = "choose group",
              #                           selected =levels(seurat$dataset),
              #                           selectize = TRUE, multiple = TRUE)
              #   )
              # ),
              shiny::uiOutput("cell.annotation.plot.ui")
            )
          ),
          
          ## ======================================================== ##
          ## Tab with individual expression levels
          ## ======================================================== ##
          list(
            shiny::tabPanel(
              "Individual gene expression",
              shiny::fluidRow(
                shiny::column(
                  6,
                  shiny::selectInput(
                    inputId = "expression.group",
                    choices = c(colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)])),
                    label = "Group by",
                    selected = "dataset",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  6,
                  shiny::selectInput(
                    inputId = "expression.genes",
                    choices = rownames(seurat),
                    label = "Selected gene",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                # shiny::column(
                #   2, shiny::downloadButton(outputId = "export_vlnPlot",
                #                            label = "export as PDF")
                # )
                
              ),
              
              # shiny::fluidRow(
              #   shiny::column(
              #     6, shiny::selectInput(inputId = "colGroupExpr",
              #                           choices = seq(from=0, to=30, by=1),
              #                           label = "choose group",
              #                           selected =levels(seurat$dataset),
              #                           selectize = TRUE, multiple = TRUE)
              #   )
              # ),
              shiny::uiOutput("expression.plot.ui")
            )
          ),
          
          ## ======================================================== ##
          ## Tab with heatmap
          ## ======================================================== ##
          list(
            shiny::tabPanel(
              "Heatmap",
              shiny::fluidRow(
                shiny::column(
                  3,
                  shiny::selectInput(
                    inputId = "heatmap.group",
                    choices = colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)]),
                    label = "Group by",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  6,
                  shiny::selectInput(
                    inputId = "heatmap.genes",
                    choices = rownames(seurat),
                    label = "Selected genes (at least 2)",
                    selectize = TRUE,
                    multiple = TRUE
                  )
                ),
                # shiny::column(
                #   2, shiny::downloadButton(outputId = "export_heatPlot",
                #                            label = "export as PDF")
                # )
              ),
              shiny::uiOutput("heatmap.plot.ui")
            )
          ),
          
          ## ======================================================== ##
          ## Tab with marker genes
          ## ======================================================== ##
          list(
            shiny::tabPanel(
              "Marker Genes",
              shiny::fluidRow(
                shiny::column(
                  8,
                  fileInput(
                    "MarkerGenes",
                    "load marker gene list",
                    multiple = FALSE,
                    accept = ".txt",
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  )
                ),
                shiny::column(
                  2,
                  shiny::downloadButton(outputId = "export_markerPlot", label = "export as PDF")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  4,
                  shiny::selectInput(
                    inputId = "cluster.group",
                    choices = levels(seurat),
                    selected = levels(seurat),
                    label = "Selected cluster",
                    selectize = TRUE,
                    multiple = TRUE
                  )
                  
                ),
                shiny::column(
                  3,
                  shiny::numericInput(
                    inputId = "adjPval_Marker_cut",
                    value = 0.01,
                    min = 0,
                    max = 1,
                    step = 0.005,
                    label = "adjusted p-value cutoff"
                  )
                  
                ),
                shiny::column(
                  2,
                  shiny::numericInput(
                    inputId = "marker.count.start",
                    value = 0,
                    min = 0,
                    max = 500,
                    label = "From"
                  )
                  
                ),
                shiny::column(
                  2,
                  shiny::numericInput(
                    inputId = "marker.count.end",
                    value = 10,
                    min = 1,
                    max = 500,
                    label = "To"
                  )
                )
              ),
              shiny::uiOutput("marker.plot.ui")
            )
          ),
          
          ## ======================================================== ##
          ## Tab with specific genes
          ## ======================================================== ##
          list(
            shiny::tabPanel(
              "Specific Genes",
              shiny::fluidRow(
                shiny::column(
                  8,
                  fileInput(
                    "specificGenes",
                    "load specific gene list",
                    multiple = FALSE,
                    accept = ".txt",
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  )
                ),
                shiny::column(
                  2,
                  shiny::downloadButton(outputId = "export_specificPlot", label = "export as PDF")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  6,
                  shiny::selectInput(
                    inputId = "cluster.group2",
                    choices = levels(seurat),
                    label = "Selected cluster",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                  
                ),
                shiny::column(
                  6,
                  shiny::sliderInput(
                    "geneRangeList",
                    label = "gene range from list",
                    min = 0,
                    max = 250,
                    value = c(0, 15)
                  )
                )
              ),
              shiny::uiOutput("specificGenes.plot.ui")
            )
          ),
          
          ## ======================================================== ##
          ## Tab with DE genes between 2 conditions
          ## ======================================================== ##
          list(
            shiny::tabPanel(
              "DE Genes",
              shiny::fluidRow(
                shiny::column(
                  8,
                  fileInput(
                    "DEgenes",
                    "load DE gene list",
                    multiple = FALSE,
                    accept = ".txt",
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  )
                ),
                shiny::column(
                  2,
                  shiny::downloadButton(outputId = "export_DEPlot", label = "export as PDF")
                )
              ),
              shiny::fluidRow(
                shiny::column(
                  3,
                  shiny::selectInput(
                    inputId = "DEgenes.group",
                    choices = c(colnames(seurat@meta.data[, apply(seurat@meta.data, 2, function(x)
                      class(x) != "numeric" && length(unique(x)) <= 50)])),
                    label = "Group by",
                    selected = "dataset",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                shiny::column(
                  3,
                  shiny::selectInput(
                    inputId = "DEgenes.cluster",
                    choices = levels(seurat),
                    label = "Choose cluster",
                    selectize = TRUE,
                    multiple = FALSE
                  )
                ),
                
                shiny::column(
                  2,
                  shiny::numericInput(
                    inputId = "adjPval_cut",
                    value = 0.01,
                    min = 0,
                    max = 1,
                    step = 0.005,
                    label = "adjusted p-value cutoff"
                  )
                  
                ),
                shiny::column(
                  2,
                  shiny::numericInput(
                    inputId = "top.count.start",
                    value = 0,
                    min = 0,
                    max = 500,
                    label = "From"
                  )
                  
                ),
                shiny::column(
                  2,
                  shiny::numericInput(
                    inputId = "top.count.end",
                    value = 10,
                    min = 1,
                    max = 500,
                    label = "To"
                  )
                )
                
              ),
              shiny::uiOutput("DE.plot.ui")
            )
          )
          
          ###### add new tabs before here
        )
      )))
    )
  }
  
  server_function <- function(input, output, session) {
    options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)
    
    ## ====================================================================== ##
    ## Dimension reduction plots: By Cluster ID
    ## ====================================================================== ##
    
    output$tsne.by.cluster.plot.ui <- shiny::renderUI({
      lev <- sort(unique(input$colGroup)) # sorting so that "things" are unambiguous
      cols <- gg_fill_hue(length(lev))
      
      tagList(
        shiny::fluidRow(lapply(seq_along(lev), function(i) {
          shiny::column(2,
                        colourInput(
                          inputId = paste0("col", lev[i]),
                          label = paste0("Choose colour for group ", lev[i]),
                          value = cols[i]
                        ))
          
        })),
        plotly::plotlyOutput(
          "tsne.by.cluster.plot",
          width = "90%",
          height = "600px"
        )   %>% shinycssloaders::withSpinner(color = "#0dc5c1", type = 6)
      )
    })
    
    output$tsne.by.cluster.plot <- plotly::renderPlotly({
      colInp <- paste0("c(", paste0("input$col", sort(input$colGroup), collapse = ", "), ")")
      cols <- eval(parse(text = colInp))
      colVec <- c(
        ggpubr::get_palette("npg", 12),
        ggpubr::get_palette("uchicago", 8),
        ggpubr::get_palette("aaas", 10)
      )
      colVec <- c(cols, colVec)
      
      p_tsne <- Seurat::DimPlot(
        object = seurat,
        reduction = input$dr.method,
        group.by = input$tsne.cluster,
        pt.size = input$pointsize,
        split.by = input$split.grp,
        cols = colVec,
        raster = FALSE
      ) +
        theme(legend.position = input$legend)
      
      # Convert to plotly
      plotly::ggplotly(p_tsne) %>% toWebGL() # improve performance
    })
    
    
    # export output to pdf
    output$export_TSNE <- downloadHandler("tsne_metavar.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 29,
          height = 22,
          units = 'cm',
          bg = 'transparent'
        )
        colInp <- paste0("c(", paste0("input$col", sort(input$colGroup), collapse = ", "), ")")
        cols <- eval(parse(text = colInp))
        colVec <- c(
          ggpubr::get_palette("npg", 12),
          ggpubr::get_palette("uchicago", 8),
          ggpubr::get_palette("aaas", 10)
        )
        colVec <- c(cols, colVec)
        
        p_tsne <- Seurat::DimPlot(
          object = seurat,
          reduction = input$dr.method,
          group.by = input$tsne.cluster,
          pt.size = 2,
          split.by = input$split.grp,
          cols = colVec
        ) +
          theme(legend.position = input$legend)
        
        print(p_tsne)
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    ## ====================================================================== ##
    ## Dimension reduction plots: By Gene Expression
    ## ====================================================================== ##
    
    ## Plotly Version:
    output$tsne.by.gene.plot.ui <- shiny::renderUI({
      plotly::plotlyOutput("tsne.by.gene.plot",
                           width = "90%",
                           height = "600px")  %>% shinycssloaders::withSpinner(color = "#0dc5c1", type =
                                                                                 6)
    })
    
    output$tsne.by.gene.plot <- plotly::renderPlotly({
      # Create FeaturePlot with split.by to generate separate ggplot objects for each specified group
      plots <- FeaturePlot(
        seurat,
        features = input$tsne.gene,
        pt.size = input$pointsize,
        reduction = input$dr.method.feature,
        split.by = input$split.grp2,
        cols = c("azure3", "red3"),
        raster = FALSE
      )
      
      # n_plots <- length(plots)
      # max_cols <- 5  # Set the desired number of columns
      # nrows <- ceiling(n_plots / max_cols)  # Calculate number of rows
      
      factor_levels <- levels(as.factor((seurat@meta.data[[input$split.grp2]])))
      
      plotly_plots <- lapply(seq_along(plots), function(i) {
        plot <- plots[[i]]  +
          ggtitle("") +
          xlab(paste("umap_1", "\n", factor_levels[i], collapse = ", ")) # cannot merge titles, so I move it to the x label
      })
      
      # Combine the plotly objects in a grid layout
      plotly::subplot(
        plotly_plots,
        shareX = T,
        shareY = T,
        titleX = T,
        titleY = T,
        margin = 0.05
      ) %>% toWebGL()
    })
    
    # export output to pdf
    output$export_TSNEgene <- downloadHandler("tsne_gene.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 29,
          height = 22,
          units = 'cm',
          bg = 'transparent'
        )
        pFeature <- Seurat::FeaturePlot(
          object = seurat,
          features = input$tsne.gene,
          cols = c("azure3", "red3"),
          reduction = input$dr.method.feature,
          pt.size = input$pointsize,
          split.by = input$split.grp2
        )
        print(pFeature)
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    
    # ### Static Version:
    # output$tsne.by.gene.plot.ui <- shiny::renderUI({
    #   shiny::plotOutput("tsne.by.gene.plot", width = "90%", height = "600px")
    # })
    
    # output$tsne.by.gene.plot <- shiny::renderPlot({
    #   Seurat::FeaturePlot(object = seurat, features = input$tsne.gene,
    #                       cols = c("azure3","red3"), reduction = input$dr.method.feature,
    #                       pt.size = 1, split.by = input$split.grp2)
    # })
    
    # # export output to pdf
    # output$export_TSNEgene <- downloadHandler("tsne_gene.pdf", function(theFile) {
    #   makePdf <- function(filename){
    #     Cairo(type = 'pdf', file = filename, width = 29, height = 22, units='cm', bg='transparent')
    #     pFeature <- Seurat::FeaturePlot(object = seurat, features = input$tsne.gene,
    #                         cols = c("azure3","red3"), reduction = input$dr.method.feature,
    #                         pt.size = input$pointsize, split.by = input$split.grp2)
    #     print(pFeature)
    #     dev.off()
    #   }
    
    #   makePdf(theFile)
    
    # })
    
    ## ====================================================================== ##
    ## Cell annotation
    ## ====================================================================== ##
    output$cell.annotation.plot.ui <- shiny::renderUI({
      lev <- sort(unique(input$colGroupAnno)) # sorting so that "things" are unambiguous
      cols <- gg_fill_hue(length(lev))
      
      tagList(
        shiny::fluidRow(lapply(seq_along(lev), function(i) {
          shiny::column(2,
                        colourInput(
                          inputId = paste0("colAnno", lev[i]),
                          label = paste0("Choose colour for group ", lev[i]),
                          value = cols[i]
                        ))
        })),
        plotly::plotlyOutput(
          "cell.annotation.plot",
          width = "90%",
          height = "600px"
        ) # updated to plotlyOutput
      )
    })
    
    output$cell.annotation.plot <- plotly::renderPlotly({
      # updated to renderPlotly
      
      colInp <- paste0("c(", paste0("input$colAnno", sort(input$colGroupAnno), collapse = ", "), ")")
      cols <- eval(parse(text = colInp))
      colVec <- c(
        ggpubr::get_palette("npg", 12),
        ggpubr::get_palette("uchicago", 8),
        ggpubr::get_palette("aaas", 10)
      )
      colVec <- c(cols, colVec)
      
      # Create ggplot object
      p_cellAnno <- if (class(seurat[[input$col.choice]]) == "numeric") {
        ggplot(
          seurat@meta.data,
          aes_string(
            x = input$col.choice,
            fill = input$color.choice
          )
        ) +
          geom_density() + theme_bw() + scale_fill_manual(values = colVec)
      } else {
        ggplot(
          seurat@meta.data,
          aes_string(
            x = input$col.choice,
            fill = input$color.choice
          )
        ) +
          geom_bar() + theme_bw() + scale_fill_manual(values = colVec)
      }
      
      # Convert ggplot to Plotly
      plotly::ggplotly(p_cellAnno) %>% toWebGL() # improve performance
    })
    
    # Export output to PDF
    output$export_cellAnno <- downloadHandler("cellAnno.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 22,
          height = 27,
          units = 'cm',
          bg = 'transparent'
        )
        
        colInp <- paste0("c(",
                         paste0("input$colAnno", sort(input$colGroupAnno), collapse = ", "),
                         ")")
        cols <- eval(parse(text = colInp))
        colVec <- c(
          ggpubr::get_palette("npg", 12),
          ggpubr::get_palette("uchicago", 8),
          ggpubr::get_palette("aaas", 10)
        )
        colVec <- c(cols, colVec)
        
        # Create ggplot object
        p_cellAnno <- if (class(seurat[[input$col.choice]]) == "numeric") {
          ggplot(
            seurat@meta.data,
            aes_string(
              x = input$col.choice,
              fill = input$color.choice
            )
          ) +
            geom_density() + theme_bw() + scale_fill_manual(values = colVec)
        } else {
          ggplot(
            seurat@meta.data,
            aes_string(
              x = input$col.choice,
              fill = input$color.choice
            )
          ) +
            geom_bar() + theme_bw() + scale_fill_manual(values = colVec)
        }
        
        print(p_cellAnno)
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    
    ## ====================================================================== ##
    ## Abundance plots
    ## ====================================================================== ##
    output$expression.plot.ui <- shiny::renderUI({
      lev <- sort(unique(input$colGroupExpr)) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))
      
      # New IDs "colX1" so that it partly coincide with input$select...
      tagList(
        shiny::fluidRow(lapply(seq_along(lev), function(i) {
          shiny::column(2,
                        colourInput(
                          inputId = paste0("colExpr", lev[i]),
                          label = paste0("Choose colour for group ", lev[i]),
                          value = cols[i]
                        ))
          
        })),
        
        plotly::plotlyOutput("expression.plot", width = "90%", height = "500px")
      )
    })
    
    output$expression.plot <- plotly::renderPlotly({
      colInp <- paste0("c(", paste0("input$colExpr", sort(input$colGroupExpr), collapse = ", "), ")")
      cols <- eval(parse(text = colInp))
      colVec <- c(
        ggpubr::get_palette("npg", 12),
        ggpubr::get_palette("uchicago", 8),
        ggpubr::get_palette("aaas", 10)
      )
      colVec <- c(cols, colVec)
      
      # Create VlnPlot object
      pVln <- Seurat::VlnPlot(
        object = seurat,
        features = input$expression.genes,
        group.by = input$expression.group,
        cols = colVec,
        pt.size = 0
      )
      
      # Convert to plotly
      plotly::ggplotly(pVln) %>% toWebGL() # improve performance
    })
    
    # export output to pdf
    output$export_vlnPlot <- downloadHandler("vlnPlot.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 27,
          height = 22,
          units = 'cm',
          bg = 'transparent'
        )
        
        colInp <- paste0("c(",
                         paste0("input$colExpr", sort(input$colGroupExpr), collapse = ", "),
                         ")")
        cols <- eval(parse(text = colInp))
        colVec <- c(
          ggpubr::get_palette("npg", 12),
          ggpubr::get_palette("uchicago", 8),
          ggpubr::get_palette("aaas", 10)
        )
        colVec <- c(cols, colVec)
        
        pVln <- Seurat::VlnPlot(
          object = seurat,
          features = input$expression.genes,
          group.by = input$expression.group,
          cols = colVec,
          pt.size = 0
        )
        print(pVln)
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    
    ## ====================================================================== ##
    ## Heatmaps
    ## ====================================================================== ##
    output$heatmap.plot.ui <- shiny::renderUI({
      plotly::plotlyOutput("heatmap.plot", width = "90%", height = "600px") # updated to plotlyOutput
    })
    
    output$heatmap.plot <- plotly::renderPlotly({
      # updated to renderPlotly
      if (length(input$heatmap.genes) > 1) {
        p_heatmapSel <- Seurat::DoHeatmap(
          object = seurat,
          features = input$heatmap.genes,
          group.by = input$heatmap.group
        ) +
          scale_fill_distiller(palette = "RdBu")
        
        # Convert to plotly
        plotly::ggplotly(p_heatmapSel)
      } else {
        NULL
      }
    })
    
    # Export output to PDF
    output$export_heatPlot <- downloadHandler("heatmap.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 24,
          height = 22,
          units = 'cm',
          bg = 'transparent'
        )
        if (length(input$heatmap.genes) > 1) {
          p_heatmapSel <- Seurat::DoHeatmap(
            object = seurat,
            features = input$heatmap.genes,
            group.by = input$heatmap.group
          ) +
            scale_fill_distiller(palette = "RdBu")
          print(p_heatmapSel)
        } else {
          NULL
        }
        
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    
    ## ====================================================================== ##
    ## Marker Genes
    ## ====================================================================== ##
    
    
    output$marker.plot.ui <- shiny::renderUI({
      shiny::plotOutput("marker.plot", width = "90%", height = "600px")
    })
    
    dataMarkers <- reactive({
      validate(need(try(!is.null(input$MarkerGenes))
                    , "Please load marker gene file"))
      inFile <- input$MarkerGenes
      seurat_markers_all <- read_tsv(inFile$datapath, col_names = T)
      return(seurat_markers_all)
    })
    
    output$marker.plot <- shiny::renderPlot({
      seurat_markers <- dataMarkers() %>% group_by(cluster) %>% top_n(input$marker.count.end, avg_logFC) %>% top_n(-(input$marker.count.end -
                                                                                                                       input$marker.count.start),
                                                                                                                   avg_logFC) %>% as.data.frame() %>% dplyr::filter(cluster %in% input$cluster.group) %>% dplyr::filter(p_val_adj <= input$adjPval_Marker_cut)
      
      Seurat::DoHeatmap(object = seurat, features = seurat_markers$gene) +
        #group.by = levels(seurat)) +
        scale_fill_distiller(palette = "RdBu")
    })
    
    # export output to pdf
    output$export_markerPlot <- downloadHandler("markerPlot.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 24,
          height = 22,
          units = 'cm',
          bg = 'transparent'
        )
        seurat_markers <- dataMarkers() %>% group_by(cluster) %>% top_n(input$marker.count.end, avg_logFC) %>% top_n(-(input$marker.count.end -
                                                                                                                         input$marker.count.start),
                                                                                                                     avg_logFC) %>% as.data.frame() %>% filter(cluster %in% input$cluster.group) %>% filter(p_val_adj <= input$adjPval_Marker_cut)
        
        p_heatmap <- Seurat::DoHeatmap(object = seurat, features = seurat_markers$gene) +
          #group.by = levels(seurat)) +
          scale_fill_distiller(palette = "RdBu")
        print(p_heatmap)
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    ## ====================================================================== ##
    ## Specific Genes
    ## ====================================================================== ##
    
    
    output$specificGenes.plot.ui <- shiny::renderUI({
      shiny::plotOutput("specificGenes.plot",
                        width = "90%",
                        height = "600px")
    })
    
    dataSpecific <- reactive({
      validate(need(try(!is.null(input$specificGenes))
                    , "Please load specific gene list"))
      inFile2 <- input$specificGenes
      specififcGeneList <- read_tsv(inFile2$datapath, col_names = T) %>% filter(cluster1 %in% input$cluster.group2)
      return(specififcGeneList)
    })
    
    
    output$specificGenes.plot <- shiny::renderPlot({
      genes <- dataSpecific()
      
      Seurat::DoHeatmap(object = seurat, features = genes$gene[input$geneRangeList[1]:input$geneRangeList[2]]) +
        #group.by = input$heatmap.group) +
        scale_fill_distiller(palette = "RdBu")
      
    })
    
    # export output to pdf
    output$export_specificPlot <- downloadHandler("specificGenesPlot.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 24,
          height = 22,
          units = 'cm',
          bg = 'transparent'
        )
        genes <- dataSpecific()
        p_specific <- Seurat::DoHeatmap(object = seurat, features = genes$gene[input$geneRangeList[1]:input$geneRangeList[2]]) +
          #group.by = input$heatmap.group) +
          scale_fill_distiller(palette = "RdBu")
        print(p_specific)
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    
    ## ====================================================================== ##
    ## DE Genes
    ## ====================================================================== ##
    
    
    output$DE.plot.ui <- shiny::renderUI({
      shiny::plotOutput("DE.plot", width = "90%", height = "600px")
    })
    
    dataDE <- reactive({
      validate(need(try(!is.null(input$DEgenes))
                    , "Please load DE gene file"))
      inFile <- input$DEgenes
      DEgenes <- read_tsv(inFile$datapath, col_names = T)
      return(DEgenes)
    })
    
    output$DE.plot <- shiny::renderPlot({
      seurat_DE <- dataDE() %>% filter(p_val_adj <= input$adjPval_cut)
      if ("cluster" %in% colnames(seurat_DE)) {
        seurat_DE <- seurat_DE %>% filter(cluster == input$DEgenes.cluster)
      }
      
      Seurat::DoHeatmap(
        object = seurat,
        features = seurat_DE$gene[input$top.count.start:input$top.count.end],
        group.by = input$DEgenes.group
      ) +
        scale_fill_distiller(palette = "RdBu")
      
      
    })
    
    # export output to pdf
    output$export_DEPlot <- downloadHandler("DEGenesPlot.pdf", function(theFile) {
      makePdf <- function(filename) {
        Cairo(
          type = 'pdf',
          file = filename,
          width = 24,
          height = 22,
          units = 'cm',
          bg = 'transparent'
        )
        seurat_DE <- dataDE() %>% filter(p_val_adj <= input$adjPval_cut)
        if ("cluster" %in% colnames(seurat_DE)) {
          seurat_DE <- seurat_DE %>% filter(cluster == input$DEgenes.cluster)
        }
        
        p_DEgenes <- Seurat::DoHeatmap(
          object = seurat,
          features = seurat_DE$gene[input$top.count.start:input$top.count.end],
          group.by = input$DEgenes.group
        ) +
          scale_fill_distiller(palette = "RdBu")
        print(p_DEgenes)
        dev.off()
      }
      
      makePdf(theFile)
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
  
  
  
  
  
  shiny::shinyApp(ui = pLayout, server = server_function)
}
