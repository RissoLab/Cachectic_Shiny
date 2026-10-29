# ===================================================================================
# CachexiaAtlas
# Interactive tool for exploring gene expression across cell populations in cachexia
# ===================================================================================

library(shiny)
library(SingleCellExperiment)
library(ggplot2)
library(viridis)
library(plotly)
library(here)

load(here("data", "palette_color.RData"))
sce <- readRDS(here("data", "shiny_app_sce.rds"))

# UMAP data configuration 

coords <- as.data.frame(reducedDim(sce, "UMAP"))
coords <- cbind(coords, colData(sce))
coords <- coords[, -6]
coords$Condition <- factor(coords$Condition, levels = c("CTRL", "C26"))

metadata <- colnames(colData(sce)[-4])
# assay_names <- names(assays(sce))
genes <- rownames(sce)

# user interface ######
ui <- fluidPage(
  titlePanel("CachexiaAtlas: visualizaion of gene expression in cachexia condition"),
  tags$p("CachexiaAtlas is an interactive tool for exploring gene expression across cell populations in cachexia.
         Visualize cell identities and gene expression patterns using UMAP and violin plots, and investigate expression differences between cell types.
         If you use CachexiaAtlas in your work, please cite: Pezzini C., Persico G., Zuin E., et al. 2026 (Cell Reports Medicine).", 
         style = paste("font-size: 16px;", "color: #555555;","line-height: 1.6;", "margin-bottom: 25px;")),
  
  sidebarPanel(
    selectInput(inputId = "group", label = "Color by:", choices = metadata, selected = "Cell Type"),
    selectizeInput(inputId = "gene", label = "Select gene:", choices = genes, selected = "Trim63"),
    # selectInput(inputId = "assay", label = "Expression assay:", choices = assay_names, selected = "logcounts"),
    
    hr(),
    h4("Plot settings"),
    sliderInput(inputId = "point_size", label = "UMAP point size:", min = 0.5, max = 2, value = 1, step = 0.5),
    sliderInput(inputId = "point_alpha", label = "Point transparency:", min = 0.1, max = 1, value = 0.8, step = 0.1),
    
    hr(),
    h4("Download figures"),
    downloadButton(outputId = "download_umap_ct", label = "Download UMAP", class = "btn-primary btn-block"),
    br(),
    downloadButton(outputId = "download_umap_expr", label = "Download UMAP (gene expression)", class = "btn-primary btn-block"),
    br(),
    downloadButton(outputId = "download_violin", label = "Download violin plot", class = "btn-primary btn-block"),
    br(),
    downloadButton(outputId = "download_all", label = "Download all figures (ZIP)", class = "btn-success btn-block"),
    width = 3),
  
  # Main panel
  
  mainPanel(
    # First row
    fluidRow(column(width = 12, plotlyOutput(outputId = "umap_celltype", height = "650px"))),
    hr(),
    # Second row
    fluidRow(column(width = 6, plotlyOutput(outputId = "umap_expression", height = "650px")),
             column(width = 6, plotlyOutput(outputId = "violin_plot", height = "650px")))
  ))

# server #####

server <- function(input, output, session) {
  
  expression_data <- reactive({
    # req(input$gene, input$assay)
    # validate(need(input$gene %in% genes, "Gene not found in the dataset."))
    # values <- assay(sce, input$assay)[input$gene, ]
    values <- assay(sce, "logcounts")[input$gene, ]
    as.numeric(values)
  })  
  
  # UMAP: cell type #####
  
  umap_ct <- reactive({
    req(input$group)
    
    df <- coords
    df$Group <- factor(as.character(df[[input$group]]))
    
    ggplot(df, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = Group,
                   text = paste0( "UMAP 1: ", round(wnnUMAP_1, 2),
                                  "<br>UMAP 2: ", round(wnnUMAP_2, 2),
                                  "<br>", input$group, ": ", Group))) + 
      geom_point(size = input$point_size, alpha = input$point_alpha) +
      scale_color_manual(values = palette_color) + facet_wrap(~Condition) +
      theme_bw() + theme(text = element_text(size = 12), legend.text = element_text(size = 12)) +
      labs(color = input$group)
  })
  
  output$umap_celltype <- renderPlotly({
    ggplotly(umap_ct(), tooltip = "text") %>%
      layout(hoverlabel = list(bgcolor = "white", font = list(size = 14)))
  })
  
  # UMAP: gene expression ####
  
  umap_exp <- reactive({
    req(input$gene)
    
    df <- coords
    df$Expression <- expression_data()
    df$CellType <- as.factor(coords$`Cell Type`)
    df$Condition <- as.factor(coords$Condition)

    ggplot(df, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = Expression,
                        text = paste0("UMAP 1: ", round(wnnUMAP_1, 2),
                                      "<br>UMAP 2: ", round(wnnUMAP_2, 2),
                                      "<br>Gene: ", input$gene,
                                      "<br>Expression: ", round(Expression, 2),
                                      "<br>Cell Type: ", CellType,
                                      "<br>Condition: ", Condition))) +
      geom_point(size = input$point_size, alpha = input$point_alpha) +
      scale_color_gradient(low = "grey90", high = "darkblue", name = paste0(input$gene, " \nExpression")) +
      theme_bw() + theme(text = element_text(size = 12), legend.text = element_text(size = 12))
  })
  
  output$umap_expression <- renderPlotly({
    
    ggplotly(umap_exp(), tooltip = "text") %>%
      layout(hoverlabel = list(bgcolor = "white", font = list(size = 14)))
  })
  
  # Violin plot #### 
  vln_plot <- reactive({
    req(input$gene)
    
    plot_df <- data.frame(
      Expression = expression_data(),
      CellType = as.factor(coords$`Cell Type`),
      Condition = coords$Condition
    )
    
    plot_df <- plot_df %>%
      dplyr::group_by(CellType, Condition) %>%
      dplyr::mutate(
        mean = mean(Expression, na.rm = TRUE),
        median = median(Expression, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
    
    ggplot(plot_df, aes(x = CellType, y = Expression, fill = CellType, 
                        text = paste0("Mean expression: ", round(mean, 2),
                                      "<br>Median expression: ", round(median, 2)))) +
      geom_violin(trim = TRUE, scale = "width") +
      theme_bw() + scale_fill_manual(values = palette_color) + facet_grid(Condition ~ .) +
      theme(text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title.x = element_text(margin = margin(t = 30)),
            legend.position = "none") + labs(title = input$gene, x = "Cell Type")
  })
  
  output$violin_plot <- renderPlotly({
    ggplotly(vln_plot(), tooltip = "text") %>%
      layout(hoverlabel = list(bgcolor = "white", font = list(size = 14)))
  })
  
  # Downaload figures ####
  
  output$download_umap_ct <- downloadHandler(
    
    filename = function() {
      "cachexiaAtlas_umap_celltype.pdf"
    },
    
    content = function(file) {
      ggsave(filename = file, plot = umap_ct(), device = "pdf", width = 12, 
             height = 8, units = "in", dpi = 300, bg = "white")
    }
  )
  
  output$download_umap_expr <- downloadHandler(
    
    filename = function() {
      paste0("cachexiaAtlas_", input$gene, "_umap_expression.pdf")
      
    },
    
    content = function(file) {
      ggsave(filename = file, plot = umap_exp(), device = "pdf", width = 10,
             height = 8, units = "in", dpi = 300, bg = "white")
    }
  )
  
  output$download_violin <- downloadHandler(
    filename = function() {
      paste0("cachexiaAtlas_", input$gene, "_violin_plot.pdf")
    },
    
    content = function(file) {
      ggsave(filename = file, plot = vln_plot(), device = "pdf", width = 12, 
             height = 8, units = "in", dpi = 300, bg = "white")
    }
  )
  
  output$download_all <- downloadHandler(
    
    filename = function() {
      paste0("cachexiaAtlas_", input$gene, "_figures.zip")
    },
    
    content = function(file) {
      temp_dir <- tempdir()
      plot_dir <- file.path(temp_dir, "cachexiaAtlas_figures")
      dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
      
      file_umap_ct <- file.path(plot_dir, "cachexiaAtlas_umap_celltype.pdf")
      file_umap_expr <- file.path(plot_dir, paste0("cachexiaAtlas_", input$gene, "_umap_expression.pdf"))
      file_violin <- file.path(plot_dir, paste0("cachexiaAtlas_", input$gene, "_violin_plot.pdf"))
      
      ggsave(filename = file_umap_ct,  plot = umap_ct(), device = "pdf",
             width = 12, height = 8, units = "in", dpi = 300, bg = "white")
      
      ggsave(filename = file_umap_expr, plot = umap_exp(), device = "pdf",
             width = 10, height = 8, units = "in", dpi = 300, bg = "white")
      
      ggsave(filename = file_violin, plot = vln_plot(), device = "pdf", width = 12,
             height = 8, units = "in", dpi = 300, bg = "white")
      
      zip(zipfile = file, files = c(file_umap_ct, file_umap_expr, file_violin),
          flags = "-j")
    }
  )
}

shinyApp(
  ui = ui,
  server = server
)

