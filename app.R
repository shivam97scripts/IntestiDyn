# Install required packages if not already installed
if (!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if (!requireNamespace("shinythemes", quietly = TRUE)) install.packages("shinythemes")
if (!requireNamespace("dtw", quietly = TRUE)) install.packages("dtw")
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

lapply(c("shiny", "shinythemes", "GEOquery", "limma", "dplyr", "tidyverse", 
         "ggplot2", "ggpubr", "stringr", "gridExtra", "dtw", 
         "matrixStats", "data.table"), library, character.only = TRUE)


# --- Data Preparation ---
gse <- getGEO("GSE22307", GSEMatrix = TRUE)
eSet <- gse[[1]]
raw_data <- exprs(eSet)
normalize_data <- normalizeBetweenArrays(raw_data, method = "quantile")
normalize_data1 <- log2(normalize_data + 1)

feature_data <- fData(eSet)
probe_ID <- feature_data$ID
geneSymbol <- feature_data$`Gene Symbol`
rownames(normalize_data1) <- geneSymbol[match(rownames(normalize_data1), probe_ID)]

pheno_data <- pData(eSet)
col_data <- pheno_data[, c("title", "time point (day):ch1", "geo_accession")] %>%
  rename(day = `time point (day):ch1`) %>%
  mutate(group = ifelse(day == 0, "Untreated", paste0(day, " Days of DSS treatment")),
         day_of_DSS_treatment = as.numeric(day))

# Z-score and long format
norm_z_score <- t(scale(t(normalize_data1), center = TRUE, scale = TRUE))
norm_z_score <- as.data.frame(na.omit(pmax(pmin(norm_z_score, 2), -2)))
norm_z_score$gene <- rownames(norm_z_score)

data_long <- norm_z_score %>%
  gather(key = 'samples', value = 'Normalized_expression', -gene) %>%
  left_join(., col_data, by = c("samples" = "geo_accession"))

# Inflammation genes
GIN_26 <- read.csv("dge_update_25jun24.txt", header = FALSE)
GIN26_vals <- str_to_title(GIN_26$V1)
Inflam_expression <- data_long[data_long$gene %in% GIN26_vals,]

# --- UI ---
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("Gene Expression Spline Curve Viewer with Statistics"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode", "Choose Input Mode:",
                   choices = c("Single Gene", "Upload Gene List"),
                   selected = "Single Gene"),
      conditionalPanel(
        condition = "input.mode == 'Single Gene'",
        selectInput("selected_gene", "Select Gene:", choices = unique(rownames(normalize_data1)))
      ),
      conditionalPanel(
        condition = "input.mode == 'Upload Gene List'",
        tagList(
          fileInput("gene_file", "Upload CSV file (first column = gene names)", accept = ".csv"),
          textInput("custom_label", "Custom label for test expression line:", value = "Average Expression")
        )
      ),
      downloadButton("downloadPlot", "Download Plot (PNG)"),
      downloadButton("downloadStats", "Download Stats Table (CSV)")
    ),
    mainPanel(
      plotOutput("splinePlot", width = "800px", height = "600px"),
      h4("Statistical Comparison:"),
      tableOutput("statsTable")
    )
  )
)

# --- Server ---
server <- function(input, output, session) {
  
  selected_data <- reactive({
    if (input$mode == "Single Gene") {
      gene <- input$selected_gene
      df <- data_long[data_long$gene == gene, ]
      df$label <- gene
    } else if (!is.null(input$gene_file)) {
      genes <- read.csv(input$gene_file$datapath, header = TRUE)[[1]]
      df <- data_long %>% filter(gene %in% genes)
      df <- df %>%
        group_by(samples, day_of_DSS_treatment) %>%
        summarise(Normalized_expression = mean(Normalized_expression), .groups = "drop")
      df$label <- input$custom_label
    } else {
      return(NULL)
    }
    return(df)
  })
  
  output$splinePlot <- renderPlot({
    req(selected_data())
    plot_data <- selected_data()
    
    ggplot() +
      geom_smooth(data = plot_data,
                  aes(x = day_of_DSS_treatment, y = Normalized_expression),
                  formula = y ~ s(x, bs = "cs", k = 4),
                  method = "gam",
                  color = "darkgreen", fill = "lightgreen", se = TRUE) +
      geom_smooth(data = Inflam_expression,
                  aes(x = day_of_DSS_treatment, y = Normalized_expression),
                  formula = y ~ s(x, bs = "cs", k = 4),
                  method = "gam",
                  color = "red", fill = "pink", se = TRUE) +
      theme_classic() +
      annotate("text", x = 6, y = 2.7, label = "Inflammation gene expression", color = "red", hjust = 1) +
      annotate("text", x = 6, y = 2.1, label = unique(plot_data$label), color = "darkgreen", hjust = 1)
  })
  
  stats_table <- reactive({
    req(selected_data())
    plot_data <- selected_data()
    
    test_spline <- plot_data %>%
      group_by(day_of_DSS_treatment) %>%
      summarise(expr = mean(Normalized_expression), .groups = "drop")
    
    infl_spline <- Inflam_expression %>%
      group_by(day_of_DSS_treatment) %>%
      summarise(expr = mean(Normalized_expression), .groups = "drop")
    
    merged <- merge(test_spline, infl_spline, by = "day_of_DSS_treatment", suffixes = c("_test", "_inflam"))
    
    mae <- mean(abs(merged$expr_test - merged$expr_inflam))
    global_var_test <- var(merged$expr_test)
    global_var_inflam <- var(merged$expr_inflam)
    dtw_dist <- dtw(merged$expr_test, merged$expr_inflam)$distance
    var_each_timepoint <- apply(merged[, c("expr_test", "expr_inflam")], 1, var)
    
    # Spearman correlation
    spearman_result <- cor.test(merged$expr_test, merged$expr_inflam, method = "spearman")
    spearman_rho <- round(spearman_result$estimate, 4)
    
    data.frame(
      Parameter = c("Mean Absolute Error", 
                    "Global Variance (Test)", 
                    "Global Variance (Inflammation)", 
                    "DTW Distance",
                    paste0("Variance (Day ", merged$day_of_DSS_treatment, ")"),
                    "Spearman Correlation (ρ)"),
      Value = c(round(mae, 4), 
                round(global_var_test, 4), 
                round(global_var_inflam, 4), 
                round(dtw_dist, 4),
                round(var_each_timepoint, 4),
                spearman_rho)
    )
  })
  
  output$statsTable <- renderTable({
    stats_table()
  })
  
  output$downloadStats <- downloadHandler(
    filename = function() {
      suffix <- if (input$mode == "Single Gene") input$selected_gene else "avg_genes"
      paste0("statistical_comparison_", suffix, ".csv")
    },
    content = function(file) {
      write.csv(stats_table(), file, row.names = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      suffix <- if (input$mode == "Single Gene") input$selected_gene else "avg_genes"
      paste0("gene_expression_plot_", suffix, ".png")
    },
    content = function(file) {
      req(selected_data())
      plot_data <- selected_data()
      
      p <- ggplot() +
        geom_smooth(data = plot_data,
                    aes(x = day_of_DSS_treatment, y = Normalized_expression),
                    formula = y ~ s(x, bs = "cs", k = 4),
                    method = "gam",
                    color = "darkgreen", fill = "lightgreen", se = TRUE) +
        geom_smooth(data = Inflam_expression,
                    aes(x = day_of_DSS_treatment, y = Normalized_expression),
                    formula = y ~ s(x, bs = "cs", k = 4),
                    method = "gam",
                    color = "red", fill = "pink", se = TRUE) +
        theme_classic() +
        annotate("text", x = 6, y = 2.7, label = "Inflammation gene expression", color = "red", hjust = 1) +
        annotate("text", x = 6, y = 2.1, label = unique(plot_data$label), color = "darkgreen", hjust = 1)
      
      ggsave(file, plot = p, device = "png", width = 8, height = 6, dpi = 300)
    }
  )
}

# Launch app
shinyApp(ui = ui, server = server)
