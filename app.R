
# Install required packages if not already installed
if(!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if(!requireNamespace("shinythemes", quietly = TRUE)) install.packages("shinythemes")
if(!requireNamespace("dtw", quietly = TRUE)) install.packages("dtw")
if(!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if(!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if(!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
if(!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if(!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

#Load Packages
lapply(c("shiny", "shinythemes","dplyr", "tidyverse", 
         "ggplot2", "ggpubr", "stringr", "gridExtra", "dtw", 
         "matrixStats", "data.table"), library, character.only = TRUE)

#Load GSE22307 dataset
data_long <- read.csv("GSE22307_data_long.csv")


# --- Load Secondary Dataset (GSE131032_data_long.csv) ---
data_long2 <- read.csv("GSE131032_data_long.csv")

# --- Load Inflammation Genes ---
GIN_26 <- read.csv("dge_update_25jun24.txt", header = FALSE)
GIN26_vals <- str_to_title(GIN_26$V1)

# Function to get inflammation data for any dataset
get_inflam_data <- function(data) {
  data[data$gene %in% GIN26_vals, ]
}





#Function to generate spline curve

Inf_dynamics<-function(gene_list,dl){
  {test_expression<-dl[dl$gene == gene_list,]}
  ggplot() +
    geom_smooth(data = test_expression,
                aes(x = day_of_DSS_treatment, y = Zscore_expression),
                formula = y ~ s(x, bs = "cs", k = 4),
                method = "gam",
                color = "darkgreen",
                fill = "lightgreen",
                se = TRUE)+
    geom_smooth(data = get_inflam_data(g2),
                aes(x = day_of_DSS_treatment, y = Zscore_expression),
                formula = y ~ s(x, bs = "cs", k = 4),
                method = "gam",
                color = "red",
                fill = "pink",
                se = TRUE)+
    theme_classic()+
    annotate("text",x =6, y = 2.7, label = "Inflammation gene expression", color = "red", hjust = 1)+
    annotate("text",x =6, y = 2.1, label = as.character(substitute(g1)), color = "darkgreen", hjust = 1)
}

plot_spline_comparison("Mmp7",data_long2)

# Plotting function (returns ggplot object)
plot_spline_comparison <- function(gene_list, dl) {
  # Expect dl to have columns: gene, samples, day_of_DSS_treatment, Zscore_expression
  if (is.null(gene_list) || length(gene_list) == 0) stop("Provide at least one gene.")
  if (!all(c("gene","samples","day_of_DSS_treatment","Zscore_expression") %in% colnames(dl))) {
    stop("Dataset must contain columns: gene, samples, day_of_DSS_treatment, Zscore_expression")
  }
  
  # Prepare test expression (average across samples if multiple genes provided)
  plot_data <- dl %>%
    filter(gene %in% gene_list) %>%
    group_by(samples, day_of_DSS_treatment) %>%
    summarise(Zscore_expression = mean(Zscore_expression, na.rm = TRUE), .groups = "drop")
  
  infl_data <- get_inflam_data(dl)
  
  # For plot label
  label_text <- if (length(gene_list) == 1) gene_list else paste0("Avg(", length(gene_list), " genes)")
  
  ggplot() +
    geom_smooth(data = plot_data,
                aes(x = day_of_DSS_treatment, y = Zscore_expression),
                formula = y ~ s(x, bs = "cs", k = 4),
                method = "gam",
                color = "darkgreen",
                fill = "lightgreen",
                se = TRUE) +
    geom_smooth(data = infl_data,
                aes(x = day_of_DSS_treatment, y = Zscore_expression),
                formula = y ~ s(x, bs = "cs", k = 4),
                method = "gam",
                color = "red",
                fill = "pink",
                se = TRUE) +
    theme_classic() +
    theme(axis.line = element_line(color = "black"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(color="black", size=14, face="bold"),
          axis.title.x = element_text(color="black", size=12, face="bold"),
          axis.title.y = element_text(color="black", size=12, face="bold")) +
    labs(x = "Day of DSS treatment", y = "Z-score expression", title = paste("Spline:", label_text)) +
    annotate("text", x = max(dl$day_of_DSS_treatment, na.rm = TRUE), y = max(dl$Zscore_expression, na.rm = TRUE) * 0.95,
             label = "Inflammation gene expression", color = "red", hjust = 1, fontface = "bold", size = 4) +
    annotate("text", x = max(dl$day_of_DSS_treatment, na.rm = TRUE), y = max(dl$Zscore_expression, na.rm = TRUE) * 0.75,
             label = label_text, color = "darkgreen", hjust = 1, fontface = "bold", size = 4)
}


#Function to calculate spline cureve stats

calculate_expression_stats <- function(gene_list,dl) {
  

  # Get inflam expression 
  inf_genes<- get_inflam_data(dl)
  
  # Filter for selected gene(s)
  plot_data <- dl %>%
    filter(gene %in% gene_list) %>%
    group_by(samples, day_of_DSS_treatment) %>%
    summarise(Zscore_expression = mean(Zscore_expression, na.rm = TRUE), .groups = "drop")
  
  # Prepare splines (mean expression per day)
  test_spline <- plot_data %>%
    group_by(day_of_DSS_treatment) %>%
    summarise(expr = mean(Zscore_expression, na.rm = TRUE), .groups = "drop")
  
  infl_spline <- inf_genes %>%
    group_by(day_of_DSS_treatment) %>%
    summarise(expr = mean(Zscore_expression, na.rm = TRUE), .groups = "drop")
  
  # Merge for comparison
  merged <- merge(test_spline, infl_spline, by = "day_of_DSS_treatment", suffixes = c("_test", "_inflam"))
  
  # Compute statistics
  mae <- mean(abs(merged$expr_test - merged$expr_inflam))
  global_var_test <- var(merged$expr_test)
  global_var_inflam <- var(merged$expr_inflam)
  
  # Dynamic Time Warping
  dtw_dist <- dtw::dtw(merged$expr_test, merged$expr_inflam)$distance
  
  # Variance per time point
  var_each_timepoint <- apply(merged[, c("expr_test", "expr_inflam")], 1, var)
  
  # Spearman correlation
  spearman_result <- suppressWarnings(cor.test(merged$expr_test, merged$expr_inflam, method = "spearman"))
  spearman_rho <- round(spearman_result$estimate, 4)
  
  # Combine into data frame
  stats_table <- data.frame(
    Parameter = c(
      "Mean Absolute Error", 
      "Global Variance (Test)", 
      "Global Variance (Inflammation)", 
      "DTW Distance",
      paste0("Variance (Day ", merged$day_of_DSS_treatment, ")"),
      "Spearman Correlation (ρ)"
    ),
    Value = c(
      round(mae, 4), 
      round(global_var_test, 4), 
      round(global_var_inflam, 4), 
      round(dtw_dist, 4),
      round(var_each_timepoint, 4),
      spearman_rho
    )
  )
  
  return(stats_table)
}



# --- Shiny UI ---
ui <- fluidPage(
  theme = shinythemes::shinytheme("cosmo"),
  titlePanel("Gene Expression Spline Viewer (Multiple datasets)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset_choice", "Choose dataset:",
                  choices = c("GSE22307" = "data_long", "GSE131032" = "data_long2"),
                  selected = "data_long"),
      radioButtons("mode", "Input mode:", choices = c("Single Gene", "Upload Gene List"), selected = "Single Gene"),
      conditionalPanel(
        condition = "input.mode == 'Single Gene'",
        uiOutput("gene_selector_ui")
      ),
      conditionalPanel(
        condition = "input.mode == 'Upload Gene List'",
        fileInput("gene_file", "Upload CSV (first column = gene names)", accept = c(".csv")),
        textInput("custom_label", "Custom label for test expression (optional)", value = "Average Expression")
      ),
      actionButton("refresh_choices", "Refresh gene list"),
      br(), br(),
      downloadButton("downloadPlot", "Download Plot (PNG)"),
      downloadButton("downloadStats", "Download Stats (CSV)")
    ),
    mainPanel(
      plotOutput("splinePlot", width = "800px", height = "600px"),
      h4("Statistical Comparison"),
      tableOutput("statsTable")
    )
  )
)

# --- Shiny Server ---
server <- function(input, output, session) {
  # reactive dataset
  ds <- reactive({
    if (input$dataset_choice == "data_long") {
      data_long
    } else {
      data_long2
    }
  })
  
  # gene selector UI (populated from chosen dataset)
  output$gene_selector_ui <- renderUI({
    df <- ds()
    genes <- unique(df$gene)
    selectInput("selected_gene", "Select gene:", choices = genes, selected = genes[1], multiple = FALSE, selectize = TRUE)
  })
  
  # allow refresh (useful if datasets changed externally)
  observeEvent(input$refresh_choices, {
    updateSelectInput(session, "selected_gene", choices = unique(ds()$gene))
  })
  
  # reactive: list of genes to analyze
  selected_genes <- reactive({
    if (input$mode == "Single Gene") {
      req(input$selected_gene)
      return(c(input$selected_gene))
    } else {
      req(input$gene_file)
      genes <- read.csv(input$gene_file$datapath, header = TRUE, stringsAsFactors = FALSE)[[1]]
      # If user supplied a custom label we won't change gene_list, the label is used for plot title
      return(as.character(genes))
    }
  })
  
  # reactive: computed plot
  plot_obj <- reactive({
    req(selected_genes())
    dl <- ds()
    # guard for missing expected columns
    if (!"Zscore_expression" %in% colnames(dl)) {
      stop("The selected dataset must contain 'Zscore_expression' column.")
    }
    p <- plot_spline_comparison(selected_genes(), dl)
    # if custom label provided for uploaded lists, replace title
    if (input$mode == "Upload Gene List" && !is.null(input$custom_label) && input$custom_label != "") {
      p <- p + ggtitle(paste("Spline:", input$custom_label))
    }
    return(p)
  })
  
  # render plot
  output$splinePlot <- renderPlot({
    p <- plot_obj()
    print(p)
  })
  
  # reactive: stats table
  stats_tbl <- reactive({
    req(selected_genes())
    dl <- ds()
    calculate_expression_stats(selected_genes(), dl)
  })
  
  output$statsTable <- renderTable({
    stats_tbl()
  })
  
  # downloads
  output$downloadStats <- downloadHandler(
    filename = function() {
      suffix <- if (input$mode == "Single Gene") input$selected_gene else "uploaded_genes"
      paste0("stats_", input$dataset_choice, "_", suffix, ".csv")
    },
    content = function(file) {
      write.csv(stats_tbl(), file, row.names = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      suffix <- if (input$mode == "Single Gene") input$selected_gene else "uploaded_genes"
      paste0("spline_", input$dataset_choice, "_", suffix, ".png")
    },
    content = function(file) {
      ggsave(filename = file, plot = plot_obj(), device = "png", width = 8, height = 6, dpi = 300)
    }
  )
  # ---- SESSION TIMEOUT CODE ----
  # Reactive timer and reset trigger
  timeout_timer <- reactiveVal(Sys.time())
  
  # Watch input for activity (reset timeout)
  observe({
    input$selected_gene
    input$gene_file
    input$custom_label
    input$mode
    
    # Reset inactivity timer
    timeout_timer(Sys.time())
  })
  
  # Session killer if 3 minutes passed without activity
  observe({
    invalidateLater(10000, session)  # Check every 10 seconds
    if (difftime(Sys.time(), timeout_timer(), units = "secs") > 60) {
      session$close()
    }
  })
}

# --- Run app ---
shinyApp(ui = ui, server = server)

