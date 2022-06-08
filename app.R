# Load libraries and source code

library(shiny)
library(shinyFiles)
library(fs)
library(pdftools)
library(purrr)
library(DT)
library(bslib)
library(shinyhelper)
source("generate_input_path.R")
source("iterate_graphs.R")
source("statistical_analysis.R")
source("helpers.R")


noyes <- c("no","yes")

# Build the User Interface

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "darkly"),
  
  #turn-off shown errors in final tool
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"),
  
  tags$head(
    tags$style(HTML("
    .shiny-output-error-validation {
    color: red;
    }
    "))
  ),

  tags$h1("IgIDivA", align = "center"),
  
  navlistPanel(
    tabPanel("Import Data",
             
      textInput("output_folder", "Please provide Folder name to save the final results.", value = "Results", placeholder = "Enter desired folder name here"),
             
      actionButton("create", "Create Results Folder", class = "btn-primary"),
             
             
      br(),
      br(),
      br(),
      br(),       
             
      shinyDirButton("dir", "Choose input directory", "Please select a folder", class = "btn-primary"),
      
      verbatimTextOutput("dir", placeholder = TRUE),
      
      br(),
      
      actionButton("upload", "Upload", class = "btn-primary"),
      
      br(),
      br(),
      
      selectInput("samples", "Select samples to analyze", choices = character(0), selected = character(0), multiple = TRUE, width = "100%"),
      
      textOutput("error"),
      
      br(),
      
      tags$b("Press Verify when finished with inputting samples"),
      
      br(),
      br(),
      
      actionButton("verify", "Verify", class = "btn-primary"),
      
      br(),
      br(),
      
      tags$b("Create a .txt file determining the grouping of the chosen samples for comparisons."),
      
      br(),
      
      tags$b("The file must contain two columns and the first row must contain the names of the columns."),
      
      br(),
      
      tags$b("First column contains the IDs of the samples to be analyzed. It must be named 'sample_id'."),
      
      br(),
      
      tags$b("Second column includes the name of the group that each sample belongs to. Default name is 'group_name' but it can be changed below."),
      
      br(),
      
      tags$b("The columns and their corresponding values must be separated by tabs."),
      
      br(),
      br(),
      
      textInput("groups_name", "Enter the name chosen for the second column:", value = "group_name", placeholder = "Enter column name"),
      
      br(),
      
      fileInput("samples_groups", "Upload .txt file with grouped samples:", multiple = FALSE, accept = ".txt", width = "100%", buttonLabel = "Browse...", placeholder = "No file has been uploaded")
    ),
    
    tabPanel("Set Parameters",
      
      numericInput("col_start", "Enter starting column [suggested: 5 (beginning of FR1 region), 23-59 (when using FR1 primers)]:", value = 5, step = 1),
      
      br(),
      
      numericInput("col_end", "Enter ending column [suggested: 313 (end of FR3 region)]:", value = 313, step = 1),
      
      br(),
      
      numericInput("min_reads", "Enter threshold minimum reads for the nodes [suggested: 10]:", value = 10, step = 1),
      
      br(),
      
      numericInput("p_thres", "Enter p-value threshold [suggested: 0.01 or 0.05]:", value = 0.05),
      
      br(),
      
      
      selectInput("adjust", "Do you want the p-values to be adjusted?", noyes),
      
      br(),
      
      textInput("highly_sim_clonos", "Clonotypes to be taken into account for the analysis: Choose rows from the highly_sim_clonos_file and enter their indexes in comma-separated values", value = "1", placeholder = "Enter comma-separated numbers of rows"),
      
      br(),
      
      helper(
        shiny::checkboxGroupInput("include", "Which should be included from the following?", choices = c("Summary tables", "Jumps between non-adjacent nodes", "Amino-acid mutations", "Size scaling of nodes proportional to reads", "Graph metrics", "Graph networks", "Metrics comparisons"), selected = c("Jumps between non-adjacent nodes", "Amino-acid mutations", "Size scaling of nodes proportional to reads", "Graph metrics", "Summary tables", "Graph networks", "Metrics comparisons")),
        colour = "red",
        type = "inline",
        content = "Choose which options to include in the analysis. Summary tables produces tables throughout the process. The jumps between non-adjacent nodes allows that nt vars with common SHMs differing by two or more SHMs to be included in the analysis. If the option of amino acids mutations is selected, replacement mutations will be shown in the analysis and summary tables about the replacement mutations will be produced. Other option includes the possibility of having the size of the nodes of the graph networks proportional to the number of reads of the respective nucleotide variant. Other options consist on including graph metrics, graph networks and metric comparison in the analysis. For more information, please consult the IgIDivA UserGuide."
      ),
      
      br(),
      
      uiOutput("graph_metrics_to_choose"),
      
      br(),
      
      uiOutput("comparison_metrics_to_choose"),
      
      br(),
      
      actionButton("start", "Start", class = "btn-lg btn-success"),
      
      br(),
      br(),
      
      tags$strong("Press Reset before starting a new analysis."),
      
      br(),
      
      tags$strong("Resets both input parameters and output results"),
      
      br(),
      br(),
      
      actionButton("reset", "Reset", class = "btn-lg btn-danger")
    ),
    
    tabPanel("Visualize Results",
      tabsetPanel(
        tabPanel("Summary Calculations",
          navlistPanel(
            id = "samples1"
          )
        ),
        
        tabPanel("Extra Mutations Calculations",
          navlistPanel(
           id = "samples2"
          )
        ),
        
        tabPanel("Less Mutations Calculations",
          navlistPanel(
           id = "samples3"
          )
        ),
        
        tabPanel("Mutations",
          navlistPanel(
            id = "samples4"
          )
        ),
        
        tabPanel("Amino-acid Mutations Main Variant",
          navlistPanel(
            id = "samples5"
          )
        ),
        
        tabPanel("Global Amino-acid Mutations Main Variant",
          dataTableOutput("aa_mutations_main_var_global")
        ),
        
        tabPanel("Amino-acid Mutations",
          navlistPanel(
            id = "samples6"
          )
        ),
        
        tabPanel("Global Amino-acid Mutations",
          dataTableOutput("aa_mutations_global")
        ),
        
        tabPanel("Graph Metrics",
          navlistPanel(
            id = "samples7"
          )
        ),
        
        tabPanel("Global Graph Metrics",
          dataTableOutput("graph_metrics_global")
        ),
        
        tabPanel("Graph Networks",
          navlistPanel(
            id = "samples8"
          )
        ),
        
        tabPanel("Metrics Comparisons",
          navlistPanel(
            id = "metrics"
          )
        ),
        
        tabPanel("Discarded Samples",
                 dataTableOutput("discard_samples"
          )
        )
        
      )
    )
  )
)

# Construct the Server
server <- function(input, output, session) {
  observe_helpers()
  
  volumes <- c(Home = path_home(), "R Installation" = R.home(), getVolumes()())
  
  #input Button
  shinyDirChoose(input, "dir", roots = volumes, filetypes = c("", "txt"))
  
  dir <- reactive(input$dir)
  
  output$dir <- renderText({
    parseDirPath(volumes, dir())
  })    

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$upload
    },
    handlerExpr = {
      samples <- generate_input_file(parseDirPath(volumes, dir()))
      updateSelectInput(inputId = "samples", choices = c(samples[[1]], samples[[2]]))
    }
  )
  
  files <- reactiveValues(
    ga_files = NULL,
    hsim_files = NULL,
    id = NULL,
    id_ga = NULL,
    id_hsim = NULL,
    ind_ga = NULL,
    ind_hsim = NULL
  )
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$samples
    },
    handlerExpr = {
      files$ga_files <- which(grepl('Grouped Alignment_nt_',input$samples))
      files$hsim_files <- which(grepl('highly_sim_all_clonotypes',input$samples))
      files$id <- strsplit(input$samples,'\\_|\\.')
      files$id <- unlist(lapply(files$id,function(x){return(x[length(x)-1])}))
      files$id_ga <- files$id[files$ga_files]
      files$id_hsim <- files$id[files$hsim_files]
      files$ind_ga <- order(files$id_ga)
      files$ind_hsim <- order(files$id_hsim)
      files$id_ga <- files$id_ga[files$ind_ga]
      files$ga_files <- input$samples[files$ga_files[files$ind_ga]]
      files$id_hsim <- files$id_hsim[files$ind_hsim]
      files$hsim_files <- input$samples[files$hsim_files[files$ind_hsim]]
    }
  )
  
  output$error <- renderText({
    req(input$samples)
    if(!(length(files$id_ga)==length(files$id_hsim) & all(files$id_ga==files$id_hsim))) {
      validate("Warning! There should be both higly similar and grouped alignment files for each sample.")
    }
  })
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$verify
    },
    handlerExpr = {
      req(input$samples)
      showNotification("Nice!", type = "message")
      chosen_samples <- data.table::fread(paste0(getwd(), "/input_files.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      chosen_samples <- chosen_samples[which(chosen_samples$sample_id %in% files$id_hsim), ]
      write.table(chosen_samples, paste0(getwd(),"/", input$output_folder, "/chosen_samples.txt"), sep = "\t", dec = ".",
                  row.names = FALSE, col.names = TRUE,quote = FALSE)
    }
  )
   
   observeEvent(
      ignoreNULL = TRUE,
      eventExpr = {
          input$create
      },
      handlerExpr = {
          req(input$output_folder)
          generate_output_path(input$output_folder)
          showNotification("Folder Created!", type = "message")
      }
  )
    
  
  output$graph_metrics_to_choose <- renderUI({
    if ("Graph metrics" %in% input$include) {
      helper(
        shiny::checkboxGroupInput("include_metrics", "Which graph metrics to show?", choices = c("Main variant identity", "Relative convergence (reads)", "Most relevant pathway score", "Most relevant pathway (nodes)", "End nodes density", "Max path length", "Max mutations path length", "Total reads", "Average degree", "Average distance"), selected = c("Sample ID", "Main variant identity", "Relative convergence (reads)", "Most relevant pathway score", "Most relevant pathway (nodes)", "End nodes density", "Max path length", "Max mutations path length", "Total reads", "Average degree", "Average distance")),
        colour = "red",
        type = "inline",
        content = "The main variant identity shows the percentage of identity of the main nt var with its respective germline. The relative convergence reads calculates the ratio of the number of sequences of the most relevant pathways to the number of sequences of the main nt var. The most relevant pathway score is the ratio o the total number of sequences of the nodes forming one block of pathways to the total number of sequences of all the nodes of the network with more SHMs than the main nt var. The most relevant pathway score nodes shows the number of nodes of the most relevant pathway. The end nodes density is the ratio of the number of end nodes to the number of nt vars with additional SHMs. The max path length is the number of levels of additional SHMs. The max mutations path length shows the maximum level of additional SHMs, allowing non-consecutive SHMs. The total reads shows the total number of reads of the sample. The average degree is the average total number of connections of each nt var. The average distance isthe average number of steps along the shortest pathways between each pair of nt vars."
      )
    }
  })
  
  output$comparison_metrics_to_choose <- renderUI({
    if ("Metrics comparisons" %in% input$include) {
      helper(
        shiny::checkboxGroupInput("compare_metrics", "Which metrics to compare?", choices = c("Convergence score", "End nodes density", "Max path length", "Max mutations length", "Average degree", "Average distance"), selected = c("Convergence score", "End nodes density", "Max path length", "Max mutations length", "Average degree", "Average distance")),
        colour = "red",
        type = "inline",
        content = "Select, among the graph metrics, which one(s) to use to perform comparisons between groups of samples."
      )
    }
  })
  
  clonotypes <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$highly_sim_clonos
    },
    valueExpr = {
      as.integer(unlist(strsplit(input$highly_sim_clonos, ",")))
    }
  )
  
  metrics_to_compare <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$compare_metrics
    },
    valueExpr = {
      req("Metrics comparisons" %in% input$include)
      metrics <- paste(input$compare_metrics, collapse = ",")
      metrics_to_compare <- choose_comparison_metrics(metrics)
      metrics_to_compare
    }
  )
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      withProgress(
        message = "Analysis in progress",
        detail = "This may take a while...",
        {
          for (i in 1:length(files$id_hsim)) {
            doGraph(save_path = input$output_folder, files$hsim_files[i], files$ga_files[i], files$id_hsim[i], include_jump = ("Jumps between non-adjacent nodes" %in% input$include), col_start = input$col_start, col_end = input$col_end, min_reads = input$min_reads, highly_sim_clonos = clonotypes(), nodes_size_scaling = ("Size scaling of nodes proportional to reads" %in% input$include), include_aa_muts = ("Amino-acid mutations" %in% input$include))
            incProgress(1 / length(files$id_hsim))
          }
          iterate_do_graph(argument_file=paste0(getwd(), '/', input$output_folder, "/chosen_samples.txt"), save_path = input$output_folder, include_jump = ("Jumps between non-adjacent nodes" %in% input$include), col_start = input$col_start, col_end = input$col_end, min_reads = input$min_reads, highly_sim_clonos = clonotypes(), nodes_size_scaling = ("Size scaling of nodes proportional to reads" %in% input$include), include_aa_muts = ("Amino-acid mutations" %in% input$include))
          #####################
          perform_statistical_analysis(groups_file = input$samples_groups$datapath,
                                       save_folder = input$output_folder,
                                       final_metric_table = paste0(getwd(), '/', input$output_folder, 
                                                                   "/Output/metric_table_all.txt"), 
                                       comparison_metrics = metrics_to_compare(), 
                                       compare_by = input$groups_name, 
                                       p_threshold = input$p_thres, 
                                       adjust = input$adjust)
        
          }
      )
    }
  )

  # summary calculations TAB
  summary_calculations_tables <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Summary tables" %in% input$include)
      summary_calculations_tables <- map2(files$id_hsim, input$output_folder, create_summary_calculations_table)
      summary_calculations_tables
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(summary_calculations_tables())
      iwalk(summary_calculations_tables(), ~{
        output_name <- paste0("summary_calculations_table_", .y)
        output[[output_name]] <- renderDataTable(
          {.x},
          caption = paste("Summary Calculations", files$id_hsim[.y]),
          extensions = "Buttons",
          options = list(
            paging = FALSE,
            searching = FALSE,
            dom = "Bfrtip",
            buttons = list(list(
              extend = "collection",
              buttons = c("csv", "excel", "pdf"),
              text = "Download"
            ))
          )
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(summary_calculations_tables())
      summary_calculations_tables_list <- imap(summary_calculations_tables(), ~{
        tagList(
          insertTab(
            "samples1",
            tabPanel(paste(files$id_hsim[.y]),
              dataTableOutput(
                outputId = paste0("summary_calculations_table_", .y)
              )
            )
          )
        )
      })
      tagList(summary_calculations_tables_list)
    }
  )

  # extra mutations TAB
  extra_mutations_calculations_tables <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Summary tables" %in% input$include)
      extra_mutations_calculations_tables <- map2(files$id_hsim, input$output_folder, create_extra_mutations_calculations_table)
      extra_mutations_calculations_tables
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(extra_mutations_calculations_tables())
      iwalk(extra_mutations_calculations_tables(), ~{
        output_name <- paste0("extra_mutations_calculations_table_", .y)
        output[[output_name]] <- renderDataTable(
          {.x},
          caption = paste("Extra Mutations Calculations", files$id_hsim[.y]),
          extensions = "Buttons",
          options = list(
            paging = FALSE,
            searching = FALSE,
            dom = "Bfrtip",
            buttons = list(list(
              extend = "collection",
              buttons = c("csv", "excel", "pdf"),
              text = "Download"
            ))
          )
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(extra_mutations_calculations_tables())
      extra_mutations_calculations_tables_list <- imap(extra_mutations_calculations_tables(), ~{
        tagList(
          insertTab(
            "samples2",
            tabPanel(paste(files$id_hsim[.y]),
              dataTableOutput(
                outputId = paste0("extra_mutations_calculations_table_", .y)
              )
            )
          )
        )
      })
      tagList(extra_mutations_calculations_tables_list)
    }
  )

  # less mutations TAB
  less_mutations_calculations_tables <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Summary tables" %in% input$include)
      less_mutations_calculations_tables <- map2(files$id_hsim, input$output_folder, create_less_mutations_calculations_table)
      less_mutations_calculations_tables
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(less_mutations_calculations_tables())
      iwalk(less_mutations_calculations_tables(), ~{
        output_name <- paste0("less_mutations_calculations_table_", .y)
        output[[output_name]] <- renderDataTable(
          {.x},
          caption = paste("Less Mutations Calculations", files$id_hsim[.y]),
          extensions = "Buttons",
          options = list(
            paging = FALSE,
            searching = FALSE,
            dom = "Bfrtip",
            buttons = list(list(
              extend = "collection",
              buttons = c("csv", "excel", "pdf"),
              text = "Download"
            ))
          )
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(less_mutations_calculations_tables())
      less_mutations_calculations_tables_list <- imap(less_mutations_calculations_tables(), ~{
        tagList(
          insertTab(
            "samples3",
            tabPanel(paste(files$id_hsim[.y]),
              dataTableOutput(
                outputId = paste0("less_mutations_calculations_table_", .y)
              )
            )
          )
        )
      })
      tagList(less_mutations_calculations_tables_list)
    }
  )

  ############################################################################   
  #  mutations TAB   
  mutations_tables <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Summary tables" %in% input$include)
      mutations_tables <- pmap(list(files$id_hsim, input$min_reads, input$output_folder), create_mutations_table)
      mutations_tables
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(mutations_tables())
      iwalk(mutations_tables(), ~{
        output_name <- paste0("mutations_table_", .y)
        output[[output_name]] <- renderDataTable(
          {.x},
          caption = paste("Mutations", files$id_hsim[.y]),
          extensions = "Buttons",
          options = list(
            paging = FALSE,
            searching = FALSE,
            dom = "Bfrtip",
            buttons = list(list(
              extend = "collection",
              buttons = c("csv", "excel", "pdf"),
              text = "Download"
            ))
          )
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(mutations_tables())
      mutations_tables_list <- imap(mutations_tables(), ~{
        tagList(
          insertTab(
            "samples4",
            tabPanel(paste(files$id_hsim[.y]),
              dataTableOutput(
                outputId = paste0("mutations_table_", .y)
              )
            )
          )
        )
      })
      tagList(mutations_tables_list)
    }
  )

  #  Amino-acide mutations main variant TAB 
  aa_mutations_main_var_tables <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Summary tables" %in% input$include, "Amino-acid mutations" %in% input$include)
      aa_mutations_main_var_tables <- map2(files$id_hsim, input$output_folder, create_aa_mutations_main_var_table)
      aa_mutations_main_var_tables
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(aa_mutations_main_var_tables())
      iwalk(aa_mutations_main_var_tables(), ~{
        output_name <- paste0("aa_mutations_main_var_table_", .y)
        output[[output_name]] <- renderDataTable(
          {.x},
          caption = paste("Amino - acid Mutations (main)", files$id_hsim[.y]),
          extensions = "Buttons",
          options = list(
            paging = FALSE,
            searching = FALSE,
            dom = "Bfrtip",
            buttons = list(list(
              extend = "collection",
              buttons = c("csv", "excel", "pdf"),
              text = "Download"
            ))
          )
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(aa_mutations_main_var_tables())
      aa_mutations_main_var_tables_list <- imap(aa_mutations_main_var_tables(), ~{
        tagList(
          insertTab(
            "samples5",
            tabPanel(paste(files$id_hsim[.y]),
              dataTableOutput(
                outputId = paste0("aa_mutations_main_var_table_", .y)
              )
            )
          )
        )
      })
      tagList(aa_mutations_main_var_tables_list)
    }
  )

  # Global Amino-acid Mutations Main Variant TAB
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Summary tables" %in% input$include, "Amino-acid mutations" %in% input$include)
      output$aa_mutations_main_var_global <- renderDataTable(
        {
          aa_mutations_main_var_global <- data.table::fread(paste0(getwd(), '/', input$output_folder, "/Output/aa_muts_weight_main_variant.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
          aa_mutations_main_var_global
        },
        caption = paste("Global Amino - acid Mutations (main)"),
        extensions = "Buttons",
        options = list(
          paging = FALSE,
          searching = FALSE,
          dom = "Bfrtip",
          buttons = list(list(
            extend = "collection",
            buttons = c("csv", "excel", "pdf"),
            text = "Download"
          ))
        )
      )
    }
  )

  aa_mutations_tables <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Summary tables" %in% input$include, "Amino-acid mutations" %in% input$include)
      aa_mutations_tables <- map2(files$id_hsim, input$output_folder, create_aa_mutations_table)
      aa_mutations_tables
     }
    
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(aa_mutations_tables())
      iwalk(aa_mutations_tables(), ~{
        output_name <- paste0("aa_mutations_table_", .y)
        output[[output_name]] <- renderDataTable(
          {.x},
          caption = paste("Amino - acid Mutations (rest)", files$id_hsim[.y]),
          extensions = "Buttons",
          options = list(
            paging = FALSE,
            searching = FALSE,
            dom = "Bfrtip",
            buttons = list(list(
              extend = "collection",
              buttons = c("csv", "excel", "pdf"),
              text = "Download"
            ))
          )
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(aa_mutations_tables())
      aa_mutations_tables_list <- imap(aa_mutations_tables(), ~{
        tagList(
          insertTab(
            "samples6",
            tabPanel(paste(files$id_hsim[.y]),
              dataTableOutput(
                outputId = paste0("aa_mutations_table_", .y)
              )
            )
          )
        )
      })
      tagList(aa_mutations_tables_list)
    }
  )

  #Global Amino-acid Mutations TAB
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Summary tables" %in% input$include, "Amino-acid mutations" %in% input$include)
      output$aa_mutations_global <- renderDataTable(
        {
          aa_mutations_global <- data.table::fread(paste0(getwd(), '/', input$output_folder,"/Output/aa_muts_weight.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
          aa_mutations_global
        },
        caption = paste("Global Amino - acid Mutations (rest)"),
        extensions = "Buttons",
        options = list(
          paging = FALSE,
          searching = FALSE,
          dom = "Bfrtip",
          buttons = list(list(
            extend = "collection",
            buttons = c("csv", "excel", "pdf"),
            text = "Download"
          ))
        )
      )
    }
  )
  
 #Discarded Samples TAB 
  observeEvent(
      ignoreNULL = TRUE,
      eventExpr = {
          input$start
      },
      handlerExpr = {
          req("Summary tables" %in% input$include, "Amino-acid mutations" %in% input$include)
          output$discard_samples <- renderDataTable(
              {
                  discard_samples <- data.table::fread(paste0(getwd(), '/', input$output_folder, "/Output/discarded_samples_table.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
                  discard_samples
              },
              caption = paste("Discarded Samples"),
              extensions = "Buttons",
              options = list(
                  paging = FALSE,
                  searching = FALSE,
                  dom = "Bfrtip",
                  buttons = list(list(
                      extend = "collection",
                      buttons = c("csv", "excel", "pdf"),
                      text = "Download"
                  ))
              )
          )
      }
  )
  
  
  
   # Graph Metrics TAB 

  graph_metrics_tables <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Summary tables" %in% input$include, "Graph metrics" %in% input$include)
      metrics <- paste(input$include_metrics, collapse = ",")
      graph_metrics_tables <- pmap(list(files$id_hsim, metrics, input$output_folder), create_graph_metrics_table)
      graph_metrics_tables
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(graph_metrics_tables())
      iwalk(graph_metrics_tables(), ~{
        output_name <- paste0("graph_metrics_table_", .y)
        output[[output_name]] <- renderDataTable(
          {.x},
          caption = paste("Graph Metrics", files$id_hsim[.y]),
          extensions = "Buttons",
          options = list(
            paging = FALSE,
            searching = FALSE,
            dom = "Bfrtip",
            buttons = list(list(
              extend = "collection",
              buttons = c("csv", "excel", "pdf"),
              text = "Download"
            ))
          )
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(graph_metrics_tables())
      graph_metrics_tables_list <- imap(graph_metrics_tables(), ~{
        tagList(
          insertTab(
            "samples7",
            tabPanel(paste(files$id_hsim[.y]),
              dataTableOutput(
                outputId = paste0("graph_metrics_table_", .y)
              )
            )
          )
        )
      })
      tagList(graph_metrics_tables_list)
    }
  )
 
   #Global Graph Metrics TAB 
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Summary tables" %in% input$include, "Graph metrics" %in% input$include)
      output$graph_metrics_global <- renderDataTable(
        {
         graph_metrics_global <- data.table::fread(paste0(getwd(), '/', input$output_folder,"/Output/metric_table_all.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
         graph_metrics_global
        },
        caption = paste("Global Graph Metrics"),
        extensions = "Buttons",
        options = list(
          paging = FALSE,
          searching = FALSE,
          dom = "Bfrtip",
          buttons = list(list(
            extend = "collection",
            buttons = c("csv", "excel", "pdf"),
            text = "Download"
          ))
        )
      )
    }
  )
  
  #Graph Networks TAB
  graph_networks <- eventReactive(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    valueExpr = {
      req("Graph networks" %in% input$include)
      graph_networks <- map(files$id_hsim, create_graph_network)
      graph_networks
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(graph_networks())
      iwalk(graph_networks(), ~{
        output_name <- paste0("graph_network_", .y)
        output[[output_name]] <- renderImage(
          {.x},
          deleteFile = FALSE
        )
        download_name <- paste0("download_graph_", .y)
        output[[download_name]] <- downloadHandler(
          filename = function() {
            paste0("plot", files$id_hsim[.y], ".png")
          },
          content = function(file) {
            file.copy(paste0(getwd(), "/", files$id_hsim[.y], "_ID.png"), file)
          },
          contentType = "image/png"
        )
      })
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req(graph_networks())
      graph_networks_list <- imap(graph_networks(), ~{
        tagList(
          insertTab(
            "samples8",
            tabPanel(paste(files$id_hsim[.y]),
              downloadButton(
                outputId = paste0("download_graph_", .y),
                label = "Download"
              ),
              br(),
              imageOutput(
                outputId = paste0("graph_network_", .y)
              )
            )
          )
        )
      })
      tagList(graph_networks_list)
    }
  )

  #Metrics Comparisons TAB
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Metrics comparisons" %in% input$include, "Convergence score" %in% input$compare_metrics)
      insertTab(
        "metrics",
        tabPanel("Convergence Score",
          downloadButton(
            outputId = "download_comparison_1",
            label = "Download"
          ),
          br(),
          imageOutput(
            outputId = "convergence_score_comparison"
          )
        )
      )
      output$convergence_score_comparison <- renderImage(
         {
           pdf_convert(paste0(getwd(), '/', input$output_folder, "/Comparisons/convergence_score_vs_", input$groups_name, ".pdf"), format = "png", filenames = paste0(getwd(), '/', input$output_folder, "/Comparisons/convergence_score_vs_", input$groups_name, ".png"), dpi = 300)
           list(src = paste0(getwd(), '/', input$output_folder,"/Comparisons/convergence_score_vs_", input$groups_name, ".png"), width = 800, height = 800)
         },
         deleteFile = FALSE
       )
       output$download_comparison_1 <- downloadHandler(
         filename = function() {
           paste0("convergence score vs ", input$groups_name, ".png")
         },
         content = function(file) {
           file.copy(paste0("convergence_score_vs_", input$groups_name, ".png"), file)
         },
         contentType = "image/png"
       )
     }
   )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Metrics comparisons" %in% input$include, "End nodes density" %in% input$compare_metrics)
      insertTab(
        "metrics",
        tabPanel("End Nodes Density",
          downloadButton(
            outputId = "download_comparison_2",
            label = "Download"
          ),
          br(),
          imageOutput(
            outputId = "end_nodes_density_comparison"
          )
        )
      )
      output$end_nodes_density_comparison <- renderImage(
         {
           pdf_convert(paste0(getwd(), "/", input$output_folder,"/Comparisons/end_nodes_density_vs_", input$groups_name, ".pdf"), format = "png", filenames = paste0(getwd(), "/",input$output_folder,"/Comparisons/end_nodes_density_vs_", input$groups_name, ".png"), dpi = 300)
           list(src = paste0(getwd(), "/", input$output_folder,"/Comparisons/end_nodes_density_vs_", input$groups_name, ".png"), width = 800, height = 800)
         },
         deleteFile = FALSE
       )
       output$download_comparison_2 <- downloadHandler(
         filename = function() {
           paste0("end nodes density vs ", input$groups_name, ".png")
         },
         content = function(file) {
           file.copy(paste0("end_nodes_density_vs_", input$groups_name, ".png"), file)
         },
         contentType = "image/png"
       )
     }
   )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Metrics comparisons" %in% input$include, "Max path length" %in% input$compare_metrics)
      insertTab(
        "metrics",
        tabPanel("Max Path Length",
          downloadButton(
            outputId = "download_comparison_3",
            label = "Download"
          ),
          br(),
          imageOutput(
            outputId = "max_path_length_comparison"
          )
        )
      )
     output$max_path_length_comparison <- renderImage(
         {
           pdf_convert(paste0(getwd(), '/', input$output_folder, "/Comparisons/max_path_length_vs_", input$groups_name, ".pdf"), format = "png", filenames = paste0(getwd(), '/', input$output_folder, "/Comparisons/max_path_length_vs_", input$groups_name,".png"), dpi = 300)
           list(src = paste0(getwd(), '/', input$output_folder, "/Comparisons/max_path_length_vs_", input$groups_name, ".png"), width = 800, height = 800)
         },
         deleteFile = FALSE
       )
       output$download_comparison_3 <- downloadHandler(
         filename = function() {
           paste0("max path length vs ", input$groups_name, ".png")
         },
         content = function(file) {
           file.copy(paste0("max_path_length_vs_", input$groups_name, ".png"), file)
         },
         contentType = "image/png"
       )
     }
   )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Metrics comparisons" %in% input$include, "Max mutations length" %in% input$compare_metrics)
      insertTab(
        "metrics",
        tabPanel("Max Mutations Length",
          downloadButton(
            outputId = "download_comparison_4",
            label = "Download"
          ),
          br(),
          imageOutput(
            outputId = "max_mutations_length_comparison"
          )
        )
      )
      output$max_mutations_length_comparison <- renderImage(
         {
           pdf_convert(paste0(getwd(), '/', input$output_folder,"/Comparisons/max_muts_length_vs_", input$groups_name, ".pdf"), format = "png", filenames = paste0(getwd(), '/', input$output_folder,"/Comparisons/max_muts_length_vs_", input$groups_name, ".png"), dpi = 300)
           list(src = paste0(getwd(), '/', input$output_folder,"/Comparisons/max_muts_length_vs_", input$groups_name, ".png"), width = 800, height = 800)
         },
         deleteFile = FALSE
       )
       output$download_comparison_4 <- downloadHandler(
         filename = function() {
           paste0("max mutations length vs ", input$groups_name, ".png")
         },
         content = function(file) {
           file.copy(paste0("max_muts_length_vs_", input$groups_name, ".png"), file)
         },
         contentType = "image/png"
       )
     }
   )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Metrics comparisons" %in% input$include, "Average degree" %in% input$compare_metrics)
      insertTab(
        "metrics",
        tabPanel("Average Degree",
          downloadButton(
            outputId = "download_comparison_5",
            label = "Download"
          ),
          br(),
          imageOutput(
            outputId = "average_degree_comparison"
          )
        )
      )
      output$average_degree_comparison <- renderImage(
         {
           pdf_convert(paste0(getwd(), '/', input$output_folder,"/Comparisons/avg_degree_vs_", input$groups_name, ".pdf"), format = "png", filenames = paste0(getwd(), '/', input$output_folder,"/Comparisons/avg_degree_vs_", input$groups_name, ".png"), dpi = 300)
           list(src = paste0(getwd(), '/', input$output_folder,"/Comparisons/avg_degree_vs_", input$groups_name, ".png"), width = 800, height = 800)
         },
         deleteFile = FALSE
       )
       output$download_comparison_5 <- downloadHandler(
         filename = function() {
           paste0("average degree vs ", input$groups_name, ".png")
         },
         content = function(file) {
           file.copy(paste0("avg_degree_vs_", input$groups_name, ".png"), file)
         },
        contentType = "image/png"
       )
     }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$start
    },
    handlerExpr = {
      req("Metrics comparisons" %in% input$include, "Average distance" %in% input$compare_metrics)
      insertTab(
        "metrics",
        tabPanel("Average Distance",
          downloadButton(
            outputId = "download_comparison_6",
            label = "Download"
          ),
          br(),
          imageOutput(
            outputId = "average_distance_comparison"
          )
        )
      )
      output$average_distance_comparison <- renderImage(
         {
           pdf_convert(paste0(getwd(), '/', input$output_folder,"/Comparisons/avg_distance_vs_", input$groups_name, ".pdf"), format = "png", filenames = paste0(getwd(), '/', input$output_folder,"/Comparisons/avg_distance_vs_", input$groups_name, ".png"), dpi = 300)
           list(src = paste0(getwd(), '/', input$output_folder,"/Comparisons/avg_distance_vs_", input$groups_name, ".png"), width = 800, height = 800)
         },
        deleteFile = FALSE
      )
      output$download_comparison_6 <- downloadHandler(
         filename = function() {
           paste0("average distance vs ", input$groups_name, ".png")
        },
         content = function(file) {
          file.copy(paste0("avg_distance_vs_", input$groups_name, ".png"), file)
         },
         contentType = "image/png"
       )
     }
   )

  observeEvent(
     ignoreNULL = TRUE,
     eventExpr = {
      input$reset
     },
     handlerExpr = {
      updateNumericInput(inputId = "col_start", label = "Enter starting column:", value = 5, step = 1)
      updateNumericInput(inputId = "col_end", label = "Enter ending column:", value = 313, step = 1)
       updateNumericInput(inputId = "min_reads", label = "Enter threshold minimum reads for the nodes:", value = 10, step = 1)
       updateTextInput(inputId = "highly_sim_clonos", label = "Clonotypes to be taken into account for the analysis: Choose rows from the highly_sim_clonos_file and enter their indexes in comma-separated values", value = "1", placeholder = "Enter comma-separated numbers of rows")
       updateCheckboxGroupInput(inputId = "include", label = "Which should be included from the following?", choices = c("Summary tables", "Jumps between non-adjacent nodes", "Amino-acid mutations", "Size scaling of nodes proportional to reads", "Graph metrics", "Graph networks", "Metrics comparisons"), selected = c("Jumps between non-adjacent nodes", "Amino-acid mutations", "Size scaling of nodes proportional to reads", "Graph metrics", "Summary tables", "Graph networks", "Comparison metrics", "Metrics comparisons"))
      updateCheckboxGroupInput(inputId = "include_metrics", label = "Which graph metrics to show?", choices = c("Main variant identity", "Relative convergence (reads)", "Most relevant pathway score", "Most relevant pathway (nodes)", "End nodes density", "Max path length", "Max mutations path length", "Total reads", "Average degree", "Average distance"), selected = c("Sample ID", "Main variant identity", "Relative convergence (reads)", "Most relevant pathway score", "Most relevant pathway (nodes)", "End nodes density", "Max path length", "Max mutations path length", "Total reads", "Average degree", "Average distance"))
       updateCheckboxGroupInput(inputId = "compare_metrics", label = "Which metrics to compare?", choices = c("Convergence score", "End nodes density", "Max path length", "Max mutations length", "Average degree", "Average distance"), selected = c("Convergence score", "End nodes density", "Max path length", "Max mutations length", "Average degree", "Average distance"))
       for(i in 1:8) {
         imap(files$id_hsim, ~{
           removeTab(
             paste0("samples", i),
             target = paste(files$id_hsim[.y])
           )
         })
       }
       tabs <- c("Convergence Score", "End Nodes Density", "Max Path Length", "Max Mutations Length", "Average Degree", "Average Distance")
      for (i in tabs) {
         removeTab(
          "metrics",
           target = i
         )
       }
       output$aa_mutations_main_var_global <- renderDataTable({})
       output$aa_mutations_global <- renderDataTable({})
      output$graph_metrics_global <- renderDataTable({})
       output$discard_samples <- renderDataTable({})
      }
   )
}
# Create the Shiny App object

shinyApp(ui, server)
