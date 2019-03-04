# This is a Shiny web application for visualization of hidden states of
# LSTM models. You can run the application by clickingthe 'Run App' button
# above if you are using RStudio
# Author: philipp.muench@helmholtz-hzi.de

library(shiny)
library(Biostrings)
library(quantmod)
library(h5)
library(keras)
library(altum)
library(pairsD3)
library(dygraphs)
library(ape)

source("utils.R")

# UI ---------------------------------------------------------------------------
ui <- fluidPage(theme = "bootstrap.css",
                
                # Application title
                titlePanel("HiddenGenome"),
                # Inputs -------------------------------------------------------
                # Sidebar with a slider input for number of bins
                sidebarLayout(
                  sidebarPanel(
                    selectInput(
                      "selectinput_dataset",
                      "Dataset",
                      c("precalculated", "generate")
                    ),
                    conditionalPanel(
                      condition = "input.selectinput_dataset == 'generate'",
                      textAreaInput(
                        "textareainput_genome",
                        "Input genome sequence (nt)",
                        width = "100%",
                        height = "200px"
                      ),
                      fileInput(
                        "fileinput_fasta",
                        "or choose file (fasta)",
                        multiple = FALSE,
                        accept = c(
                          "text/fasta",
                          "text/comma-separated-values,text/plain",
                          ".csv"
                        )
                      ),
                      fileInput(
                        "fileinput_gff3",
                        "upload gff3 annotation file",
                        multiple = FALSE,
                        accept = c(
                          "text/gff3",
                          "text/comma-separated-values,text/plain",
                          ".gff3"
                        )
                      ),
                      # let user select the model stored as .Rdata in data/model/
                      selectInput('selectinput_model', 'Select model',
                                  choice = list.files('data/model/'))
                    ),
                    conditionalPanel(
                      condition = "input.selectinput_dataset == 'precalculated'",
                      selectInput(
                        'selectinput_states',
                        'Select precalculated cell response',
                        choice = list.files('data/states/')
                      )
                    ),
                    tags$hr(),
                    selectInput(
                      "selectinput_cell",
                      "cell number(s)",
                      selected = 1,
                      choices = 1:100 ,
                      multiple = TRUE
                    )
                  ),
                  # Output -----------------------------------------------------
                  mainPanel(tabsetPanel(
                    tabPanel("position", dygraphOutput("dygraph")),
                    tabPanel(
                      "correlation",
                      pairsD3Output("pairs",
                                    width = "80%",
                                    height = 1300)
                    )
                    
                  ))
                ))
# SERVER -----------------------------------------------------------------------
server <- function(input, output, session) {
  # get the input sequence either by textareainput or fileinput
  sequence  <- reactive({
    if (isTruthy(input$textareainput_genome)) {
      input$textareainput_genome
      output <- input$textareainput_genome
      output
    } else if (isTruthy(input$fileinput_fasta)) {
      progress <- shiny::Progress$new()
      progress$set(message = "loading fasta file", value = 0)
      records <- readDNAStringSet(input$fileinput_fasta$datapath)
      first_seq <- paste(records)[1]
      on.exit(progress$close())
      first_seq
    } else {
      NULL
    }
  })
  
  # get the hidden states of the input sequence
  dataset <- reactive({
    req(input$selectinput_dataset)
    if (input$selectinput_dataset == "generate"){
      req(sequence())
      progress <- shiny::Progress$new()
      progress$set(message = "preprocess", value = 0)
      preprocessed_text <- preprocess(sequence(),
                                      vocabulary = c("\n", "a", "c", "g", "t"))
      
      #generate hdf5 file with state information
      progress$set(message = "computing states", value = 1)
      states <-
        getstates(paste0("data/model/", input$selectinput_model),
                  preprocessed_text$X,
                  type = "csv")
      on.exit(progress$close())
      states
    } else {
      req(input$selectinput_states)
      progress <- shiny::Progress$new()
      progress$set(message = "loading states", value = 1)
      states <-
        read.table(
          paste0("data/states/", input$selectinput_states),
          header = F,
          stringsAsFactors = F,
          sep = ";"
        )
      on.exit(progress$close())
      states
    }
  })
  
  # load the annotation information
  annotation <- reactive({
    if (isTruthy(input$fileinput_gff)) {
      progress <- shiny::Progress$new()
      progress$set(message = "loading gff file", value = 0)
      gff <- read.gff(input$fileinput_gff$datapath)
      # filter by direct_repeat
      gff <- gff[which(gff$type == "direct_repeat"), ]
      if (nrow(gff) > 0) {
        gff
      }
      on.exit(progress$close())
    } else {
      NULL
    }
  })
  
  output$dygraph <- renderDygraph({
    cell_num <- as.numeric(input$selectinput_cell)
    states_df <- dataset()[, cell_num]
    head(states_df)
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          states_df)
    dygraph(cell_df, main = "Cell response") %>%
      dyLegend(show = "onmouseover", hideOnMouseOut = FALSE)  %>%
      dyOptions(stepPlot = TRUE) %>% dyRangeSelector()
    # add annotation to plot
    #   if (isTruthy(annotation())) {
    #      dyShading(from = annotation()$start, to = annotation()$end, color = "#FFE6E6")
    #   }
    
  })
  
  output$pairs <- renderPairsD3({
    cell_num <- as.numeric(input$selectinput_cell)
    states_df <- dataset()[, cell_num]
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          states_df)
    is_crispr <- c(rep("A", nrow(cell_df) - 100), rep("B", 100))
    pairsD3(
      cell_df,
      group = is_crispr,
      cex = 1,
      theme = "bw",
      big = T,
      opacity = 0.5,
      leftmar = 10,
      topmar = 0,
      width = 800
    )
  })
  
  
}
# Run the application
shinyApp(ui = ui, server = server)
