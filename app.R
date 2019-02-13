# This is a Shiny web application for visualization of hidden states of 
# LSTM models. You can run the application by clickingthe 'Run App' button 
# above if you are using RStudio
# Author: philipp.muench@helmholtz-hzi.de

library(shiny)
library(Biostrings)
library(quantmod)
library(h5)
library(keras)
library(HiddenGenome)
library(dygraphs)

source("utils.R")

ui <- fluidPage(theme = "bootstrap.css",
                
                # Application title
                titlePanel("HiddenGenome"),
                
                # Sidebar with a slider input for number of bins
                sidebarLayout(
                  sidebarPanel(
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
                    tags$hr(),
                    
                    # let user select the model stored as .Rdata in data/model/
                    selectInput('selectinput_model', 'Select model',
                                choice = list.files('data/model/')),
                    
                    selectInput(
                      "selectinput_cell",
                      "cell number(s)",
                      selected = 1,
                      choices = 1:100 ,
                      multiple = TRUE
                    )
                  ),
                  
                  # Show a plot of the generated distribution
                  mainPanel(dygraphOutput("dygraph")
                  )))

# Define server logic required to draw a histogram
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
    req(sequence())
    progress <- shiny::Progress$new()
    progress$set(message = "preprocess", value = 0)
    preprocessed_text <- preprocess(sequence(),
                                    vocabulary = c("\n", "a", "c", "g", "t"))
    
    #generate hdf5 file with state information
    progress$set(message = "computing states", value = 1)
    getstates(paste0("data/model/", input$selectinput_model),
              preprocessed_text$X,
              type = "csv")
    # load hdf5
    states <- read.table(
      "output_states.csv",
      header = F,
      stringsAsFactors = F,
      sep = ";"
    )
    on.exit(progress$close())
    states
  })
  
  output$dygraph <- renderDygraph({
    cell_num <- as.numeric(input$selectinput_cell)
    states_df <- dataset()[, cell_num]
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          states_df)
    dygraph(cell_df, main = "Cell response") %>%
      dyLegend(show = "onmouseover", hideOnMouseOut = FALSE)  %>%
      dyOptions(stepPlot = TRUE) %>% dyRangeSelector()
  })
}
# Run the application
shinyApp(ui = ui, server = server)

