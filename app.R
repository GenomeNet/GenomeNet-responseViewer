# This is a Shiny web application for visualization of hidden states of
# LSTM models. You can run the application by clickingthe 'Run App' button
# above if you are using RStudio
# Author: philipp.muench@helmholtz-hzi.de

library(shiny)
library(Biostrings)
library(readr)
library(quantmod)
library(h5)
library(keras)
library(deepG)
library(DT)
library(pairsD3)
library(ggplot2)
library(dygraphs)
library(ape)
library(ggvis)

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
                        choice = list.files('data/ncbi_data/states/')
                      )
                    ),
                    tags$hr(),
                    selectInput(
                      "selectinput_cell",
                      "cell number(s)",
                      selected = 1,
                      choices = 1:125 ,
                      multiple = TRUE
                    ),
                    numericInput(
                      "numericInput_start",
                      "start position",
                      800,
                      min = 1,
                      max = 4000),
                    numericInput(
                      "numericInput_end",
                      "end position",
                      1500,
                      min = 1,
                      max = 4000)
                  ),
                  # Output -----------------------------------------------------
                  mainPanel(tabsetPanel(
                    tabPanel("position", dygraphOutput("dygraph"), dygraphOutput("dygraph2", height = 50)),
                    tabPanel(
                      "correlation",
                      pairsD3Output("pairs",
                                    width = "80%",
                                    height = 1300)
                    ),
                    tabPanel(
                      "repeat_correlation",
                      DT::dataTableOutput("table1"),
                      plotOutput("plot1")
                #      ggvisOutput("plot2")
                    )
                  ))
                ))
# SERVER -----------------------------------------------------------------------
server <- function(input, output, session) {
  # get the input sequence either by textareainput or fileinput
  sequence  <- reactive({
    if (input$selectinput_dataset == "generate" && 
        isTruthy(input$textareainput_genome)) {
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
          paste0("data/ncbi_data/states/", input$selectinput_states),
          header = F,
          stringsAsFactors = F,
          sep = ";"
        )
      on.exit(progress$close())
      states
    }
  })
  
  # load coordinates of CRISPR information if input$selectinput_states
  metadata <- reactive({
    # todo: check if file exists
    print("load metadata")
    req(input$selectinput_states)
    progress <- shiny::Progress$new()
    progress$set(message = "loading array annotation", value = 1)
    metadata <-
      read.table(
        paste0("data/ncbi_data/position/", input$selectinput_states),
        header = F,
        stringsAsFactors = F,
        sep = ";"
      )
    on.exit(progress$close())
    colnames(metadata) <- c("from", "to")
    metadata
  })
  
  # load coordinates of CRISPR information if input$selectinput_states
  taxadata <- reactive({
    # todo: check if file exists
    print("load taxonomic data")
    req(input$selectinput_states)
    progress <- shiny::Progress$new()
    progress$set(message = "loading taxon annotation", value = 1)
    taxadata <-
      read.table(
        paste0("data/ncbi_data/meta/", input$selectinput_states),
        header = T,
        stringsAsFactors = F,
        sep = "\t"
      )
    on.exit(progress$close())
    taxadata
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
  
  within_repeat_data <- reactive({
    if (!is.null(dataset())) {
      responses <- dataset()
      datalist = list()
      for (i in 1:nrow(metadata())) {
        datalist[[i]] <- metadata()[i,1]:metadata()[i,2]
      }
      repeat_posistions <- do.call(c, datalist)
      repeat_responses <- colMeans(responses[repeat_posistions,])
      non_repeat_responses <- colMeans(responses[-repeat_posistions,])
      diff <- repeat_responses - non_repeat_responses
      data <- rbind(repeat_responses, non_repeat_responses, diff)
      as.data.frame(t(data))
    }
  })
  
  #### plotting
  output$dygraph <- renderDygraph({
    cell_num <- as.numeric(input$selectinput_cell)
    states_df <- dataset()[, cell_num]
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          states_df)
    dy <- dygraph(cell_df, group = "a", main = paste0(taxadata()$taxa)) %>%
      dyLegend(show = "onmouseover", hideOnMouseOut = FALSE)  %>%
      dyOptions(stepPlot = TRUE) %>% dyRangeSelector(dateWindow = c(input$numericInput_start, input$numericInput_end))
    if (!is.null(metadata()))
      dy <- vec_dyShading(dy, metadata()$from, metadata()$to, metadata()$color)
    dy    
  })
  
  output$dygraph2 <- renderDygraph({
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          rep(0,nrow(dataset())))
 
    # load nucleotide text
    genome <- read_lines(paste0("data/ncbi_data/genome/", input$selectinput_states))
    ribbonData <- nucleotide2value(genome)
  
    dy <- dygraph(cell_df, group = "a") %>%
      dyRibbon(data = ribbonData, top = 1, bottom = 0, palette=c("red", "blue", "green", "orange"))
    if (!is.null(metadata()))
      dy <- vec_dyShading(dy, metadata()$from, metadata()$to, metadata()$color)
    dy    
  })
  
  output$pairs <- renderPairsD3({
    cell_num <- as.numeric(input$selectinput_cell)
    states_df <- dataset()[, cell_num]
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          states_df)
    
    # get the positions of repeat
    datalist = list()
    for (i in 1:nrow(metadata())) {
      datalist[[i]] <- metadata()[i,1]:metadata()[i,2]
    }
    repeat_posistions <- do.call(c, datalist)
    is_crispr <- rep("A", nrow(cell_df))
    is_crispr[repeat_posistions] <- "B"
  
    pairsD3(
      cell_df,
      group = is_crispr,
      cex = 1,
 #     theme = "bw",
      big = T,
      opacity = 0.5,
      leftmar = 10,
      topmar = 0,
      width = 800
    )
  })
  
  output$table1 <- DT::renderDataTable({
    dat <- within_repeat_data()
    DT::datatable(dat)
  })

  
  output$plot1 <- renderPlot({
    dat <- as.data.frame(within_repeat_data())
    p <- ggplot(dat, aes(x = repeat_responses, y = non_repeat_responses, size=diff)) + geom_point()
    print(p)
  }, height = 700)

  
#  vis <- reactive({
#    if(!is.null(within_repeat_data())){
#      dat <- as.data.frame(within_repeat_data())
#      dat 
#    }
#    })
 # within_repeat_data %>% bind_shiny("plot2")
  
  
}
# Run the application
shinyApp(ui = ui, server = server)
