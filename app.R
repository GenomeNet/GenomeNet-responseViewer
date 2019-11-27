######utlis.R - unused functions #####

#' #' Allowing reset of zooming in dygraph
#' #' @param dygraph
#' #' @export
#' dyUnzoom <- function(dygraph) {
#'   dyPlugin(
#'     dygraph = dygraph,
#'     name = "Unzoom",
#'     path = system.file("plugins/unzoom.js", package = "dygraphs")
#'   )
#' }

#' 
#' #' Dygraph
#' #' @param dygraph
#' #' @param x
#' #' @param text
#' #' @export
#' presAnnotation <- function(dygraph, x, text) {
#'   dygraph %>%
#'     dyAnnotation(x, text, attachAtBottom = TRUE, width = 60, height= 20)
#' }
#####utlis.R - used function #####
#' Dygraph
#' @param dy
#' @param from
#' @param to
#' @param color
#' @export
vec_dyShading <- function(dy, from, to, color){
  for (i in seq_along(from)) {
    dy <- dyShading(dy, from = from[i], 
                    to = to[i]) #,
                    #color = color[i])
  }
  dy
}

#' Change nucleotides to values
#' 
#' A has the value 0, C 0.33, G 0.66 and T 1
#' @param nucleotides nucleotides from the genomic sequence
#' @export
nucleotide2value <- function(nucleotides){
  require(plyr) #remove in deepG
  nuc <- unlist(strsplit(nucleotides, "", fixed = TRUE))
  values <- mapvalues(nuc, c("A", "C", "G", "T"), c(0, 0.33, 0.66, 1))
  return(values)
}

######generate.states - Unused###### 

# library(deepG)
# library(readr)
# 
# maxlen <- 30
# 
# model_path <- "data/model/example_run_full_model.hdf5"
# # 
# dir <- "ncbi_data/genome"
# files <- list.files(dir, full.names = T)
# 
# for (file in files) {
#   base <- basename(substr(file, 1, nchar(file) - 4))
#   genome <- read_lines(file)
#   genome_pre <- preprocessSemiRedundant(genome,
#                                         maxlen = maxlen,
#                                         vocabulary = c("\n", "a", "c", "g", "t"))
#   
#   states <- getStates(model_path,
#                       genome_pre$X,
#                       maxlen = maxlen,
#                       type = "csv")
#   
#   write.table(states, file = paste0("ncbi_data/states/" , base ,".csv"),
#               col.names = F,
#               row.names = F,
#               quote = F, sep = ";")
#   
# }
######Shiny app########

#' This is a Shiny web application for the visualization of hidden states of
#' LSTM models. Default will show the states from data: GCF_000008365.1_
#' ASM836v1_genomic and GCF_000006605.1_ASM660v1_genomic. It also comes with 
#' the possibility to use trained models.
#' 
#' @param sample  preprocess$x -> getStates
#' @param model.path path to model
#' @param maxlen maxlen 
#' @param batch.size batch.size
#' @param fasta.path path to fasta files -> preprocessFasta -> getStatesfromFasta
#' @export
visualizePrediction <- function(sample = "",
                                model.path = "data/model/example_run_full_model.hdf5",
                                maxlen = 30,
                                batch.size = 100,
                                fasta.path = "example_files/fasta"){

library(shiny) #remove in deepG
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

# Generate the user interface
# UI ---------------------------------------------------------------------------
ui <- fluidPage(theme = "www/bootstrap.css",
                
                # Application title
                titlePanel("GenomeNet - deepG"),
                # Inputs -------------------------------------------------------
                # Sidebar
                sidebarLayout(
                  sidebarPanel(
                    selectInput(
                      "selectinput_dataset",
                      "Dataset:",
                      c("calculated hidden states","examples")
                    ),
                    
                    conditionalPanel(
                      condition = "input.selectinput_dataset == 'calculated hidden states'",
                      selectInput('selectinput_model', 'Select model:',
                                  choice = c(list.files('data/model/'), model.path))), #Add the models right 
                    
                    conditionalPanel(
                      condition = "input.selectinput_dataset == 'examples'",
                      selectInput(
                        'selectinput_states',
                        'Select precalculated cell response:',
                        choice = list.files('data/ncbi_data/states/'))),
                    
                    tags$hr(),
                    selectInput(
                      "selectinput_cell",
                      "Cell number(s):",
                      selected = 1,
                      choices = 1:125 ,
                      multiple = TRUE
                    ),
                    numericInput(
                      "numericInput_start",
                      "Start position:",
                      800,
                      min = 1,
                      max = 4000),
                    numericInput(
                      "numericInput_end",
                      "End position:",
                      1500,
                      min = 1,
                      max = 4000)
                  ),
                  # Output -----------------------------------------------------
                  mainPanel(tabsetPanel(
                    tabPanel("Position:", 
                             dygraphOutput("dygraph"), 
                             dygraphOutput("dygraph2", height = 50)),
                    tabPanel("Correlation:",
                      pairsD3Output("pairs",
                                    width = "80%",
                                    height = 1200)),
                    tabPanel("Repeat-Responses:",
                      DT::dataTableOutput("table1")),
                    
                    tabPanel("Repeat-Responses-Graph:",
                             plotOutput("plot1"))
                  ))
                ))
# SERVER -----------------------------------------------------------------------
server <- function(input, output, session) {
  
  sequence  <- reactive({
    if (input$selectinput_dataset == "calculated hidden states" &&
        isTruthy(input$textareainput_genome)) {
      input$textareainput_genome
      output <- input$textareainput_genome
      output
    } else if (isTruthy(input$fileinput_fasta)) {
      progress <- shiny::Progress$new()
      progress$set(message = "Loading fasta file", value = 0)
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
    if (input$selectinput_dataset == "calculated hidden states"){
      req(sequence())
      progress <- shiny::Progress$new()
      progress$set(message = "Preprocessing ...", value = 0)
      preprocessed_text <- preprocessSemiRedundant(char = sequence(),
                                      vocabulary = c("\n", "a", "c", "g", "t"), maxlen = maxlen)

      #generate hdf5 file with state information
      progress$set(message = "Computing states ...", value = 1)

        states <-
        getStates(paste0("data/model/", input$selectinput_model),
                  preprocessed_text$X,
                  type = "csv")
      on.exit(progress$close())
    } else {
      req(input$selectinput_states)
      progress <- shiny::Progress$new()
      progress$set(message = "Loading states ...", value = 1)
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
 ##### Add: the sane with the data from fasta.files #####
  
  
  ##### Examples ######
  # Load coordinates of CRISPR information if input$selectinput_states
  metadata <- reactive({
    # todo: check if file exists
    print("Load metadata ...")
    req(input$selectinput_states)
    progress <- shiny::Progress$new()
    progress$set(message = "Loading array annotation ...", value = 1)
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
  
  # Load coordinates of CRISPR information if input$selectinput_states
  taxadata <- reactive({
    print("Load taxonomic data ...")
    req(input$selectinput_states)
    progress <- shiny::Progress$new()
    progress$set(message = "Loading taxon annotation ...", value = 1)
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
      progress$set(message = "Loading gff file ...", value = 0)
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
  
  #### Plotting #####
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

}
# Run the application
shinyApp(ui = ui, server = server)}
