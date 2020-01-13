#' #' Dygraph
#' presAnnotation <- function(dygraph, x, text) {
#'   dygraph %>%
#'     dyAnnotation(x, text, attachAtBottom = TRUE, width = 60, height= 20)
#' }


#' Dygraph shading (interesting in terms of coloring interesting places in the genome)
#' @param dy
#' @param from
#' @param to
#' @param color
#' @export
vec_dyShading <- function(dy, from, to, color){
  for (i in seq_along(from)) {
    dy <- dyShading(dy, from = from[i], 
                    to = to[i])
                    #,color = color[i])
  }
  dy
}

#' Change nucleotides to values
#' 
#' A has the value 0, C 0.33, G 0.66 and T 1
#' @param nucleotides nucleotides from the genomic sequence
#' @export
nucleotide2value <- function(nucleotides){
  require(plyr)
  nuc <- unlist(strsplit(nucleotides, "", fixed = TRUE))
  values <- mapvalues(nuc, c("A", "C", "G", "T"), c(0, 0.33, 0.66, 1))
  return(values)
}

######Shiny app########
#' This is a Shiny web application for the visualization of hidden states of
#' LSTM models. Default will show the states from data: GCF_000008365.1_
#' ASM836v1_genomic and GCF_000006605.1_ASM660v1_genomic. It also comes with 
#' the possibility to use own trained models. The sample is processed before with 
#' \code{preprocessSemiRedudant()} or a Fasta file with \code{preprocessFasta()}. Then 
#' the result $x will be used for the functions \code{getStatus()} and 
#' \code{getStatusfromFasta()}.
#' 
#' @param sample  character input string of text
#' @param model.path path to trained model, default "data/models"
#' @param maxlen length of one sample
#' @param batch.size number of samples
#' @param fasta.path path to fasta files, default "data/fasta"
#' @param vocabulary used vocabulary
#' @example 
#' \dontrun{
#' visualizePrediction(strrep("ATGTAGTAGTAGTAGTAGATGATGATAGATGCACACACAGATACATAGCATGCTGCT",1000))
#' visualizePrediction(strrep("ATGTAGTAGTAGT",1000), model.path = "data/models/cpu_model.hdf5", maxlen = 50, vocabulary = c("l","a","g","c","t"))
#' visualizePrediction(strrep("ATGTAGTAGTAGT",1000))
#' visualizePrediction(paste(c(rep("ATGTAGTAGTAGTAGTAGATGATGATAGATGCACACACAGATACATAGCATGCTGCT",500), rep("TGTGTGTGTGCCCAGTACGATCAGTAGTAGGACAGAGACGAG",500)), collapse = ""))
#' visualizePrediction(fasta.path = "data/fasta/a.fasta") # all with example_model_cpu_new_full_model.hdf5}
#' @export
visualizePrediction <- function(sample = "",
                                model.path = "data/models/example_model_cpu_new_full_model.hdf5",
                                fasta.path = "",
                                maxlen = 80,
                                start_position = 800,
                                end_position = 1500,
                                batch.size = 100,
                                vocabulary = c("l","p","a","g","c","t")){

library(shiny)
library(shinythemes)
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
ui <- fluidPage(theme = shinytheme("flatly"),
                
                # Application title
                titlePanel("Response viewer for the GenomeNet", windowTitle = 'GenomeNet'),
                # Inputs -------------------------------------------------------
                # Sidebar
               # sidebarLayout(
                #  sidebarPanel(
                
                tags$p(tags$b("GenomeNet:"),
                      "For deep neural LSTM networks and their response for genomic modeling, please see our", 
                      tags$a("Wiki", href="https://hiddengenome.github.io/deepG/")),
                   

                  selectInput(
                      "selectinput_dataset",
                      "Dataset:",
                      c("Calculated Response","Examples"), width = "45%"),
                    
                    conditionalPanel(
                      condition = "input.selectinput_dataset == 'Examples'",
                      selectInput(
                        'selectinput_states',
                        'Select example for the Cell Response:',
                        choice = list.files('data/ncbi_data/states/'), width = "45%")), 
                    
                    selectInput(
                      "selectinput_cell",
                      "Cell number(s):",
                      selected = 1,
                      choices = 1:125 ,
                      multiple = TRUE, width = "45%"),

                  # Output -----------------------------------------------------
                  
                tabsetPanel(tabPanel(
                    "Cell Response:", 
                              dygraphOutput("dygraph"), 
                              dygraphOutput("dygraph2", height = 80)
                    ),
                    tabPanel("Correlation:",
                      pairsD3Output("pairs",
                                    width = "100%",
                                    height = 500)),
                    tabPanel("Correlation from the Repeat-Responses:",
                      DT::dataTableOutput("table1")),

                    tabPanel("Repeat-Response-Graph:",
                             plotOutput("plot1"))
                  ))
# SERVER -----------------------------------------------------------------------
server <- function(input, output, session) {

  sequence  <- reactive({
    if (input$selectinput_dataset == "Calculated Response") {
      if (sample == "") 
        {output <- fasta.path}
      else 
        {output <- sample}
      output
    }
    else {
      NULL
    }
  })

  # get the hidden states of the input sequence
  dataset <- reactive({
    req(input$selectinput_dataset)
    if (input$selectinput_dataset == "Calculated Response"){
      req(sequence())
      if (sample == ""){
          progress <- shiny::Progress$new()
          progress$set(message = "Computing states ...", value = 1)
          
          model <- keras::load_model_hdf5(model.path)
          states <-
            getStatesFromFasta(model,
                               fasta.path = fasta.path,
                               maxlen = maxlen,
                               batch.size = batch.size)
          on.exit(progress$close())
          states
      }
      else{
          preprocessed_text <- preprocessSemiRedundant(char = sequence(), 
                                                       maxlen = maxlen, 
                                                       vocabulary = vocabulary)
          progress <- shiny::Progress$new()
          progress$set(message = "Computing states ...", value = 1)
          states <-
            getStates(model.path = model.path,
                      preprocessed_text$X,
                      maxlen = maxlen)
          on.exit(progress$close())
        }
      states
    } else {
      req(input$selectinput_states)
      progress <- shiny::Progress$new()
      progress$set(message = "Loading states ...", value = 1)
      states <-
        read.table(
          paste0("data/ncbi_data/states/", input$selectinput_states),
          header = F,
          stringsAsFactors = F,
          sep = ";")
      on.exit(progress$close())
      states
    }
  })
  
  ##### Examples ######
  # Load coordinates of CRISPR information if input$selectinput_states
  metadata <- reactive({
    # todo: check if file exists
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
    
    dy <- dygraph(cell_df, group = "a", main = "" ) %>%
      dyLegend(show = "onmouseover", hideOnMouseOut = FALSE)  %>%
      dyOptions(stepPlot = TRUE) %>% dyRangeSelector(dateWindow = c(start_position, end_position))
    
    if (input$selectinput_dataset == "Examples"){
    dy <- dygraph(cell_df, group = "a", main = paste0(taxadata()$taxa)) %>%
      dyLegend(show = "onmouseover", hideOnMouseOut = FALSE)  %>%
      dyOptions(stepPlot = TRUE) %>% dyRangeSelector(dateWindow = c(start_position, end_position))}
    
    if (!is.null(metadata()))
      if (input$selectinput_dataset == "Examples"){
      dy <- vec_dyShading(dy, metadata()$from, metadata()$to, metadata()$color)}
      dy <- dyUnzoom(dy) 
    dy    
  })
  
  output$dygraph2 <- renderDygraph({
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          rep(0,nrow(dataset())))
 
    # load nucleotide text
    genome <- read_lines(paste0("data/ncbi_data/genome/", input$selectinput_states))
    
    ribbonData <- nucleotide2value(genome)
  
    dy <- dygraph(cell_df, group = "a") %>%
      dyRibbon(data = ribbonData, top = 1, bottom = 0, palette=c("red", "blue", "green", "yellow"))
               # A = red, c = blue, g = green, t = yellow
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
      width = 750
    )
  })
  
  output$table1 <- DT::renderDataTable({
    dat <- within_repeat_data()
    DT::datatable(dat)
  })
  
  output$plot1 <- renderPlot({
    dat <- as.data.frame(within_repeat_data())
    p <- ggplot(dat, aes(x = repeat_responses, y = non_repeat_responses, color = diff, size = diff)) + geom_point()
    print(p)}, height = 500) 

}
# Run the application
shinyApp(ui = ui, server = server)}
