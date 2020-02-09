#' Dygraph shading to color interesting places in the genome
#' 
#' @param dy the dygraph
#' @param from begin from the shading
#' @param to end from the shading
#' @param colour color 
#' @export
vec_dyShading <- function(dy, from, to, colour){
  for (i in seq_along(from)) {
    dy <- dyShading(dy, from = from[i], 
                    to = to[i],
                    color = colour) #for different colours colors[i]
  }
  dy
}

#' Change nucleotides to values
#' 
#' @param nucleotides nucleotides from the genomic sequence
#' @export
nucleotide2value <- function(nucleotides){
  require(plyr)
  nuc <- unlist(strsplit(nucleotides, "", fixed = TRUE))
  values <- mapvalues(nuc, c("A", "C", "G", "T"), c(0, 0.33, 0.66, 1))
  return(values)
}

# Shiny app --------------------------------------------------------------------

#' This is a Shiny web application for the visualization of hidden states of
#' LSTM models from the package deepG and additionally there are four examples. It also comes with 
#' the possibility to use own trained models.
#' 
#' @param sample character input string of text
#' @param states if states = TRUE, then load states from states_path
#' @param model.path path to trained model
#' @param fasta.path path to fasta files
#' @param states.path path to the .h5 file with the calculated states
#' @param start_position start from the dygraph
#' @param end_position end from the dygraph
#' @param batch.size number of samples to evaluate at once
#' @param vocabulary used vocabulary
#' @param show_correlation when false it does not show the information about the correlation 
#' @param layer_depth depth of layer to evaluate
#' @param cell_number cell_number showed in the dygraph
#' @param step frequency of sampling steps
#' @param padding TRUE/FALSE generate states for the first maxlen nucleotides
#' @export
visualizePrediction <- function(sample = "",
                                states = FALSE,
                                model.path = "",
                                fasta.path = "",
                                states.path = "",
                                start_position = 0,
                                end_position = 800,
                                batch.size = 200,
                                vocabulary = c("l","a","g","c","t"),
                                show_correlation = TRUE,
                                layer_depth = 1,
                                cell_number = 1,
                                step = 1,
                                padding = TRUE){

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

# UI (User inferace) -----------------------------------------------------------
ui <- fluidPage(titlePanel("Response Viewer for deepG", windowTitle = 'GenomeNet'),
                
                # Inputs -------------------------------------------------------
                tags$p(tags$b("GenomeNet:"),
                    "For training deep neural LSTM networks for genomic modeling and visualising their hidden states,
                    please see our", tags$a("Wiki", href="https://github.com/hiddengenome/deepG/wiki"), "for help."),
                   
                selectInput(
                    "selectinput_dataset",
                    "Dataset:",
                    c("Calculated Response","Examples"), width = "45%"),
                    
                conditionalPanel(
                    condition = "input.selectinput_dataset == 'Examples'",
                    selectInput(
                        'selectinput_states',
                        'Select example for the Cell Response:',
                        choice = list.files('data/ncbi_data/genome/'), width = "45%")), 
                    
                selectInput(
                    "selectinput_cell",
                    "Cell number(s):",
                    selected = cell_number,
                    choices = 1:125 ,
                    multiple = TRUE, 
                    width = "45%"),

                # Output -------------------------------------------------------
                tabsetPanel(tabPanel(
                    "Cell Response:", 
                    dygraphOutput("dygraph"), 
                    dygraphOutput("dygraph2", height = 80)),
                    
                tabPanel("Correlation:",
                    pairsD3Output("pairs",
                    width = "100%",
                    height = 500)),
                
                tabPanel("Correlation from the Repeat-Responses:",
                    DT::dataTableOutput("table1")),

                tabPanel("Repeat-Response-Graph:",
                    plotOutput("plot1"))))
# SERVER -----------------------------------------------------------------------
server <- function(input, output, session) {

  # load maxlen from the model
  
  if(is.null(states.path)){
  model <- keras::load_model_hdf5(model.path) 
  maxlen <- model$input$shape[1] }

  # get the hidden states of the input sequence
  dataset <- reactive({
    req(input$selectinput_dataset)
    if (input$selectinput_dataset == "Calculated Response"){
      if (states == TRUE){
        progress <- shiny::Progress$new()
        progress$set(message = "Load generated states ...", value = 1)
        states <- readRowsFromH5(h5_path = states.path, complete = TRUE)
        on.exit(progress$close())
        states}
      else{
       if (sample == ""){
        progress <- shiny::Progress$new()
        progress$set(message = "Preprocessing and computing states from Fasta file...", value = 1)
        writeStatesByFastaEntries(fasta.path = fasta.path, model.path = model.path, step = step, layer_depth = layer_depth,
                                  vocabulary = vocabulary, batch.size = batch.size, filename = "fasta_states", 
                                  padding = padding, h5_path = "data")
        states <- readRowsFromH5(h5_path = "data/Nr1fasta_states.h5", complete = TRUE)
        on.exit(progress$close())}
      else{
        progress <- shiny::Progress$new()
        progress$set(message = "Preprocessing and computing states ...", value = 1)
        states <- writeStatesToH5(model.path = model.path, sequence = sample, batch.size = batch.size, 
                                  layer_depth = layer_depth, h5_filename = "states", vocabulary = vocabulary, 
                                  step = step, returnStates = TRUE, padding = padding)
        on.exit(progress$close())}
      states}}
    else {
      req(input$selectinput_states)
      progress <- shiny::Progress$new()
      progress$set(message = "Loading states from the examples ...", value = 1)
      states <-
        read.table(
          paste0("data/ncbi_data/states/", input$selectinput_states),
          header = F,
          stringsAsFactors = F,
          sep = ";")
      on.exit(progress$close())
      states}
  })
  
  # Load coordinates of CRISPR information from the examples
  metadata <- reactive({
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
  
  # Create data with the repeat responses for the table and the plot 
  within_repeat_data <- reactive({
    if (!is.null(dataset())) {
      responses <- dataset()
      datalist = list()
      for (i in 1:nrow(metadata())) {
        datalist[[i]] <- metadata()[i,1]:metadata()[i,2]}
      
      repeat_posistions <- do.call(c, datalist)
      repeat_responses <- colMeans(responses[repeat_posistions,])
      non_repeat_responses <- colMeans(responses[-repeat_posistions,])
      diff <- repeat_responses - non_repeat_responses
      data <- rbind(repeat_responses, non_repeat_responses, diff)
      as.data.frame(t(data))
    }
  })
  
  # Plotting -------------------------------------------------------------------
  output$dygraph <- renderDygraph({
    cell_num <- as.numeric(input$selectinput_cell)
    states_df <- dataset()[, cell_num]
    cell_df <- data.frame(pos = 1:nrow(dataset()), states_df)
    
    dy <- dygraph(cell_df, group = "a", main = "" ) %>%
      dyLegend(show = "onmouseover", hideOnMouseOut = FALSE)  %>%
      dyOptions(stepPlot = TRUE) %>% dyRangeSelector(dateWindow = c(start_position, end_position))

    if (!is.null(metadata()))
      if (input$selectinput_dataset == "Examples"){
      dy <- vec_dyShading(dy, metadata()$from, metadata()$to, "lightgrey")} 
      dy <- dyUnzoom(dy) 
    dy    
  })
  
  output$dygraph2 <- renderDygraph({
    cell_df <- data.frame(pos = 1:nrow(dataset()),
                          rep(0,nrow(dataset())))
    
    if (input$selectinput_dataset == "Examples"){
      genome <- read_lines(paste0("data/ncbi_data/genome/", input$selectinput_states))
      ribbonData <- nucleotide2value(genome)} # nucleotides into values [A = 0; C = 0.33; G = 0.66: T = 1]
    
    else { 
      genome <- sample
      
      if(states == FALSE){
        if (sample == ""){
          genome <- Biostrings::readDNAStringSet(fasta.path)
          genome <- paste0(paste(genome, collapse = "p"))}}
        else{print("No data about sequence, shows only the calculated states")
          genome <- ""}
      
      ribbonData <- nucleotide2value(genome)} # nucleotides into values [A = 0; C = 0.33; G = 0.66: T = 1]}
    
    dy <- dygraph(cell_df, group = "a") %>%
      dyRibbon(data = ribbonData, top = 1, bottom = 0, palette=c("yellow", "red", "green", "blue"))
      # [A = yellow; C = red; G = green: T = blue]
    dy    
  })
  
  
  if (show_correlation == TRUE){
  output$pairs <- renderPairsD3({
    cell_num <- as.numeric(input$selectinput_cell)
    states_df <- dataset()[, cell_num]
    cell_df <- data.frame(pos = 1:nrow(dataset()), states_df)
    
    # get the positions of repeat
    datalist = list()
    for (i in 1:nrow(metadata())) {
      datalist[[i]] <- metadata()[i,1]:metadata()[i,2]}
    
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
    p <- ggplot(dat, aes(x = repeat_responses, y = non_repeat_responses, color = diff)) + geom_point() +
      scale_color_gradient(low = "blue", high = "red")
    print(p)}, height = 500) 
  }} 

# Run the application ----------------------------------------------------------
shinyApp(ui = ui, server = server)}