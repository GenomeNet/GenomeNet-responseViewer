
# allowing reset of zooming in dygraph
dyUnzoom <- function(dygraph) {
  dyPlugin(
    dygraph = dygraph,
    name = "Unzoom",
    path = system.file("plugins/unzoom.js", package = "dygraphs")
  )
}


presAnnotation <- function(dygraph, x, text) {
  dygraph %>%
    dyAnnotation(x, text, attachAtBottom = TRUE, width = 60, height= 20)
}

vec_dyShading <- function(dy, from, to, color){
  for (i in seq_along(from)) {
    dy <- dyShading(dy, from = from[i], 
                    to = to[i]) #, 
                #    color = color[i])
  }
  dy
}


nucleotide2value <- function(nucleotides){
  require(plyr)
  nuc <- unlist(strsplit(nucleotides, "",fixed = TRUE))
  values <- mapvalues(nuc, c("A", "C", "G", "T"), c(0, 0.33, 0.66, 1))
  return(values)
}
