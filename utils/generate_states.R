# this script generates hidden state responses to preprocessed text and saves it in the required format

library(deepG)
library(readr)

maxlen <- 30

model_path <- "data/model/example_run_full_model.hdf5"
# 
dir <- "ncbi_data/genome"
files <- list.files(dir, full.names = T)

for (file in files) {
  base <- basename(substr(file, 1, nchar(file) - 4))
  genome <- read_lines(file)
  genome_pre <- preprocess(genome,
                           maxlen = maxlen,
                           vocabulary = c("\n", "a", "c", "g", "t"))
  
  states <- getstates(model_path,
                      genome_pre$X,
                      maxlen = maxlen,
                      type = "csv")
  
  write.table(states, file = paste0("ncbi_data/states/" , base ,".csv"),
              col.names = F,
              row.names = F,
              quote = F, sep = ";")

}
