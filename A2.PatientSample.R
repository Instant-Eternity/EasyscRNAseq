#!/usr/bin/env Rscript
rm(list = ls())
gc()

# Load necessary libraries
library(getopt)
library(dplyr)
library(Seurat)
library(ggplot2)
library(SoupX)

# Define command line options
options <- getopt(
  command.line = TRUE,
  options = list(
    'work-path=' ,
    'data-path=',
    'sample-id=',
    'output-file=',
    'help'
  ),
longnames = c("work-path", "data-path", "sample-id", "output-file", "help"),
shortnames = c("w", "d", "s", "o", "h"),	
usage = "Usage: script.R --work-path=<work_path> --data-path=<data_path> --sample-id=<sample_id> --output-file=<output_file>"
)

# Function to print help document
print_help <- function() {
  cat("Usage: script.R\n")
  cat("-w, --work-path=<work_path>   : Path to working directory\n")
  cat("-d, --data-path=<data_path>   : Path to data files\n")
  cat("-s, --sample-id=<sample_id>   : Sample ID\n")
  cat("-o, --output-file=<output_file>: Output file\n")
  cat("-h, --help                 : Print this help message\n")
  quit(save = "no")
}

# Check if help option is provided
if ("help" %in% names(options)) {
  print_help()
}

# Check if required options are provided
if (length(options$work.path) == 0 || length(options$data.path) == 0 || length(options$sample.id) == 0 || length(options$output.file) == 0) {
  stop("Error: All parameters are required. Use -h or --help for usage information.")
}

# Access the parameters
work_path <- options$work.path
data_path <- options$data.path
sample_id <- options$sample.id
output_file <- options$output.file

# Rest of your script logic goes here
# ...


# Define command line options
options <- getopt(
	command.line = TRUE,
	"output-dir=",
	"work-dir=",
	options = list(
		longnames = c("output-dir", "work-dir"),
		shortnames = c("o", "w"),
		usage = "Usage: script.R --output-dir=<output_dir> --work-dir=<working_dir>"
	)
)

# Check if required options are provided
if (!is.null(options$output.dir) && !is.null(options$work.dir)) {
  output_dir <- options$output.dir
  setwd(options$work.dir)
} else {
  stop("Error: Both output directory (--output-dir) and working directory (--work-dir) are required.")
}

setwd(options$work.dir)

Clean_Soup <- function(path){
  out <- load10X(path) %>%
    autoEstCont() %>%
    adjustCounts()
  return (out)
}

File_path <- c('SA1', 'SA2', 'SA3', 'SA4', 'SA5', 'SA6')
Sample_info <- c('1.SLE.1', '1.SLE.2', '1.SLE.3', 
                 '2.SEP.1', '2.SEP.2', '2.SEP.3')

for(num in 1:length(File_path)){
  sample <- Clean_Soup(File_path[num]) %>%
    CreateSeuratObject(project = Sample_info[num], min.cells = 3, min.features = 200)
  saveRDS(sample, paste0("G://ZXL_221212//Data//RowData//",Sample_info[num],".rds"))
}
