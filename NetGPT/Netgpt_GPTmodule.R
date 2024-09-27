#!/usr/bin/env Rscript

require(optparse)

source("network_propagation.R")
source("RDPN.R")
source("GC.R")
source("GC_projection.R")

# option
option_list <- list(
  make_option(c("--expression"), type="character", default=NULL, 
              help="Path to the expression data CSV file", metavar="character"),
  make_option(c("--mutation"), type="character", default=NULL, 
              help="Path to the mutation data CSV file", metavar="character"),
  make_option(c("--output"), type="character", default= "final_res/", 
              help="Path to final result", metavar="character")
)

# option parsing
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# confirm the input
if (is.null(opt$expression) || is.null(opt$mutation)) {
  print_help(opt_parser)
  stop("Both --expression and --mutation files must be provided.", call.=FALSE)
}

expression_data <- read.csv(opt$expression, row.names = 1)
mutation_data <- read.csv(opt$mutation, row.names = 1)

if (all.equal(colnames(expression_data), colnames(mutation_data))) {
  print("samples are equal between expression data and mutation data")
} else {
  stop("please, make exp and mut data have the same columns", call.=FALSE)
}

if (all.equal(rownames(expression_data), rownames(mutation_data))) {
  print("genes are equal between expression data and mutation data")
} else {
  stop("please, make exp and mut data have the same rows", call.=FALSE)
}

# network propagation 
result_networkpropa = network_propagation(expression_data, mutation_data)
result_rdpn = RNPN(expression_data, mutation_data,network_propagation)

GC_projection(result_rdpn,  output_path = opt$output,random.num = 150)

cat(paste("Network propagation completed. Results saved to",opt$output,"\n"))

