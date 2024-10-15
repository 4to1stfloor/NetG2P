#!/usr/bin/env Rscript

packages <- c("optparse") 

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

sapply(packages, install_if_missing)

source("./Network_propagation.R")
source("./RDPN.R")
source("./GC_projection.R")

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


if (isTRUE(all.equal(colnames(expression_data), colnames(mutation_data)))) {
  print("samples are equal between expression data and mutation data")
} else {
  print("Colnames between exp and mut are not matched. NetG2P will process the rowname by interaction.")
  inter_col = intersect(colnames(expression_data), colnames(mutation_data))
  expression_data = expression_data[,inter_col]
  mutation_data = mutation_data[,inter_col]
}

if (isTRUE(all.equal(rownames(expression_data), rownames(mutation_data)))) {
  print("genes are equal between expression data and mutation data")
} else {
  print("Rownames between exp and mut are not matched. NetG2P will process the rowname by interaction.")
  inter_row = intersect(rownames(expression_data), rownames(mutation_data))
  expression_data = expression_data[inter_row,]
  mutation_data = mutation_data[inter_row,]
}

if (isTRUE(all.equal(rownames(expression_data), rownames(mutation_data))) & isTRUE(all.equal(colnames(expression_data), colnames(mutation_data)))) {
  print("pass rownames and colnames")
} else {
  stop("please, check the rownames or colnames in exp and mut data", call.=FALSE)
}


# network propagation 
result_networkpropa = Network_propagation(expression_data, mutation_data)
result_rdpn = RDPN(expression_data, mutation_data,result_networkpropa)

GC_projection(result_rdpn,  output_path = opt$output,random.num = 150)

cat(paste("Network propagation completed. Results saved to",opt$output,"\n"))

