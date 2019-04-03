suppressPackageStartupMessages(suppressWarnings(library(optparse)))
library(tictoc)
source("retrieval_evaluation.R")

option_list = list(
  make_option(c("-t", "--test_set"), type = "character", default = NULL,
              help="test set directory ", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "output",
              help="output directory name [default= %default]", metavar="character"),
  make_option(c("-f", "--file_format"), type = "character", default = "en",
              help="(Language) format for output files [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$test_set)){
  stop("You have to provide a directory for the test set.")
}
set.seed(666)
tic()

query_evaluation(test_set = opt$test_set, outdir = opt$outdir, file_format = opt$file_format)
toc()