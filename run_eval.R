source("ngram_metrics.R")
source("similarity.R")

option_list = list(
  make_option(c("-t", "--test_set"), type = "character", default = NULL,
              help="test set directory ", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "output",
              help="output directory name [default= %default]", metavar="character"),
  make_option(c("-m", "--max_n"), type = "integer", default = 10,
              help="max. N-gram length [default= %default]", metavar = "character"),
  make_option(c("-f", "--file_format"), type = "character", default = "en",
              help="(Language) format for output files [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$test_set)){
  stop("You have to provide a directory for the test set.")
}
set.seed(666)
do_all(test_set = opt$test_set, max_n =  opt$max_n, outdir = opt$outdir, file_format = opt$file_format)
