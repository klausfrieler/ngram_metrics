source("./ngram_evaluation.R")
suppressPackageStartupMessages(suppressWarnings(library(optparse)))


write_stats <- function(stats_df, outdir = ".", name, file_format = "en"){
  sep <- ","
  dec <- "."
  if(file_format == "de"){
    sep <- ";"
    dec <- ","
  }
  write.table(stats_df,
              file.path(outdir, name),
              sep = sep,
              dec = dec,
              row.names = F,
              col.names = T,
              quote = F)

}
do_all <- function(test_set_dir, max_n = 10, outdir = ".", file_format = "en", recalc = T){
  if(recalc){
    messagef("Setting up workspace...")
    setup_workspace(test_set_dir)
    messagef("...done.")
    #rm(list = sprintf("ngram_stats"), globalenv())
    #rm(list = sprintf("ngram_analysis"), globalenv())

    messagef("Prepraring n-grams...")
    prepare_ngrams(max_n = max_n)
    messagef("...done.")
  }

  messagef("Analyzing n-grams...")
  ngram_analysis <- get_ngram_analysis(max_n = max_n, thresholds = seq(0.03, 0.07, 0.02))
  assign("ngram_analysis", ngram_analysis, globalenv())
  messagef("...done.")

  messagef("Calculating stats...")
  ngram_stats <- get_ngram_stats_retrieval(ngram_analysis)
  assign("ngram_stats", ngram_stats, globalenv())
  messagef("...done.")

  messagef("Writing results to %s...", outdir)
  if(!dir.exists(outdir)){
    messagef("Creating output directory: %s", outdir)
    dir.create(outdir)
  }
  write_stats(ngram_stats[[1]], outdir = outdir, name = "ngram_stats_raw.csv", file_format = file_format)
  write_stats(ngram_stats[[2]], outdir = outdir, name = "ngram_stats_solo.csv", file_format = file_format)
  write_stats(ngram_stats[[3]], outdir = outdir, name = "ngram_stats_sum.csv", file_format = file_format)
  messagef("...done.")

  messagef("Producing figures...")
  q <- ngram_stats[[2]] %>% ggplot(aes(x = factor(target_n), y = F1, fill = factor(threshold)))
  q <- q + geom_boxplot() + geom_point(alpha = .02) + geom_violin(alpha = .1)
  #q <- q + facet_wrap(~threshold, ncol = 1)
  q <- q + theme_bw() + labs(x = "N")

  ggsave(file.path(outdir, "ngram_stats_solo.png"), width = 9.03, heigh = 5.73, units = "in")
  messagef("...done.")
  ngrs <- ngram_stats[[3]] %>% ungroup() %>% rename(n = target_n, thresh = threshold)
  vars <- c("thresh", "n", sort(names(ngrs)[grep("F1_mean", names(ngrs))]))
  print(ngrs %>% select(vars) %>% arrange(n, desc(F1_mean)))
  vars <- c("thresh", "n", sort(names(ngrs)[grep("prec", names(ngrs))]))
  print(ngrs %>% select(vars) %>% arrange(n, desc(prec_mean)))
  vars <- c("thresh", "n", sort(names(ngrs)[grep("^rec", names(ngrs))]))
  print(ngrs %>% select(vars) %>% arrange(n, desc(rec_mean)))
  #q
}
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
do_all(test_set = opt$test_set, max_n =  opt$max_n, outdir = opt$outdir, file_format = opt$file_format)