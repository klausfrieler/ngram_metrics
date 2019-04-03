source("./ngram_evaluation.R")
source("./similarity.R")
source("./utils.R")

library(tictoc)
do_all <- function(test_set_dir, max_n = 10, outdir = ".", file_format = "en", recalc = T, thresholds = seq(0.03, 0.07, 0.02)){
  #tic()
  if(recalc){
    #tic()
    messagef("Setting up workspace...")
    setup_workspace(test_set_dir)
    messagef("...done.")
    #toc()
    #rm(list = sprintf("ngram_stats"), globalenv())
    #rm(list = sprintf("ngram_analysis"), globalenv())
    #tic()
    messagef("Prepraring n-grams...")
    prepare_ngrams(max_n = max_n)
    messagef("...done.")
    #toc()
  }

  #tic()
  messagef("Analyzing n-grams...")
  ngram_analysis <- get_ngram_analysis(max_n = max_n, thresholds = thresholds)
  assign("ngram_analysis", ngram_analysis, globalenv())
  messagef("...done.")
  #toc()

  #tic()
  messagef("Calculating stats...")
  ngram_stats <- get_ngram_stats_retrieval(ngram_analysis)
  assign("ngram_stats", ngram_stats, globalenv())
  messagef("...done.")
  #toc()

  #tic()
  messagef("Writing results to %s...", outdir)
  if(!dir.exists(outdir)){
    messagef("Creating output directory: %s", outdir)
    dir.create(outdir)
  }
  write_stats(ngram_stats[[1]], outdir = outdir, name = "ngram_stats_raw.csv", file_format = file_format)
  write_stats(ngram_stats[[2]], outdir = outdir, name = "ngram_stats_solo.csv", file_format = file_format)
  write_stats(ngram_stats[[3]], outdir = outdir, name = "ngram_stats_sum.csv", file_format = file_format)
  messagef("...done.")
  #toc()

  #tic()
  messagef("Producing figures...")
  #def_colour <- "#1f77b4"
  def_colour <- "black"

  q <- ngram_stats[[2]] %>% ggplot(aes(x = factor(target_n), y = F1, fill = factor(threshold)))
  q <- q + geom_boxplot(fill = def_colour) + geom_point(alpha = .1, color = "black") + geom_violin(alpha = .1, fill = def_colour)
  #q <- q + facet_wrap(~threshold, ncol = 1)
  q <- q + theme_bw() + labs(x = "N")
  q <- q + ggtitle("F1 Score")

  ggsave(file.path(outdir, "ngram_stats_solo_F1.pdf"), width = 9.03, heigh = 5.73, units = "in")
  messagef("F1...done.")
  q <- ngram_stats[[2]] %>% ggplot(aes(x = factor(target_n), y = prec, fill = factor(threshold)))
  q <- q + geom_boxplot(fill = def_colour) + geom_point(alpha = .1) + geom_violin(alpha = .1, fill = def_colour)
  #q <- q + facet_wrap(~threshold, ncol = 1)
  q <- q + theme_bw() + labs(x = "N")
  q <- q + ggtitle("Precision")
  ggsave(file.path(outdir, "ngram_stats_solo_prec.pdf"), width = 9.03, heigh = 5.73, units = "in")
  messagef("Precision...done.")
  q <- ngram_stats[[2]] %>% ggplot(aes(x = factor(target_n), y = rec, fill = factor(threshold)))
  q <- q + geom_boxplot(fill = def_colour) + geom_point(alpha = .1) + geom_violin(alpha = .1, fill = def_colour)
  #q <- q + facet_wrap(~threshold, ncol = 1)
  q <- q + theme_bw() + labs(x = "N")
  q <- q + ggtitle("Recall")
  ggsave(file.path(outdir, "ngram_stats_solo_rec.pdf"), width = 9.03, heigh = 5.73, units = "in")
  messagef("Recall...done.")
  #toc()

  #tic()
  ngrs <- ngram_stats[[3]] %>% ungroup() %>% rename(n = target_n, thresh = threshold)
  vars <- c("thresh", "n", sort(names(ngrs)[grep("F1_mean", names(ngrs))]))
  print(ngrs %>% select(vars) %>% arrange(n, desc(F1_mean)))
  vars <- c("thresh", "n", sort(names(ngrs)[grep("prec", names(ngrs))]))
  print(ngrs %>% select(vars) %>% arrange(n, desc(prec_mean)))
  vars <- c("thresh", "n", sort(names(ngrs)[grep("^rec", names(ngrs))]))
  print(ngrs %>% select(vars) %>% arrange(n, desc(rec_mean)))
  #pattern similarity
  #toc()

  #tic()
  pat_sim_eval <- cv_pattern_sim(num_folds = 10, size = 10)
  write_stats(pat_sim_eval[[1]], outdir = outdir, name = "pat_sim_eval_raw.csv", file_format = file_format )
  write_stats(pat_sim_eval[[2]], outdir = outdir, name = "pat_sim_eval_sum.csv", file_format = file_format )
  print(pat_sim_eval[[2]])
  #toc()
  #toc()
  #q
}
