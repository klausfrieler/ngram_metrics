suppressPackageStartupMessages(library(tidyverse))
source("./utils.R")

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
len_pat <- Vectorize (
  function(x){
    length(strsplit(x, ",")[[1]])
  }, SIMPLIFY = T, USE.NAMES = F)

read_evaluations <- function(fname_db ="output_test/eval_v0a_db.csv", fname_wjd = "output_test/eval_v0a_wjd.csv"){
  ev_db <- read.csv(fname_db, sep = ";", stringsAsFactors = F)
  names(ev_db) <- c("value", "len", "min_sim", "max_diff", "n_db", "melid")
  #ev_db <-
  #  ev_db %>%
  #  distinct(value, min_sim, max_diff, len, n_db, melid, .keep_all = T)

  ev_wjd <- read.csv(fname_wjd, sep = ";", stringsAsFactors = F)
  names(ev_wjd) <- c("value", "len", "min_sim", "max_diff", "n_wjd", "melid")
  #browser()
  #ev_wjd <-
  #  ev_wjd %>%
  #  distinct(value, min_sim, max_diff, len, melid, n_wjd, .keep_all = T)

  exact_db <-
    ev_db %>%
    filter(min_sim == 1 & max_diff == 0) %>%
    select(-min_sim, -max_diff) %>%
    #group_by(value, len) %>%
    #summarize(n_db = sum(n_db)) %>%
    #ungroup %>%
    as_tibble()

  exact_wjd <-
    ev_wjd %>%
    filter(min_sim == 1 & max_diff == 0) %>%
    select(-min_sim, -max_diff) %>%
    #group_by(value, len) %>%
    #summarize(n_db = sum(n_wjd)) %>%
    ungroup %>%
    as_tibble()
  ev_wjd <-
    ev_wjd %>%
    filter(min_sim != 1) %>%
    group_by(value, min_sim, max_diff, len) %>%
    summarize(n_wjd = sum(n_wjd, na.rm = T)) %>%
    ungroup() %>%
    as_tibble()
  assign("ev_wjd", ev_wjd, globalenv())

  ev_db <-
    ev_db %>%
    filter(min_sim != 1) %>%
    #distinct(value, min_sim, max_diff, len, melid, .keep_all = T) %>%
    group_by(value, min_sim, max_diff, len) %>%
    summarize(n_db = sum(n_db, na.rm = T)) %>%
    ungroup() %>%
    as_tibble()
  assign("ev_db", ev_wjd, globalenv())

  ev_comb <-
    ev_wjd %>%
    left_join(ev_db, by = c("value", "min_sim", "max_diff", "len")) %>%
    mutate(n_db = replace(n_db, is.na(n_db), 0)) %>%
    mutate(n_diff = n_db - n_wjd,
           rel_diff = (n_db - n_wjd)/n_wjd,
           abs_rel_diff = abs(rel_diff))

  exact_comb <-
    exact_wjd %>%
    left_join(exact_db, by = c("value", "len", "melid")) %>%
    mutate(n_db = replace(n_db, is.na(n_db), 0)) %>%
    distinct(value, len, melid, .keep_all = T) %>%
    mutate(n_diff = n_db - n_wjd,
           rel_diff = (n_db - n_wjd)/n_wjd,
           abs_rel_diff = abs(rel_diff))

  list(exact = exact_comb, collateral = ev_comb)
}
safe_prec <- Vectorize(function(TP, FP){
  if(TP == 0 && FP == 0){
    return(0)
  }
  TP/(TP + FP)
})
get_retrieval_scores <- function(data){
  scores <- data %>%
    group_by(value, melid, len) %>%
    mutate(TP = n_db,
           FP = n_diff * theta(n_diff),
           FN  = abs(n_diff) * theta(-n_diff)) %>%
    summarise( TP = sum(TP), FP = sum(FP), FN = sum(FN)) %>%
    ungroup() %>%
    group_by(value, len) %>%
    summarise( TP = sum(TP), FP = sum(FP), FN = sum(FN)) %>%
    mutate(prec = safe_prec(TP, FP), rec = TP/(TP + FN), F1 = 2 * TP/(2 * TP + FP + FN)) %>%
    ungroup()
  scores
}

query_evaluation <- function(test_set_dir, outdir = test_set_dir, file_format = "en"){
  fname_db <- file.path(test_set_dir, "retrieval_eval_db.csv")
  fname_wjd <- file.path(test_set_dir, "retrieval_eval_wjd.csv")
  ev <- read_evaluations(fname_db = fname_db, fname_wjd = fname_wjd)
  #browser()
  assign("ev", ev, globalenv())
  retrieval_eval <-
    get_retrieval_scores(ev$exact) %>%
    summarise_at(vars(prec:F1), list(~mean, ~median, ~sd))
  write_stats(retrieval_eval, outdir = outdir, name = "retrieval_eval.csv", file_format = file_format )
  collateral_eval <-
    ev$collateral %>%
    group_by(len) %>%
    summarise_at(vars(rel_diff:abs_rel_diff), list(~mean, ~median, ~sd, ~min, ~max))
  write_stats(collateral_eval, outdir = outdir, name = "collateral_eval.csv", file_format = file_format )
  q <- ev$collateral %>% ggplot(aes(x = factor(len), y = rel_diff))
  q <- q + geom_boxplot(alpha = .2, fill = "black") + geom_point(alpha = .1)
  q <- q + labs(x = "Query Length", y = "Mean Rel. Result Set Size")
  q <- q + theme_bw()
  ggsave(file.path(outdir, "rel_diff_box.pdf"), width = 9.03, heigh = 5.73, units = "in")
  q
}