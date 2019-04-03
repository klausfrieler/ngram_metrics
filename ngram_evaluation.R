suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
suppressPackageStartupMessages(suppressWarnings(library(tictoc)))
source("similarity.R")

printf <- function(...) print(sprintf(...))
messagef <- function(...) message(sprintf(...))

read_mcsv2 <- function(fname, vars = c("onset", "pitch")){
  #messagef("Reading %s ", fname)
  f <- read.csv(fname, header = T, sep=";")
  as_tibble(f[, vars])
}
read_transcription <- function(fname){
  #messagef("Reading %s ", fname)
  f <- read.csv(fname, header = T, sep=",")
  as_tibble(f[, c("onset", "pitch")])
}

setup_workspace <- function(test_set_dir = "data/note_basis/test_set_v2"){
  #read MCSV2
  mcsv2 <- lapply(list.files("data/mcsv2", full.names = TRUE), read_mcsv2)
  names(mcsv2) <- list.files("data/mcsv2", full.names = FALSE)
  mcsv2[["LouisArmstrong_CornetChopSuey_FINAL.csv"]] <- NULL
  names(mcsv2) <- gsub("=", "-", names(mcsv2))
  assign("mcsv2", mcsv2, envir = globalenv())
  messagef("Read %d WJD solos", length(mcsv2))
  #read transcriptions
  note_basis <- lapply(list.files(test_set_dir, full.names = TRUE), read_transcription)
  names(note_basis) <- gsub("Solo_note", "FINAL",
                            list.files(test_set_dir, full.names = FALSE))
  assign("note_basis", note_basis, envir = globalenv())
  messagef("Read %d extracted melodies", length(note_basis))

  wjd <- read.csv("data/WJD_metadata.csv", header = T, sep = ";", stringsAsFactors = F) %>% as_tibble()
  wjd$id2 <- gsub("\\.sv", "\\.csv", wjd$id)
  wjd <- wjd[, c("avgtempo", "instrument", "performer", "recordingyear", "rhythmfeel", "style", "tempoclass", "tonality_type", "id2")]
  wjd$instrument_red <- fct_collapse(wjd$instrument, sax = c("ss", "as", "ts", "ts-c", "bs", "cl", "bcl"), brass = c("tp", "tb", "cor"), g = "g", vib = "vib", p = "p")
  assign("wjd", wjd, envir = globalenv())
  assign("common", intersect(names(mcsv2), names(note_basis)), globalenv())
}

test_get_ngram <- function(x, n){
  ret <- as.character(x)
  for(i in 1:(n-1)){
    ret <- cbind(ret, lead(x, i))
  }
  ret
}
test_get_ngram2 <- function(x, n){
  x
  map_dfc(0:n, function(i) lead(x, i))
}

test_row_access <-function(ret){
  nrow <- nrow(ret)
  x <- sample(1:nrow, 1)
  if(!any(is.na(ret[x, ])))paste(ret[x,], collapse=",")}

get_ngram <- function(x, n, id = "NA"){
  ret <- as.character(x)
  if(n <= 0){
    stop("n must be a positive integer")
  }
  if(n == 1){
    ret <- tibble(id = id, pos = 1:length(ret), n = n, value = ret)
    return(ret)
  }
  #tic()
  for(i in 1:(n-1)){
      ret <- cbind(ret, lead(x, i))
  }
  #toc()
  #tic()
  #ret <- map(0:(n-1), function(i) lead(x, i)) %>% bind_cols()
  #toc()
  ret <- map_chr(1:(nrow(ret) - n + 1), function(x) if(!any(is.na(ret[x, ])))paste(ret[x,], collapse=","))
  ret <- tibble(id = id, pos = 1:length(ret), n = n, value = ret)
  return(ret)
  #ret <- as.tibble(ret)
  #names(ret) <- sprintf("lag%d", 0:(n-1))
  #ret <- map_chr(1:(nrow(ret)-n), function(x) if(!any(is.na(ret[x, ])))paste(ret[x,], collapse=","))
  #ret
}
get_all_ngrams <- function(x, max_n, id = "NA"){
  map_df(1:max_n, function(n) get_ngram(x, n, id))
}
get_all_ngrams2 <- function(x, max_n, id = "NA"){
  ret <- NULL
  base <- NULL
  for(i in 0:(max_n-1)){
    base <- cbind(base, lead(x, i))
  }
  l <- nrow(base)
  map_dfr(1:max_n, function(n){
        pos <- 1:(l - n + 1)
        tmp <- apply(as.matrix(base[pos, 1:n]), 1, function(row) paste(row, collapse = ","))
        #printf("n = %d, max pos: %d, len(tmp) = %d, class(tmp) = %s", n, l - n + 1, length(tmp), class(tmp))
        tibble(id = id, pos = pos, n = n, value = tmp)
  })
}
prepare_ngram_analysis <- function(note_list, max_n = 10, subset = NULL, global_name = NULL, mean_onset = F){
  if(!is.null(subset)){
    note_list <- note_list[subset]
  }
  n <- length(note_list)
  tic()
  tic()
  tmp <- map_df(1:n, function(i) get_all_ngrams(note_list[[i]]$pitch, max_n, id = names(note_list)[i]))
  toc()
  tic()
  if(mean_onset){
    tmp$onset <- map_dbl(1:nrow(tmp), function(i){
      row <- tmp[i,]
      mean(note_list[[row$id]]$onset[seq(row$pos, row$pos + row$n - 1)])
    } )
  } else {
    tmp$onset <- map_dbl(1:nrow(tmp), function(i) note_list[[tmp[i,]$id]]$onset[tmp[i,]$pos])

  }
  toc()
  tmp$solo <- tmp$id
  tmp$id  <- sprintf("%d-%d-%d", as.integer(factor(tmp$solo)), tmp$n, tmp$pos)
  tic()
  if(!is.null(global_name) && nchar(global_name)> 0){
    tmp$source <- global_name
    assign(x = global_name, value = tmp, globalenv())
  }
  toc()
  toc()
  tmp
}
prepare_ngrams <- function(max_n = 10){
  common <- intersect(names(mcsv2), names(note_basis))
  load("data/wjd_ngrams_all.rda")
  wjd_ngrams <- wjd_ngrams_all %>% filter(solo %in% common, n <= max_n)
  if(nrow(wjd_ngrams) > 0 ){
    assign("wjd_ngrams", wjd_ngrams, globalenv())
    messagef("Loaded %d WJD ngrams", nrow(wjd_ngrams))
  }
  else {
    prepare_ngram_analysis(mcsv2, max_n, subset = common, global_name = "wjd_ngrams", mean_onset = T)
  }
  prepare_ngram_analysis(note_basis, max_n, subset = common, global_name = "db_ngrams", mean_onset = T)
  messagef("Created %d WJD ngrams and %d ME ngrams", nrow(wjd_ngrams), nrow(db_ngrams))
}

find_in_range <- function(x, y, threshold = .05, as_value = F){
  if(length(y) > 1){
    stop()
  }
  pos <- which(abs(x - y) <= threshold)
  if(length(pos) == 0){
    pos <- NA
  }
  if(as_value){
    ret <- tibble(pos = pos, value = x[pos])
    return(x[pos])
  }
  pos
}

find_closest_elements <- function(source, target, threshold = .05){
  r <- 1:length(target$onset)
  #r <- 1:10
  #printf("Calling find_closest_elements with len source %d target %d", nrow(source), nrow(target))
  ret <- map_dfr(r, function(x) {
    #printf("testing %d", r)
    pos <- find_in_range(source$onset, target$onset[x], threshold = threshold, as_value = F);
    if(length(pos) == 0 || all(is.na(pos))){
      #return(NULL)
      pos <- NA
      tmp1 <- source[1, c("solo", "id", "pos", "onset", "value", "n")]
      tmp1$pos <- NA
      tmp1$pitch <- NA
      tmp1$onset <- NA
      tmp1$value <- NA
    }
    else {
      tmp1 <- source[pos, c("solo", "id", "pos", "onset", "value", "n")]
    }
    tmp2 <- target[x,  c("solo", "id", "pos", "onset", "value", "n")]
    names(tmp1) <- sprintf("source_%s", names(tmp1))
    names(tmp2) <- sprintf("target_%s", names(tmp2))


    cbind(tmp1, tmp2)
  })
  ret$d_onset <- round((ret$source_onset - ret$target_onset)*1000)
  #ret$d_n <- ret$source_n - ret$target_n
  ret$equal <- ret$source_value == ret$target_value
  ret$equal[is.na(ret$equal)] <- FALSE
  ret
}

get_ngram_analysis <- function(recalc = F,
                               max_n = 10,
                               thresholds = seq(.03, .07, .02),
                               subset = NULL,
                               add_inverse = F){
  common <- intersect(names(mcsv2), names(note_basis))
  if(recalc){
    prepare_ngrams(max_n = max_n)
  }
  if(!is.null(subset)){
    common <- intersect(common, subset)
  }
  library(furrr)
  plan(multiprocess)

  cmp <- list()
  for(n in 1:max_n){
    db_tmp <- db_ngrams[db_ngrams$n == n, ]
    wjd_tmp <- wjd_ngrams[wjd_ngrams$n == n,]
    if(nrow(db_tmp) == 0){
      printf("Data does not contain ngrams of length %d", n)
      break
    }
    for (t in thresholds){
      tic()
      tmp <- future_map_dfr(1:length(common), function(x){
        find_closest_elements(source = db_tmp[db_tmp$solo == common[x],],
                              target = wjd_tmp[wjd_tmp$solo == common[x], ],
                              threshold = t)}, .progress = T)
      toc()
      tmp$threshold <- t
      tmp$direction <- "target"
      cmp[[sprintf("target_%s_%d", as.character(t), n)]] <- tmp
      if(add_inverse || n == 1){
        tmp <- map_df(1:length(common), function(x){
          find_closest_elements(source = wjd_tmp[wjd_tmp$solo == common[x], ],
                                target = db_tmp[db_tmp$solo == common[x], ],
                                threshold = t)})
        tmp$threshold <- t
        tmp$direction <- "source"
        cmp[[sprintf("source_%s_%d", as.character(t), n)]] <- tmp

      }
    }
  }
  bind_rows(cmp)
}
get_dtw_distance <- function(solo, N = 1){
  library(dtw)

  solo <- enquo(solo)
  query <- db_ngrams %>% filter(n == N, solo == !!solo) %>% pull(value) %>% as.integer()
  target <- wjd_ngrams %>% filter(n == N, solo == !!solo) %>% pull(value) %>% as.integer()
  alignment <- dtw(query, target, keep.internals = F, distance.only = T)
  alignment$normalizedDistance
}
get_ngram_stats <- function(ngram_match){
  ngram_raw <- ngram_match %>%
    group_by(threshold, direction, target_n, target_solo, target_pos) %>%
    summarise(accuracy1 = sum(equal),
              accuracy2 = mean(equal),
              missing = sum(is.na(source_pos)),
              num_candidates = sum(!is.na(source_pos)))
  ngram_sum <- ngram_raw %>%
    group_by(threshold, direction, target_solo) %>%
    summarise(accuracy1 = mean(accuracy1),
              accuracy2 = mean(accuracy2),
              missing = mean(missing),
              num_candidates = mean(num_candidates),
              multiples = mean(num_candidates > 1))
  ngram_total <- ngram_sum %>%
    group_by(threshold, direction) %>%
    select(-target_solo) %>%
    summarise_all(list(mean = mean, sd = sd))
  list(ngram_raw, ngram_sum, ngram_total)
}

get_ngram_stats_retrieval <- function(ngram_match){
  ngram_raw <- ngram_match %>%
    group_by(threshold, direction, target_solo, target_n, target_pos) %>%
    summarise(TP = sum(equal),
              FP = sum(!equal & !is.na(source_pos)),
              FN = sum(is.na(source_pos)))
  assign("ngram_raw", ngram_raw, globalenv())
  extra_FP <-
    ngram_raw %>%
    filter(direction == "source", target_n == 1) %>%
    group_by(threshold, target_solo) %>%
    summarise(FP = sum(FN))
  assign("extra_FP", extra_FP, globalenv())
  ngram_raw <- ngram_raw %>% filter(direction == "target")
  ngram_sum <- ngram_raw %>%
    group_by(threshold, direction, target_n, target_solo) %>%
    summarise(TP = sum(TP == 1),
              FP = sum(FP > 0),
              FN = sum(FN))
  assign("ngram_sum", ngram_sum, globalenv())
  #ngram_sum[ngram_sum$]
  ngram_sum[ngram_sum$target_n == 1,]$FP <- ngram_sum[ngram_sum$target_n == 1,]$FP + extra_FP[["FP"]]
  ngram_sum <-  ngram_sum %>%
    group_by(threshold, direction, target_n, target_solo) %>%
    summarise(
      prec = TP/(TP + FP),
      rec = TP/(TP + FN),
      F1 = 2*TP/(2*TP + FP + FN))
  assign("ngram_sum2", ngram_sum, globalenv())
  ngram_total <- ngram_sum %>%
    group_by(threshold, direction, target_n) %>%
    select(-target_solo) %>%
    summarise_all(list(mean = mean, median = median, sd = sd, max = max, min = min))
  list(ngram_raw, ngram_sum, ngram_total)

}
pitch_match_plot <- function(ngram_match, solo = common[1], threshold = .05, direction = "target"){
  threshold <- enquo(threshold)
  direction <- enquo(direction)
  tmp <-
    ngram_match %>%
    filter(target_n == 1, threshold == !!threshold, target_solo == solo, direction == !!direction) %>%
    select(target_pos, target_onset, target_value, source_onset, source_value, equal) %>%
    mutate(target_value = as.integer(target_value), source_value = as.integer(source_value)) %>%
    mutate(matched = !is.na(source_value))
  aux <-  ngram_stats[[1]]  %>%
    filter(target_n == 1, threshold == !!threshold, target_solo == solo, direction == !!direction)
  point_class_levels <- c("TP", "FP", "TP + FP", "FN")
  tmp <-
    tmp %>%
    left_join(aux, by = c("target_pos")) %>%
    mutate(point_class = factor((TP>0) + 2*(FP>0)  + 4*FN, labels = point_class_levels))
  print(table(tmp$point_class))
  #print(tmp %>% select(TP, FP, FN, point_class))
  q <-
    tmp %>%
    ggplot(aes(x = target_onset, y = target_value, color = point_class, shape = point_class)) + geom_point() + theme_bw() + geom_line(aes(group=1)) + ggtitle(solo)

  return(q)
  #multiples <- tmp %>% group_by(target_pos) %>%  summarise(n = n()) %>% filter(n > 1) %>% pull(target_pos)
  #tmp <- tmp %>% filter(target_pos %in% multiples[1:5])
  #assign("tmp", tmp, globalenv())
  q <- tmp %>%  ggplot() + geom_segment(aes(x = target_onset,
                           xend = source_onset,
                           y = target_value,
                           yend = source_value))
  q <- q + geom_point(aes(x = target_onset,
                          y = target_value,
                          colour = matched)
                      )
  q <- q + geom_point(aes(x = source_onset, y = source_value),
                      colour = "black",
                      alpha = .2, size = 3)
  q <- q + theme_bw()
  q <- q + labs(x = "Onset", y= "Pitch")
  q
}
pitch_to_interval <- function(x){
    map_chr(strsplit(x, ","), function(y) y %>% as.integer() %>% diff() %>% paste(collapse=","))
}
pitch_ngrams_to_int <- function(data){
  data %>%
  filter(n > 1) %>%
  mutate(value  = pitch_to_interval(data$value)) %>%
  mutate(n  = n - 1)
}
global_pattern_sim_cmp <- function(data1, data2, N_range = 1:10, sim_measure = "total_variation"){
    sim_measure <- match.arg(sim_measure, c("total_variation", "mi", "jsd", "jaccard"))
    N_range <- setdiff(1:max(data1$n), N_range)
    N_range <- setdiff(1:max(data2$n), N_range)
    ret <- list()
    for(N in N_range){
      printf("Checking %d", N)
      set1 <- data1[data1$n == N,]$value
      set2 <- data2[data2$n == N,]$value
      sim_val <- distribution_similarity(set1, set2, type = sim_measure)
      ret[[as.character(N)]] <- tibble(N = N,
                                       wjd_ngrams = n_distinct(set1),
                                       db_ngrams = n_distinct(set2),
                                       common = length(intersect(set1, set2)),
                                       total = length(union(set1, set2)),
                                       in_wjd_not_db = length(setdiff(set1, set2)),
                                       in_db_not_wjd = length(setdiff(set2, set1)),
                                       jaccard_sim =  jaccard_sim(set1, set2),
                                       sim_val = sim_val)
  }
  bind_rows(ret)
}

pattern_sim_cmp <- function(subset = common, N_range = 1:10, summary = T){
  sim_db <-
    pattern_sim(db_ngrams %>% filter(solo %in% subset),
                "solo",
                val1 = unique(subset),
                N_range = N_range) %>%
    rename(common_db = common)

  sim_wjd <-
    pattern_sim(wjd_ngrams %>% filter(solo %in% subset),
                "solo",
                val1 = unique(subset),
                N_range = N_range) %>%
    rename(common_wjd = common)
  sim_cmp <-
    sim_wjd %>%
    left_join(sim_db, by = c("N", "id1", "id2")) %>%
    mutate(d_sim = common_wjd - common_db,
           d_sim_rel = (common_wjd - common_db)/common_wjd
           ) %>%
    filter(id1 != id2)
  #print(summary(sim_cmp$d_sim_rel))
  #return(sim_cmp)

  ret <-
    sim_cmp %>%
    filter(id1 != id2)

  if(summary){
    q <-
      sim_cmp %>%
      ggplot(aes(x = d_sim, y = ..count.. ))
    q <- q + geom_histogram() + facet_wrap(~N, scale = "free")
    r <- sim_cmp %>%
      #filter(common_wjd > 1e-3) %>%
      ggplot(aes(x = common_wjd, y = d_sim, colour = factor(N)))
    r <- r  + geom_point()
    r <- r + geom_smooth(method = "lm", formula = "y~poly(x,1)")

    #print(r)
    ret <-
      ret %>%
      group_by(N) %>%
      summarise(cor = cor(common_wjd, common_db),
                mean_common_wjd = mean(common_wjd),
                mean_d = mean(d_sim),
                sd_d = sd(d_sim),
                min_d = min(d_sim),
                max_d = max(d_sim))

  }
  ret
}
cv_pattern_sim <- function(num_folds = 10, size = 10){
  ret <- list()
  ret <- map(1:num_folds,  function(i) {
    messagef("Calculating pattern similarity batch #%d", i)
    subset <- sample(common, size =  size)
    tmp <- pattern_sim_cmp(subset, summary = F)
    tmp$batch <- i
    tmp
  })
  ret <- bind_rows(ret) %>% filter(N <= 6)

  ret_sum <-
    ret %>%
    filter(N <= 6) %>%
    group_by(N) %>%
    mutate(z_common_wjd = scale(common_wjd), z_common_db = scale(common_db), d_z = z_common_wjd - z_common_db) %>%
    summarise(cor = cor(common_wjd, common_db),
              mean_common_wjd = mean(common_wjd),
              mean_d = mean(d_sim),
              rel_mean_d = mean_d /mean_common_wjd,
              sd_d = sd(d_sim),
              min_d = min(d_sim),
              max_d = max(d_sim),
              mean_abs_dz = mean(abs(d_z)),
              sd_abs_dz = sd(abs(d_z))) %>%
    ungroup()
  print(mean(ret_sum$sd_abs_dz))
  ret_pooled <-
    ret_sum  %>%
    summarise(mdz = mean(mean_abs_dz), sdz = mean(sd_abs_dz))
  ret_sum <- bind_rows(ret_sum, tibble(N = -1,
                                   cor = NA,
                                   mean_common_wjd = NA,
                                   mean_d = NA,
                                   rel_mean_d = NA,
                                   sd_d = NA,
                                   min_d = NA,
                                   max_d = NA,
                                   mean_abs_dz = ret_pooled$mdz,
                                   sd_abs_dz = ret_pooled$sdz)) %>%
    select(N, cor, mean_abs_dz, mean_d, rel_mean_d, mean_common_wjd, sd_d:max_d)
  list(raw = ret, summary = ret_sum)
}