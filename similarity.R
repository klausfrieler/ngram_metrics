library(tidyverse)
#source("utf8_test.R", encoding = "UTF-8")
jaccard_sim <- function(x, y){
  #pat1 <- set1$value
  #pat2 <- set2$value
  lu <- length(union(x, y))
  if(lu == 0){
    return(0)
  }
  length(intersect(x, y))/lu
}
weighted_jaccard_sim <- function(x, y, wx, wy){
  tx <- tibble( x = x, wx = wx)
  ty <- tibble( x = y, wy = wy)
  common <- tx %>% inner_join(ty, by = "x") %>% distinct(x, .keep_all = T)
  all <- tx %>% full_join(ty, by = "x") %>% distinct(x, .keep_all = T)
  if(!identical(common$wx, common$wy)){
    stop("Weights do not match for common elements")
  }
  all[is.na(all$wx),]$wx <- all[is.na(all$wx),]$wy

  #print(common)
  #print(all)
  sum(common$wx)/sum(all$wx)
}
compare_sims <- function(data, N, melid1, melid2, weight_var = "freq"){
  d1 <- data[data$N %in% N & data$melid == melid1 & data$independent, ]
  d2 <- data[data$N %in% N & data$melid == melid2 & data$independent, ]
  wsim <- weighted_jaccard_sim(d1$value, d2$value, d1[, weight_var], d2[, weight_var])
  sim <- jaccard_sim(d1$value, d2$value)
  tibble(id1 = melid1, id2 = melid2, N = paste(N), sim = sim, wsim = wsim)
}
distribution_similarity <- function(x, y, type = "total_variation"){
  type <- match.arg(type, c("total_variation", "mi", "jsd", "jaccard", "all"))
  if(length(x) == 0 || length(y) == 0){
    return(0)
  }
  joint <- c(x, y) %>% table()
  tx <- factor(x, levels= names(joint)) %>% table()
  ty <- factor(y, levels = names(joint)) %>% table()
  if(type %in% c("mi", "jsd")){
    #avoid zero counts
    tx[tx == 0] <- 1/sum(tx)
    ty[ty == 0] <- 1/sum(ty)

  }
  #print(joint[joint > 1])
  joint <- joint/sum(joint)
  tx <- tx/sum(tx)
  ty <- ty/sum(ty)

  #print(tx)
  #print(tx[names(joint)[joint > 0]])
  mi <- 1 - mean(joint * log(joint) - tx*log(tx) - ty*log(ty))
  m <- .5*(tx + ty)
  jsd <-  1 - .5 * mean(tx * log(tx / m) + ty * log( ty / m))
  tv <- 1 - .5 * sum(abs(tx - ty))
  jaccard <- jaccard_sim(x, y)
  if(type == "all"){
    return(tibble(jsd = jsd, total_variation = tv, mi = mi, jaccard = jaccard))
  }
  if(type == "mi"){
    return(mi)
  }
  if(type == "jsd"){
    return(jsd)
  }
  if(type == "jaccard"){
    return(jaccard)
  }
  tv

}
mutual_pattern_share <- function(data, melid1, melid2, N_range ){
  set1 <- data[data$solo == melid1,]
  set2 <- data[data$solo == melid2,]
  set1$n <- NULL
  set2$n <- NULL

  print(nrow(set1))
  ret <- list()
  for(n in N_range){
    pat1 <- set1 %>% filter(n == n) %>% pull(value)
    pat2 <- set2 %>% filter(n == n) %>% pull(value)
    printf("n = %d, len pat1: %d pat2: %d", n, length(pat1), length(pat2))
    ret[[n]]<- tibble(N = n,
                      id1 = melid1,
                      id2 = melid2,
                      common = length(intersect(pat1, pat2))/length(union(pat1, pat2)))

  }
  bind_rows(ret)
}
pattern_sim <- function(data,
                        group_var,
                        val1,
                        val2 = NULL,
                        N_range,
                        sim_measure = "total_variation",
                        only_independent = T){
  sim_measure <- match.arg(sim_measure, c("total_variation", "mi", "jsd", "jaccard"))
  add_diagonal <- F
  if(is.null(val2)){
    #full matrix
    val2 <- val1
    add_diagonal <- T
  }
  #cross_indices <- cross2(val1, val2, .filter = function(x, y) ifelse(x == y, TRUE, FALSE))
  cross_indices <- cross2(val1, val2, .filter = function(x, y) ifelse(x <= y, TRUE, FALSE))
  #print(cross_indices)
  #if(length(cross_indices) == 0){
  #  return(1)
  #}
  N_range <- intersect(1:max(data$n), N_range)
  ret <- list()
  for(N in N_range){
    #printf("Checking %d", N)
    tmp <- data[data$n == N,]
    #if(only_independent){
    #   tmp <- tmp[tmp$independent,]
    #}
    sim_func <- function(x){
      set1 <- tmp[tmp[, group_var] == x[[1]][1],]$value
      set2 <- tmp[tmp[, group_var] == x[[2]][1],]$value
      #printf("Length: %d %d", length(set1), length(set2))
      distribution_similarity(set1, set2, type = "total_variation")
    }
    sim_val <- map_dbl(cross_indices, sim_func)
    #printf("Sin vaL. %f", sim_val)
    #stop()
    ret[[as.character(N)]] <- tibble(N = N,
                                     id1 = map_chr(cross_indices, function(x) as.character(x[[1]][1])),
                                     id2 = map_chr(cross_indices, function(x) as.character(x[[2]][1])),
                                     common = sim_val)
    ret[[sprintf("%d_sym",N)]] <- tibble(N = N,
                                         id1 = map_chr(cross_indices, function(x) as.character(x[[2]][1])),
                                         id2 = map_chr(cross_indices, function(x) as.character(x[[1]][1])),
                                         common = sim_val)

    if( add_diagonal ) {
      ret[[sprintf("diag_%d", N)]] <- tibble(N = N, id1 = as.character(val1), id2 = as.character(val2), common = 1)
      #print( ret[[sprintf("diag_%d", N)]] )
      #stop()
    }
  }
  bind_rows(ret) #%>%
  #distinct(id1, id2, N, .keep_all =  T)
  #ret
}
augment_solo_sim_mat <- function(solo_sim_mat){
  solo_sim_mat$id1 <- as.integer(solo_sim_mat$id1)
  solo_sim_mat$id2 <- as.integer(solo_sim_mat$id2)
  ret <-solo_sim_mat %>%
    mutate(melid = id1) %>%
    left_join(year_id_map) %>%
    rename(year1 = recordingyear) %>%
    left_join(performer_id_map) %>%
    rename(performer1 = performer) %>%
    mutate(melid = id2) %>%
    left_join(year_id_map) %>%
    rename(year2 = recordingyear) %>%
    left_join(performer_id_map) %>%
    rename(performer2 = performer)
  ret$melid <- NULL
  ret
}

aggregate_sim_df <- function(data,  as_distance = F, method = "select", params){
  if(as_distance){
    data$common <- 1 - data$common
  }
  if(method == "select"){
    data <- data[data$N == params["N"],]
  }
  if(method == "average"){
    data <- data %>% group_by(id1, id2) %>% summarise(common = mean(common * N ^ {params["poly_exp"]}, na.rm = T)) %>% ungroup()
  }
  if(method == "log_average"){
    data <- data %>% group_by(id1, id2) %>% summarise(common = mean(log(1+common), na.rm = T)) %>% ungroup()
  }
  if(method == "scale"){
    data <-
      data %>%
      filter(id1 != id2) %>%
      group_by(N) %>%
      mutate(common = scale(common)[,1]) %>%
      ungroup() %>%
      group_by(id1, id2) %>%
      summarise(common = mean(common)) %>%
      ungroup()
    data <- bind_rows(data, tibble(id1 = unique(data$id1), id2 = unique(data$id1), common = 0) )
  }
  if(method == "power_law"){
    if(as_distance){
      stop("Power law not combinable with distance")
    }
    else{
      data <-
        data %>%
        group_by(id1, id2) %>%
        summarise(common = log(mean(common, na.rm = T) + .00001)/mean(N)) %>%
        ungroup() %>%
        mutate(common = linear_norm(common))
    }
  }
  data
}
sim_df_to_mat <- function(data, sim_val = "common", as_distance = FALSE, norm = F, FUN = NULL, N = NULL){
  #assumes N, id1, id2, and sim values
  data <-  data %>% arrange(id1, id2)
  ids1 <- unique(data$id1)
  ids2 <- unique(data$id2)
  if(length(intersect(ids1, ids2)) != length(union(ids1, ids2))){
    stop("Bad format")
  }
  n <-length(ids1)
  ret <- matrix(as.matrix(data[, sim_val]), nrow = n, ncol = n)
  colnames(ret) <- ids1
  rownames(ret) <- ids1
  if(!is_symmetric(ret)){
    #print(ret)
    stop("Sim matrix not symmetric")
  }
  if(as_distance){
    ret <- 1 - ret
  }
  attr(ret, "distance_matrix") <- as_distance
  ret
}
sim_matrix <- function(data, as_matrix = FALSE, as_distance = F){
  melids <- unique(data$melid)
  n <- length(melids)
  ret_df <- list()
  if(as_distance){
    ret <- matrix(0, l, l)
  }
  else{
    ret <- matrix(1, l, l)
  }
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      sim <- mutual_pattern_share(data, i, j, 1:max(data$N))
      if(as_distance){
        sim <- 1 - sim
      }
      ret_df[[sprintf("%d_%d", i, j)]] <- sim
      ret_df[[sprintf("%d_%d", j, i)]] <- sim
      ret[i, j] <- sim
      ret[j, i] <- sim
    }
  }
  if(as_matrix){
    attr(ret, "distance_matrix") <- as_distance
    return(ret)
  }
  return(bind_rows(ret_df))
}
filter_sim_mat <- function(sim_mat, def_set, inverse = F){
  if(!inverse){
    return(sim_mat[sim_mat$id1 %in% def_set & sim_mat$id2 %in% def_set,])
  }
  sim_mat[!(sim_mat$id1 %in% def_set & sim_mat$id2 %in% def_set),]
}

compare_performer <- function(sim_mat, data, N, id1 = NULL, id2 = NULL){
  if(is.null(id1)){
    id1 <- sim_mat[1,]$id1
  }
  if(is.null(id2)){
    id2 <- sim_mat[1,]$id2
  }
  printf("Comparing %s <-> %s", id1, id2)
  if(!is.null(N) & all(N > 0)){
    data <- data[data$N %in% N,]
  }
  x <- data[data$performer == id1,]$value
  y <- data[data$performer == id2,]$value
  #print(table(x))
  #print((intersect(x, y)))
  print(distribution_similarity(x, y, type = "all"))
  data[data$value %in% intersect(x, y) & data$performer %in% c(id1, id2),]
}#