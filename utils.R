write_stats <- function(stats_df, outdir = ".", name, file_format = "en", RDS_copy = T){
  sep <- ","
  dec <- "."
  if(file_format == "de"){
    sep <- ";"
    dec <- ","
  }
  if(RDS_copy){
    rds_name <- file.path(outdir, gsub("\\.csv", ".RDS", name))
    saveRDS(stats_df, file = rds_name)
  }
  write.table(stats_df,
              file.path(outdir, name),
              sep = sep,
              dec = dec,
              row.names = F,
              col.names = T,
              quote = F)

}
theta <- function(x){ ifelse( x > 0, 1, 0)}

multi_ggsave <- function(outdir, fname, formats = c("png", "pdf"), width  = 9.03, height = 5.73){
  map(formats, function(x){
    ggsave(file.path(outdir, sprintf("%s.%s", fname, x)), width = width, heigh = height, units = "in")
  }) %>% invisible()
}