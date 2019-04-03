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
theta <- function(x){ ifelse( x > 0, 1, 0)}