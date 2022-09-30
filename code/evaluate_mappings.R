libpath<-snakemake@input[[2]]
.libPaths(c(libpath, .libPaths()))

suppressMessages(suppressWarnings(library(dplyr,warn.conflicts = F, quietly = T)))

krakendatapath<-snakemake@input[[1]]

read_lengths<-read_tsv(file = input[[4]], 
                        col_names = c("read", "length")) %>%
separate(read, into = c("read", NA), sep=" ")

read_sam <- function(sample, read_lengths) {
  file<-sample, ".sam", sep = "")
  sam_read<-read_tsv(file = file, 
           col_names = c("barcode", "flag", "read", "start_position", "mapq", "cigar", "drop1", "drop2", "drop3", "drop4", "drop5", "align_metrics")) %>%
    select(barcode, flag, start_position, mapq, cigar, align_metrics) %>%
    filter(flag != 4) %>%
    separate(align_metrics, into = c(NA, NA, "align_score", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), sep = ":") %>%
    separate(align_score, into = c("align_score", NA)) %>%
    mutate(sample = sample)%>%
    mutate(end_position = start_position+24) # differentiate for adapters and barcodes!
  
  
  return(sam_read)
}


sam_read<-sam_read %>%
    left_join(read_lengths, by = "read")
  