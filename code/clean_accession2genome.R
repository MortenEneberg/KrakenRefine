
##Snakemake inputs
WD<-snakemake@input[[1]] 
setwd(WD)
libpath<-snakemake@input[[3]]
accession2genome<-snakemake@input[[2]]

.libPaths(libpath)

suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readr)))

##Loading accession2genome file
accession2genome<- read_tsv(accession2genome, col_names = F) %>%
  unite("accession", X2:X4, remove = TRUE, sep=":") %>%
  separate(accession, into = c("accession", "drop"), sep = "\\s", extra = "merge")%>%
  mutate(V1 = as.character(X1),
         genome = basename(X1)) %>%
  separate(genome, into = c("genome", "drop2"), sep = ".fna") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_genomic") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_ASM") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_Viral") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_Genome") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_GRCh38") %>%
  select(accession, genome)
accession2genome["genome"][accession2genome["genome"] == "UniVec_core"] <- "UniVec_core.fna.gz"

write_tsv(accession2genome, paste0(WD, "data/accession2genome_cleaned.tsv"), col_names = T)
