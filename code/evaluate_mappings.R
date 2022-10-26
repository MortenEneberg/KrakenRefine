setwd("/user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter/")

library(tidyr)
library(RSQLite)


libpath<-"/user_data/men/sepseq/R_lib/x86_64-pc-linux-gnu-library/3.6"
#libpath<-snakemake@input[[2]]
.libPaths(c(libpath, .libPaths()))

suppressMessages(suppressWarnings(library(dplyr,warn.conflicts = F, quietly = T)))

#krakendatapath<-snakemake@input[[3]]
krakendatapath<-"data/kraken2_classification/reads_90/reads_90.report"


kraken_report<-read.csv2(paste(krakendatapath, sep=""), header = F, comment.char = "#", sep = "\t")  
colnames(report)<-c("percentage", "cladeReads", "taxonReads", "n_kmers", "n_unique_kmers",
                    "taxRank", "taxID", "name")

#### DEFINING FUNCTIONS ####
read_sam <- function(file) {
  sam_read<-read.csv2(file = file, header = F, sep = "\t")%>%
    select(V1, V2, V3, V4, V5, V6)
  colnames(sam_read)<-c("read", "flag", "accession", "start_position", "mapq", "cigar")
  sam_read<-sam_read%>%  
  filter(flag != 4) %>%
  mutate(kraken_genus_taxid = tools::file_path_sans_ext(basename(file)))
  
  return(sam_read)
}

combine_sam_files<-function(file_list){
  combined_sam<-data.frame()
  for (i in file_list) {
    sam_batch<-read_sam(i)
    combined_sam<-combined_sam %>%
      bind_rows(
        sam_batch
      )
  }
  return(combined_sam)
}


#sam_folder = input[[5]]
sam_folder<-"data/mapped_reads/reads_90/"
sam_list<-paste(sam_folder, list.files(sam_folder, pattern = ".sam"), sep = "")

sam_files<-combine_sam_files(sam_list)
  
#### LOADING FLEXTAXD SQLITE DATABASE ####
#dbfile<-snakemake@input[[6]]
dbfile<-"/user_data/men/sepseq/databases/2022_10_18_kraken2_EUPATH_database/databases/NCBI_GTDB_merge.db"

sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = dbfile)

## Loading the different database components into dataframes
genomes2taxid <- dbReadTable(db, "genomes")
nodes<-dbReadTable(db, "nodes")
rank<-dbReadTable(db, "rank")
tree<-dbReadTable(db, "tree")
