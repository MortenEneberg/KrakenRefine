setwd("/user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter/")

libpath<-"/user_data/men/sepseq/R_lib/x86_64-pc-linux-gnu-library/3.6"
#libpath<-snakemake@input[[2]]
.libPaths(c(libpath, .libPaths()))

library(tidyr)
library(dplyr)
library(RSQLite)

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


read_bed <- function(file) {
  bed_read<-read.csv2(file = file, header = F, sep = "\t")
  colnames(bed_read)<-c("accession", "coverage", "cov_bases", "bases", "cov_perc")
  bed_read<-bed_read%>%
    mutate(kraken_genus_taxid = tools::file_path_sans_ext(basename(file)))
  
  return(sam_read)
}

combine_bed_files<-function(file_list){
  combined_bed<-data.frame()
  for (i in file_list) {
    bed_batch<-read_bed(i)
    combined_bed<-combined_bed %>%
      bind_rows(
        bed_batch
      )
  }
  return(combined_bed)
}

#### LOADING FLEXTAXD SQLITE DATABASE ####
#dbfile<-snakemake@input[[6]]
dbfile<-"/user_data/men/sepseq/databases/2022_10_18_kraken2_EUPATH_database/databases/NCBI_GTDB_merge.db"

sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = dbfile)

genomes2taxid <- dbReadTable(db, "genomes")
nodes<-dbReadTable(db, "nodes")
rank<-dbReadTable(db, "rank")
tree<-dbReadTable(db, "tree")

## Loading the different database components into dataframes
accession2genome<- read.csv2(file = "data/mapped_reads/accession2genome.tsv", header = F, sep = "\t") %>%
  unite("accession", V2:V4, remove = TRUE, sep=":") %>%
  separate(accession, into = c("accession", "drop"), sep = "\\s", extra = "merge")%>%
  mutate(V1 = as.character(V1),
  genome = basename(V1)) %>%
  separate(genome, into = c("genome", "drop2"), sep = ".fna") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_genomic") %>%
  select(accession, genome)
accession2genome["genome"][accession2genome["genome"] == "UniVec_core"] <- "UniVec_core.fna.gz"

accession2genome2taxid<-accession2genome %>% 
  left_join(genomes2taxid, by = "genome") %>%
  left_join(nodes, by = "id") 



#sam_folder = input[[5]]
sam_folder<-"data/mapped_reads/reads_90/"
sam_list<-paste(sam_folder, list.files(sam_folder, pattern = ".sam"), sep = "")

sam_files<-combine_sam_files(sam_list) %>%
  left_join(accession2genome2taxid, by = "accession")
  group_by(kraken_genus_taxid) 
#%>%
 #   summarise(length(unique(accession)))
  
  
bed_folder<-"data/coverage/reads_90/"
bed_list<-paste(sam_folder, list.files(sam_folder, pattern = ".bed"), sep = "")

bed_files<-combine_bed_files(bed_list) %>%
  left_join(accession2genome2taxid, by = "accession")
