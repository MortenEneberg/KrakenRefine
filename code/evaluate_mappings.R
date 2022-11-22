
WD<-snakemake@input[[5]]
setwd(WD)

libpath<-snakemake@input[[1]]
.libPaths(c(libpath, .libPaths()))

suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(RSQLite)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(svglite)))

krakendatapath<-snakemake@input[[2]]

sample = tools::file_path_sans_ext(basename(krakendatapath))

kraken_report<-read.csv2(paste(krakendatapath, sep=""), header = F, comment.char = "#", sep = "\t")  
colnames(kraken_report)<-c("percentage", "cladeReads", "taxonReads", "n_kmers", "n_unique_kmers",
                    "taxRank", "kraken_genus_taxid", "name")

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



#### LOADING FLEXTAXD SQLITE DATABASE ####
dbfile<-snakemake@input[[6]]

sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = dbfile)

## Loading the different database components into dataframes
genomes2taxid <- dbReadTable(db, "genomes")
nodes<-dbReadTable(db, "nodes")
rank<-dbReadTable(db, "rank")
tree<-dbReadTable(db, "tree")

accession2genome<-snakemake@input[[4]]
##Loading accession2genome file
accession2genome<- read.csv2(file = accession2genome, header = F, sep = "\t") %>%
  unite("accession", V2:V4, remove = TRUE, sep=":") %>%
  separate(accession, into = c("accession", "drop"), sep = "\\s", extra = "merge")%>%
  mutate(V1 = as.character(V1),
  genome = basename(V1)) %>%
  separate(genome, into = c("genome", "drop2"), sep = ".fna") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_genomic") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_ASM") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_Viral") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_Genome") %>%
  separate(genome, into = c("genome", "drop3"), sep = "_GRCh38") %>%
  select(accession, genome)
accession2genome["genome"][accession2genome["genome"] == "UniVec_core"] <- "UniVec_core.fna.gz"

accession2genome2taxid<-accession2genome %>% 
  left_join(genomes2taxid, by = "genome") %>%
  left_join(nodes, by = "id") 


#### Loading in Sam files ####
sam_folder<-snakemake@params[[1]]
sam_list<-paste(sam_folder, list.files(sam_folder, pattern = ".sam"), sep = "/")


sam_files<-combine_sam_files(sam_list) %>%
  left_join(accession2genome2taxid, by = "accession") 
  
accession2length="data/accession2length.tsv"

contig_length<-read.csv2(accession2length, header = F, sep="\t")
colnames(contig_length)<- c("accession", "contig_length")

genome2krakentaxid<-sam_files %>%
  filter(flag == 0 | flag == 16) %>% #take only the primary mapped read
  distinct(genome, kraken_genus_taxid)

genomes<-unique(genome2krakentaxid$genome)
contig_length_subset <- contig_length %>%
  left_join(accession2genome2taxid, by = "accession") %>%
  filter(genome %in% genomes) %>%
  left_join(genome2krakentaxid, by = "genome") %>%
  rename(endpos = contig_length) %>%
  group_by(genome, kraken_genus_taxid) %>%
  mutate(end_pos = cumsum(endpos+1),
         start_pos = end_pos - endpos)

windowsize = 8000

sam_files_grouped <- sam_files %>%
  distinct(read, .keep_all = T) %>%
  ungroup() %>%
  left_join(
    contig_length_subset %>% select(kraken_genus_taxid, genome, accession, start_pos, end_pos), by = c("genome", "accession", "kraken_genus_taxid")
    ) %>%
  distinct(accession, read, kraken_genus_taxid, .keep_all = T) %>%
  mutate(start_genome = start_position + start_pos) %>%
  mutate(window = start_genome %/% windowsize + 1) %>%
  group_by(genome, window) %>%
  summarise(reads = n(),
            kraken_genus_taxid = kraken_genus_taxid,
            genome = genome) %>%
  distinct(genome, window, .keep_all = T)



genome2n_window<-contig_length_subset %>% 
  ungroup() %>% 
  group_by(genome) %>% 
  summarise(genome = genome,
            size = max(end_pos),
            kraken_genus_taxid = kraken_genus_taxid) %>% 
  distinct(genome, size, .keep_all = T) %>% 
  mutate(n_windows = size %/% windowsize + 1) %>% 
  expand(genome, window = 1:n_windows) %>%
  mutate(reads = 0) %>%
  left_join(
    contig_length_subset %>% select(genome, kraken_genus_taxid) %>% distinct(genome, .keep_all = T), by = c("genome")
  )

genome_windows<-sam_files_grouped %>% 
  bind_rows(
    genome2n_window
  ) %>%
  arrange(genome, window) %>%
  distinct(genome, window, .keep_all = T) %>%
  group_by(genome, kraken_genus_taxid) 

genome_windows_sum<-sam_files_grouped %>% 
  bind_rows(
    genome2n_window
  ) %>%
  arrange(genome, window) %>%
  distinct(genome, window, .keep_all = T) %>%
  group_by(genome, kraken_genus_taxid) %>%
  summarise(genome = unique(genome),
            kraken_genus_taxid = unique(kraken_genus_taxid),
            reads_tot = sum(reads),
            Indicator = if (max(window)<=65) {
              if (max(reads)/sum(reads)<0.08) {"TP"} else {"FP"}
            } else #Genomes above 65*windowsize
              if (sum(reads)<10) {if (max(reads) > 1) {"FP"} else {"TP"}
            } else if (sum(reads)<20) {if (max(reads) > 2) {"FP"} else {"TP"}
            } else if (sum(reads)<40) {if (max(reads) > 3) {"FP"} else {"TP"}
            } else if (sum(reads)<80) {if (max(reads) > 4) {"FP"} else {"TP"}
            } else if (sum(reads)<120) {if (max(reads) > 5) {"FP"} else {"TP"}
            } else if (sum(reads)<160) {if (max(reads) > 6) {"FP"} else {"TP"}
            } else if (sum(reads)<200) {if (max(reads) > 7) {"FP"} else {"TP"}
            } else if (sum(reads)<240) {if (max(reads) > 8) {"FP"} else {"TP"}
            } else if (sum(reads)<280) {if (max(reads) > 9) {"FP"} else {"TP"}
            } else if (sum(reads)<320) {if (max(reads) > 10) {"FP"} else {"TP"}
            } else if (max(reads)/sum(reads)<0.08) {"TP"} else {"FP"}) %>%
  group_by(kraken_genus_taxid) %>%
  filter(reads_tot == max(reads_tot))

kraken_report_indicator <- kraken_report %>%
  mutate(name = as.character(name)) %>%
  mutate_if(is.character, trimws) %>%
  filter(taxRank == "G") %>%
  filter(name != "Homo") %>%
  filter(cladeReads>3) %>%
  mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)) %>%
  left_join(genome_windows_sum, by = "kraken_genus_taxid") %>%
  filter(!is.na(Indicator)) %>%
  mutate(tick_color = case_when((Indicator == "TP" ~"darkgreen"),
                                (Indicator == "FP"~"red"),
                                (Indicator == "Small genome"~"yellow3")),
         tick_name = glue("<i style='color:{tick_color}'>{name}</i>")) 


kraken_report_indicator$tick_name <- factor(kraken_report_indicator$tick_name, levels=(kraken_report_indicator$tick_name)[order(kraken_report_indicator$cladeReads)])
krakenrefineplot<-ggplot(kraken_report_indicator, aes("",tick_name, fill = cladeReads)) +
  geom_tile() +
  theme_bw() +
  geom_text(aes(label = cladeReads), size = 3) +
  theme(axis.text.y = element_markdown(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colours = brewer.pal(8, "YlOrRd"),
                       values=rescale(c(0, 5, 250, 1000, 2500, 5000, 10000, 20000))) + 
  labs(x="", y = "Genus", fill = "Reads", title = "Grouped by KrakenRefine") +
  scale_color_manual(name = "KrakenRefine")

krakenrefineplot

kraken_report_plot <- kraken_report %>%
  mutate(name = as.character(name)) %>%
  mutate_if(is.character, trimws) %>%
  filter(taxRank == "G") %>%
  filter(name != "Homo") %>%
  filter(cladeReads>3) %>%
  mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)) %>%
  mutate(tick_color = "black",
         tick_name = glue("<i style='color:{tick_color}'>{name}</i>")) 

kraken_report_plot$tick_name <- factor(kraken_report_plot$tick_name, levels=(kraken_report_plot$tick_name)[order(kraken_report_plot$cladeReads)])


krakennotrefineplot<-ggplot(kraken_report_plot, aes("",tick_name, fill = cladeReads)) +
  geom_tile() +
  theme_bw() +
  geom_text(aes(label = cladeReads), size = 3) +
  theme(axis.text.y = element_markdown(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colours = brewer.pal(8, "YlOrRd"),
                       values=rescale(c(0, 5, 250, 1000, 2500, 5000, 10000, 20000))) + 
  labs(x="", y = "Genus", fill = "Reads", title = "Original classification",
       tick_color = "Legend") +
  scale_color_manual(name = "KrakenRefine")


image<-plot_grid(krakennotrefineplot, krakenrefineplot, labels = "AUTO")

ggsave(file = paste0("data/KrakenRefine/", sample, "/KrakenRefine_", sample, ".svg"), plot = image, width = 20, height = 20)

plot_df<-genome_windows %>%
  mutate(genome = as.character(genome)) %>%
  filter(genome %in% (genome_windows_sum %>% pull(genome))) %>%
  left_join(nodes %>% rename(kraken_genus_taxid = id) %>% mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)), by = "kraken_genus_taxid")%>%
  ungroup() %>%
  group_by(name) %>%
  do(
    plots = ggplot(data = .) + aes(window, reads) +
  geom_point()+ ggtitle(.$name) +
  theme_bw() + 
  labs(x = "Genome window", y = "Read count") 
  )

#pdf(paste("data/KrakenRefine_", sample, ".pdf", sep = ""), onefile = T)
#for(i in seq(length(plot_df$plots))) {
#  do.call("grid.arrange", plot_df$plots[[i]])
#}
#dev.off()