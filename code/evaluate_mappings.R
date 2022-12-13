## ----Snakemake inputs-----------------------------------------------------------------------------------------
#mode<-"exp"
mode<-"snakemake"
if (mode == "snakemake"){
WD<-snakemake@input[[5]] 
setwd(WD)
libpath<-snakemake@input[[1]]
krakendatapath<-snakemake@input[[2]]
sam_folder<-snakemake@params[[1]] #Folder with Sam files
accession2length=snakemake@input[[3]]
accession2genome<-snakemake@input[[4]]
dbfile<-snakemake@input[[6]]} else if (mode == "exp") {
  WD<-"/user_data/men/sepseq/KrakenRefine/"
  libpath<-"/user_data/men/sepseq/R_lib/x86_64-pc-linux-gnu-library/4.2"
  dbfile<-"/user_data/men/sepseq/databases/2022_05_25_new_kraken_database/databases/NCBI_GTDB_merge.db"
  
  sample<-"GCF_000960005.1_ASM96000v1_genomic_3"
  krakendatapath<-paste0(WD, "data/kraken2_classification/", sample, "/", sample, ".report")
  sam_folder<-paste0(WD, "data/mapped_reads/", sample)
  accession2genome<-"/user_data/men/sepseq/KrakenRefine/data/accession2genome_cleaned.tsv"
  accession2length<-"/user_data/men/sepseq/KrakenRefine/data/accession2length.tsv"
}



## ----Read libraries-------------------------------------------------------------------------------------------
.libPaths(libpath)

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(RSQLite)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(svglite)))
suppressMessages(suppressWarnings(library(data.table)))


## ----Read Kraken2 reports for a sample------------------------------------------------------------------------
sample = tools::file_path_sans_ext(basename(krakendatapath))
kraken_report<-read_tsv(krakendatapath, col_names = F, comment = "#")  
colnames(kraken_report)<-c("percentage", "cladeReads", "taxonReads", "n_kmers", "n_unique_kmers",
                    "taxRank", "kraken_genus_taxid", "name")
kraken_report<- kraken_report %>%
  mutate(name = as.character(name)) %>%
  mutate_if(is.character, trimws) %>%
  filter(taxRank == "G") %>%
  filter(name != "Homo") %>%
  filter(cladeReads>3)


## ----Defining functions used in this Markdown-----------------------------------------------------------------
read_sam <- function(file) {
  sam_read<-read_tsv(file = file, col_names = F)%>%
    select(X1, X2, X3, X4, X5, X6)
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


## ----Loading flextaxd SQLite database-------------------------------------------------------------------------
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = dbfile)

## Loading the different database components into dataframes
genomes2taxid <- dbReadTable(db, "genomes")
nodes<-dbReadTable(db, "nodes")
rank<-dbReadTable(db, "rank")
tree<-dbReadTable(db, "tree")


## ----Loading accession2genome and accession2length------------------------------------------------------------
##Loading accession2genome file
accession2genome<- read_tsv(accession2genome, col_names = T) 

accession2length<-read_tsv(accession2length, col_names = F)
colnames(accession2length)<- c("accession", "contig_length")


## ----Merging genome2taxid and accession2genome dataframes-----------------------------------------------------
accession2genome2taxid<-accession2genome %>% 
  left_join(genomes2taxid, by = "genome") %>%
  left_join(nodes, by = "id") 


## ----Loading sam files----------------------------------------------------------------------------------------
sam_list<-paste(sam_folder, list.files(sam_folder, pattern = ".sam"), sep = "/") #Make a list of all the Sam files to be loaded

sam_files<-combine_sam_files(sam_list) %>%            #Using the function defined in top to read sam files
  left_join(accession2genome2taxid, by = "accession") #Adding information from the accession2genome2taxid dataframe


## ----Subset to primary mappings and get the genomeIDs and taxIDs of these mappings----------------------------
genome2krakentaxid<-sam_files %>%
  filter(flag == 0 | flag == 16) %>% #take only the primary mapped read
  distinct(genome, kraken_genus_taxid)

genomes<-unique(genome2krakentaxid$genome) #Genomes that this sample has mapped to
accession2length_subset <- accession2length %>% #Get a subset of the accession2length2genome2taxid df by filtering to the genomes that were mapped to
  left_join(accession2genome2taxid, by = "accession") %>%
  filter(genome %in% genomes) %>%
  left_join(genome2krakentaxid, by = "genome") %>%
  rename(endpos = contig_length) %>%
  group_by(genome, kraken_genus_taxid) %>%
  mutate(end_pos = cumsum(endpos+1),
         start_pos = end_pos - endpos) #Adding together all the contigs of a particular genome to create a "complete sequence"


## ----Count reads in each window-------------------------------------------------------------------------------
windowsize = 8000

genome_window_counts <- sam_files %>%
  distinct(read, .keep_all = T) %>%
  ungroup() %>%
  left_join(
    accession2length_subset %>% select(kraken_genus_taxid, genome, accession, start_pos, end_pos), by = c("genome", "accession", "kraken_genus_taxid")
    ) %>%
  distinct(accession, read, kraken_genus_taxid, .keep_all = T) %>%
  mutate(start_genome = start_position + start_pos) %>% #Adding the start position of a read on a contig to the start position of the contig in the "complete sequence" genome
  mutate(window = start_genome %/% windowsize + 1) %>%  #Finding out which genome window the read is placed in
  group_by(genome, kraken_genus_taxid, window) %>%
  summarise(reads = n(),    #Counting the number of reads in each window of each genome
            kraken_genus_taxid = kraken_genus_taxid,
            genome = genome) %>%
  distinct(genome, window, kraken_genus_taxid, .keep_all = T)


## ----Getting rows for all the empty windows where no reads were mapped----------------------------------------
genome_empty_windows<-accession2length_subset %>% 
  ungroup() %>% 
  group_by(genome, kraken_genus_taxid) %>% 
  summarise(genome = genome,
            size = max(end_pos),
            kraken_genus_taxid = kraken_genus_taxid) %>% 
  distinct(genome, kraken_genus_taxid, size, .keep_all = T) %>% 
  mutate(n_windows = size %/% windowsize + 1) %>% 
  expand(window = 1:n_windows) %>%
  mutate(reads = 0) %>%
  left_join(
    accession2length_subset %>% select(genome, kraken_genus_taxid) %>% distinct(genome, .keep_all = T), by = c("genome", "kraken_genus_taxid")
  )


## ----Merging the dataframes of windows with and without reads-------------------------------------------------
genome_windows<-genome_window_counts %>% 
  bind_rows(
    genome_empty_windows
  ) %>%
  arrange(genome, window) %>%
  distinct(genome, window, kraken_genus_taxid, .keep_all = T) %>%
  group_by(genome, kraken_genus_taxid) 


## ----Summarising the counts of reads in each window and applying filters to discriminate TP and FP identifications----
genome_windows_sum<-genome_windows %>%
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
            } else if (sum(reads)<120) {if (max(reads) > 6) {"FP"} else {"TP"}
            } else if (sum(reads)<160) {if (max(reads) > 8) {"FP"} else {"TP"}
            } else if (sum(reads)<200) {if (max(reads) > 10) {"FP"} else {"TP"}
            } else if (sum(reads)<240) {if (max(reads) > 12) {"FP"} else {"TP"}
            } else if (sum(reads)<280) {if (max(reads) > 14) {"FP"} else {"TP"}
            } else if (sum(reads)<320) {if (max(reads) > 15) {"FP"} else {"TP"}
            } else if (max(reads)/sum(reads)<0.08) {"TP"} else {"FP"}) %>%
  group_by(kraken_genus_taxid) %>%
  filter(reads_tot == max(reads_tot)) #Filtering the dataframe to only include 

## ----Filter by read mapping taxonomy--------------------------------------------------------------------------
tax_indicator<-sam_files %>% 
  filter(flag == 0 | flag == 16) %>%
  left_join(kraken_report %>% 
              rename(kraken_name = name) %>%
              select(kraken_genus_taxid, kraken_name) %>%
              mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)), by = c("kraken_genus_taxid")) %>%
  separate(name, into = c("genus_name", "species_name"), sep = "\\s", extra = "merge") %>%
  mutate(genus_name = as.character(genus_name),
         kraken_name = as.character(kraken_name),
         tax_indicator = case_when(kraken_name == genus_name ~ "Same",
                                   kraken_name != genus_name ~ "Different")) %>%
  group_by(kraken_genus_taxid, tax_indicator) %>%
  summarise(
    kraken_genus_taxid = unique(kraken_genus_taxid),
    tax_indicator = unique(tax_indicator),
    count = n()
      ) %>%
  group_by(kraken_genus_taxid) %>%
  filter(sum(count)>3) %>%
  filter(!is.na(tax_indicator)) %>%
  pivot_wider(names_from = tax_indicator, values_from = count) %>%
  replace(is.na(.), 0) %>%
  mutate(percent_agreement = Same/(Same + Different),
         tax_indicator = case_when(percent_agreement > 0.8 ~ "TP",
                                   percent_agreement <= 0.8 ~ "FP"))


## ----Merging the results to the original kraken report--------------------------------------------------------
kraken_report_indicator <- kraken_report %>%
  mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)) %>%
  left_join(genome_windows_sum, by = "kraken_genus_taxid") %>%
  left_join(tax_indicator, by = "kraken_genus_taxid") %>%
  filter(!is.na(Indicator)) %>%
  distinct(name, .keep_all = T) %>%
  mutate(tick_color = case_when((Indicator == "FP" | tax_indicator == "FP" ~ "red"),
                                (Indicator == "TP" | tax_indicator == "TP" ~ "darkgreen"),
                                (Indicator == "Small genome"~"yellow3")),
         tick_name = glue("<i style='color:{tick_color}'>{name}</i>")) 


## ----Creating the KrakenRefine Plot object--------------------------------------------------------------------
kraken_report_indicator$tick_name <- factor(kraken_report_indicator$tick_name, levels=(kraken_report_indicator$tick_name)[order(kraken_report_indicator$cladeReads)])
kraken_refine_plot<-ggplot(kraken_report_indicator, aes("",tick_name, fill = cladeReads)) +
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


## ----Creating the Normal Kraken results Plot object-----------------------------------------------------------
kraken_report_plot <- kraken_report %>%
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


## ----Creating the Kraken2Unique results Plot object-----------------------------------------------------------
kraken_unique_report_plot <- kraken_report %>%
  mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)) %>%
  mutate(Indicator = case_when((n_unique_kmers/n_kmers < 0.8 ~ "FP"),
                                (n_unique_kmers/n_kmers >= 0.8 ~ "TP"))) %>%
  mutate(tick_color = case_when((Indicator == "TP" ~"darkgreen"),
                                (Indicator == "FP"~"red")),
         tick_name = glue("<i style='color:{tick_color}'>{name}</i>"))

kraken_unique_report_plot$tick_name <- factor(kraken_unique_report_plot$tick_name, levels=(kraken_unique_report_plot$tick_name)[order(kraken_unique_report_plot$cladeReads)])

krakenuniqueplot<-ggplot(kraken_unique_report_plot, aes("",tick_name, fill = cladeReads)) +
  geom_tile() +
  theme_bw() +
  geom_text(aes(label = cladeReads), size = 3) +
  theme(axis.text.y = element_markdown(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_gradientn(colours = brewer.pal(8, "YlOrRd"),
                       values=rescale(c(0, 5, 250, 1000, 2500, 5000, 10000, 20000))) + 
  labs(x="", y = "Genus", fill = "Reads", title = "Grouped by Kraken2Uniq",
       tick_color = "Legend") +
  scale_color_manual(name = "KrakenRefine")


## ----Saving the image to a .svg file and coverage plots to multiple .svg files--------------------------------
if(nrow(kraken_report)>0) { #Checkpoint before printing plots since some datasets may not have >3 microbial reads
image<-plot_grid(krakennotrefineplot, kraken_refine_plot, krakenuniqueplot, labels = "AUTO")

title <- ggdraw() + 
  draw_label(
    sample,
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

image_title<-plot_grid(
  title, image,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave(file = paste0("data/KrakenRefine/", sample, "/KrakenRefine_", sample, ".svg"), plot = image_title, width = 20, height = 20)

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


for(i in seq(length(plot_df$plots))) {
  print(plot_df$plots[i])
  name<-plot_df$name[i]
  ggsave(file=paste0(WD, "/data/KrakenRefine/", sample, "/", sample, "_", name, ".svg"))
}
} else{
  print("No microbial genera w. > 3 reads")
}


## ----Write tsv files w. Refined report files------------------------------------------------------------------
write.table(kraken_report_indicator, file = paste0(WD, "/data/KrakenRefine/", sample, "/", sample, "_refined", ".tsv"), row.names = F, col.names = T, sep = "\t")

write.table(kraken_unique_report_plot, file = paste0(WD, "/data/KrakenRefine/", sample, "/", sample, "_uniqued", ".tsv"), row.names = F, col.names = T, sep = "\t")


## ----Species level classifications using bowtie2 alignments for the genera that were not flagged as false positive----
species_classification_df<-sam_files %>%
  mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)) %>%
  filter(kraken_genus_taxid %in% (kraken_report_indicator %>% filter(tick_color == "darkgreen") %>% select(kraken_genus_taxid) %>% mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)) %>% pull(kraken_genus_taxid))) %>%
     group_by(kraken_genus_taxid) %>%
  count(name) %>%
  left_join(nodes %>% rename(kraken_genus_taxid = id, genera = name) %>% mutate(kraken_genus_taxid = as.character(kraken_genus_taxid)), by="kraken_genus_taxid") %>%
  group_by(kraken_genus_taxid) %>%
  arrange(n) %>%
  top_n(15, n) %>%
  mutate(name=factor(name, levels=name)) 

species_classification<-species_classification_df%>%
  group_by(kraken_genus_taxid) %>%
  do(
  name = unique(as.character(.$genera)),
    plots = ggplot(data = ., aes(x = n, y = name, fill = name)) +
  geom_col() + 
  theme_bw() +
  labs(x = "Number of reads",
       y = "Species") +
    ggtitle(paste0("Top 15 species for ", unique(.$genera), ", id: ", unique(.$kraken_genus_taxid))) +
  theme(legend.position = "off")
)


for(i in seq(length(species_classification$plots))) {
  genera<-species_classification$name[i]
  ggsave(file=paste0(WD, "/data/KrakenRefine/", sample, "/", sample, "_", genera, ".svg"), plot = species_classification$plots[[i]], device = "svg")
}

write.table(species_classification_df, file = paste0(WD, "/data/KrakenRefine/", sample, "/", sample, "_species", ".tsv"), row.names = F, col.names = T, sep = "\t")

