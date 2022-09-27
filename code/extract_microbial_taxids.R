libpath<-snakemake@input[[2]]
.libPaths(c(libpath, .libPaths()))

suppressMessages(suppressWarnings(library(dplyr,warn.conflicts = F, quietly = T)))

krakendatapath<-snakemake@input[[1]]

read_kraken2_report<-function(data) {
    report<-read.csv2(data, header = F, comment.char = "#", sep = "\t", 
                      col.names = c("percentage", "cladeReads", 
                                    "taxonReads", "n_kmers", "n_unique_kmers",
                                    "taxRank", "taxID", "name"))  
return(report)
}

filter_genus_tax<-function(report_file) {
    genus_tax_list<-filter(report_file, taxRank == "G" 
    #& taxID != "9605"
    )$taxID
}



report<-read_kraken2_report(krakendatapath)

report<-filter_genus_tax(report)

write.table(report, file=snakemake@output[[1]], 
            row.names=FALSE, col.names=FALSE, sep="")
