#install package if missing
required_package <- c('Biostrings', 'ggplot2')
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check[1]){
   source("http://bioconductor.org/biocLite.R")
   biocLite("Biostrings")
}
if(!my_check[2]){
   install.packages("ggplot2", repos='http://cran.us.r-project.org')
}

#load package
library(Biostrings)
library(ggplot2)

#store fasta
my_mirna <- readDNAStringSet('mature_thymine.fa', 'fasta')

#human miRNAs
organism <- 'Homo sapien'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
mirna_length <- factor(width(organism_mirna))
postscript('image/human_mirbase_length.eps')
qplot(mirna_length, xlab='miRNA length', ylab='Frequency')
dev.off()

#mouse miRNAs
organism <- 'Mus musculus'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
mirna_length <- factor(width(organism_mirna))
postscript('image/mouse_mirbase_length.eps')
qplot(mirna_length, xlab='miRNA length', ylab='Frequency')
dev.off()

#celegans miRNAs
organism <- 'Caenorhabditis elegans'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
mirna_length <- factor(width(organism_mirna))
postscript('image/celegans_mirbase_length.eps')
qplot(mirna_length, xlab='miRNA length', ylab='Frequency')
dev.off()

#zebrafish
organism <- 'Danio rerio'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
mirna_length <- factor(width(organism_mirna))
postscript('image/zebrafish_mirbase_length.eps')
qplot(mirna_length, xlab='miRNA length', ylab='Frequency')
dev.off()
