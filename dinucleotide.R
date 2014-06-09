#install package if missing
required_package <- 'Biostrings'
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check) source("http://bioconductor.org/biocLite.R"); biocLite("Biostrings")

#load package
library(Biostrings)

#store fasta
my_mirna <- readDNAStringSet('mature_thymine.fa.gz', 'fasta')

#get dinucleotide frequencies
di_freq <- dinucleotideFrequency(my_mirna)
colSums(di_freq)

#human miRNAs
my_mirna[grep('Homo sapien',names(my_mirna))]
human_di_freq <- dinucleotideFrequency(my_mirna[grep('Homo sapien',names(my_mirna))])
colSums(human_di_freq)
