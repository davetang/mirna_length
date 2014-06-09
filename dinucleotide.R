#install package if missing
required_package <- 'Biostrings'
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check){
   source("http://bioconductor.org/biocLite.R")
   biocLite("Biostrings")
}

#load package
library(Biostrings)

#store fasta
my_mirna <- readDNAStringSet('mature_thymine.fa.gz', 'fasta')

#get dinucleotide frequencies
di_freq <- dinucleotideFrequency(my_mirna)
colSums(di_freq)

#human miRNAs
human_mirna <- my_mirna[grep('Homo sapien',names(my_mirna))]
human_di_freq <- dinucleotideFrequency(human_mirna)
colSums(human_di_freq)

#probabilities of dinucleotides
colSums(human_di_freq)*100/sum(human_di_freq)

#frequency of first bases
colSums(alphabetFrequency(subseq(human_mirna,1,1))[,1:4])

#probabilites of first bases
colSums(alphabetFrequency(subseq(human_mirna,1,1))[,1:4])*100/length(human_mirna)
