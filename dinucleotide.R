#install package if missing
required_package <- 'Biostrings'
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check){
   source("http://bioconductor.org/biocLite.R")
   biocLite("Biostrings")
}

#load package
library(Biostrings)

nucleotides <- c('A', 'C', 'G', 'T')

#store fasta
my_mirna <- readDNAStringSet('mature_thymine.fa', 'fasta')

#get dinucleotide frequencies
di_freq <- dinucleotideFrequency(my_mirna)
colSums(di_freq)

#human miRNAs
organism <- 'Homo sapien'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
di_freq <- dinucleotideFrequency(organism_mirna)
#probabilities of dinucleotides
di_freq_summary <- colSums(di_freq)*100/sum(di_freq)
mytransitionmatrix <- matrix(di_freq_summary, nrow=4, byrow=T)
rownames(mytransitionmatrix) <- nucleotides
colnames(mytransitionmatrix) <- nucleotides
#frequency of first bases
#colSums(alphabetFrequency(subseq(human_mirna,1,1))[,1:4])
#probabilites of first bases
myinitialprobs <- unname(colSums(alphabetFrequency(subseq(organism_mirna,1,1))[,1:4])*100/length(organism_mirna))
save(myinitialprobs, file='hg38_init_prob.Robject')
save(mytransitionmatrix, file='hg38_trans_mat.Robject')

#mouse miRNAs
organism <- 'Mus musculus'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
di_freq <- dinucleotideFrequency(organism_mirna)
#probabilities of dinucleotides
di_freq_summary <- colSums(di_freq)*100/sum(di_freq)
mytransitionmatrix <- matrix(di_freq_summary, nrow=4, byrow=T)
rownames(mytransitionmatrix) <- nucleotides
colnames(mytransitionmatrix) <- nucleotides
#frequency of first bases
#colSums(alphabetFrequency(subseq(human_mirna,1,1))[,1:4])
#probabilites of first bases
myinitialprobs <- unname(colSums(alphabetFrequency(subseq(organism_mirna,1,1))[,1:4])*100/length(organism_mirna))
save(myinitialprobs, file='mm10_init_prob.Robject')
save(mytransitionmatrix, file='mm10_trans_mat.Robject')

#celegans miRNAs
organism <- 'Caenorhabditis elegans'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
di_freq <- dinucleotideFrequency(organism_mirna)
#probabilities of dinucleotides
di_freq_summary <- colSums(di_freq)*100/sum(di_freq)
mytransitionmatrix <- matrix(di_freq_summary, nrow=4, byrow=T)
rownames(mytransitionmatrix) <- nucleotides
colnames(mytransitionmatrix) <- nucleotides
#frequency of first bases
#colSums(alphabetFrequency(subseq(human_mirna,1,1))[,1:4])
#probabilites of first bases
myinitialprobs <- unname(colSums(alphabetFrequency(subseq(organism_mirna,1,1))[,1:4])*100/length(organism_mirna))
save(myinitialprobs, file='ce10_init_prob.Robject')
save(mytransitionmatrix, file='ce10_trans_mat.Robject')

#zebrafish
organism <- 'Danio rerio'
organism_mirna <- my_mirna[grep(organism,names(my_mirna))]
di_freq <- dinucleotideFrequency(organism_mirna)
#probabilities of dinucleotides
di_freq_summary <- colSums(di_freq)*100/sum(di_freq)
mytransitionmatrix <- matrix(di_freq_summary, nrow=4, byrow=T)
rownames(mytransitionmatrix) <- nucleotides
colnames(mytransitionmatrix) <- nucleotides
#frequency of first bases
#colSums(alphabetFrequency(subseq(human_mirna,1,1))[,1:4])
#probabilites of first bases
myinitialprobs <- unname(colSums(alphabetFrequency(subseq(organism_mirna,1,1))[,1:4])*100/length(organism_mirna))
save(myinitialprobs, file='danRer7_init_prob.Robject')
save(mytransitionmatrix, file='danRer7_trans_mat.Robject')
