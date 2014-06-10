#!/bin/env Rscript

#store arguments in object
args <- commandArgs(TRUE)
#check for 1 arguments
args_length <- length(args)
if (args_length != 1){
   stop("Requires size of fragment as input")
}
my_size <- args[1]

#install package if missing
required_package <- 'Biostrings'
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check){
   source("http://bioconductor.org/biocLite.R")
   biocLite("Biostrings")
}

#load package
library(Biostrings)

#code adapted from http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter10.html
#generate sequence using Markov chain

#define the alphabet of nucleotides
nucleotides <- c('A', 'C', 'G', 'T')

#frequency in hg38
a <- 898285419
c <- 623727342
g <- 626335137
t <- 900967885
acgt <- sum(a, c, g, t)

#set probabilities
my_prob <- c(a*100/acgt, c*100/acgt, g*100/acgt, t*100/acgt)

#generate a sequence
#sample(nucleotides, my_size, rep=T, prob=my_prob)

gen_seq <- function(x){
   my_seq <- sample(nucleotides, my_size, rep=T, prob=my_prob)
   my_seq <- paste(my_seq, sep='', collapse='')
   return(my_seq)
}

my_number <- 1000000
my_outfile <- paste('my_random_seq_gen_freq_', my_size, '.fa', sep='')
my_random_seq <- sapply(1:my_number, gen_seq)
my_random_seq <- DNAStringSet(my_random_seq)
names(my_random_seq) <- 1:my_number
writeXStringSet(my_random_seq, my_outfile, format="fasta")
