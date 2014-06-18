#!/bin/env Rscript

#store arguments in object
args <- commandArgs(TRUE)
#check arguments
args_length <- length(args)
if (args_length != 3){
   stop("Requires size of fragment as input, number to generate, and directory to write files")
}
my_size <- args[1]
my_number <- args[2]
my_dir <- args[3]

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

#set probabilities
my_prob <- c(0.25, 0.25, 0.25, 0.25) 

#generate a sequence
#sample(nucleotides, my_size, rep=T, prob=my_prob)

gen_seq <- function(x){
   my_seq <- sample(nucleotides, my_size, rep=T, prob=my_prob)
   my_seq <- paste(my_seq, sep='', collapse='')
   return(my_seq)
}

my_outfile <- paste(my_dir, '/my_random_seq_equal_', my_size, '_', my_number, '.fa', sep='')
my_random_seq <- sapply(1:my_number, gen_seq)
my_random_seq <- DNAStringSet(my_random_seq)
names(my_random_seq) <- 1:my_number
writeXStringSet(my_random_seq, my_outfile, format="fasta")
