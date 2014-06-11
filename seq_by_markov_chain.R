#!/bin/env Rscript

#store arguments in object
args <- commandArgs(TRUE)
#check arguments
args_length <- length(args)
if (args_length != 5){
   stop("Requires size of fragment, transition matrix object, initial probabilities object, number to generate, and directory to write to")
}
my_size <- args[1]
mytransitionmatrix <- args[2]
myinitialprobs <- args[3]
my_number <- args[4]
my_dir <- args[5]

load(mytransitionmatrix)
load(myinitialprobs)

#install package if missing
required_package <- 'Biostrings'
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check){
   source("http://bioconductor.org/biocLite.R")
   biocLite("Biostrings")
}

#load package
library(Biostrings)

mytransitionmatrix

generatemarkovseq <- function(x, transitionmatrix, initialprobs, seqlength){
   nucleotides     <- c("A", "C", "G", "T") # Define the alphabet of nucleotides
   mysequence      <- character()           # Create a vector for storing the new sequence
   # Choose the nucleotide for the first position in the sequence:
   firstnucleotide <- sample(nucleotides, 1, rep=TRUE, prob=initialprobs)
   mysequence[1]   <- firstnucleotide       # Store the nucleotide for the first position of the sequence
   for (i in 2:seqlength){
      prevnucleotide <- mysequence[i-1]     # Get the previous nucleotide in the new sequence
      # Get the probabilities of the current nucleotide, given previous nucleotide "prevnucleotide":
      probabilities  <- transitionmatrix[prevnucleotide,]
      # Choose the nucleotide at the current position of the sequence:
      nucleotide     <- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
      mysequence[i]  <- nucleotide          # Store the nucleotide for the current position of the sequence
   }
	mysequence <- paste(mysequence, sep='', collapse='')
	return(mysequence)
}

my_outfile <- paste(my_dir, '/my_random_seq_markov_', my_size, '_', my_number, '.fa', sep='')

my_random_seq <- sapply(1:my_number, generatemarkovseq, transitionmatrix=mytransitionmatrix, initialprobs=myinitialprobs, seqlength=my_size)

my_random_seq <- DNAStringSet(my_random_seq)
names(my_random_seq) <- 1:my_number
writeXStringSet(my_random_seq, my_outfile, format="fasta")
