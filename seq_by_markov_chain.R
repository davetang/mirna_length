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

#human dinucleotide probabilities
#AA       AC       AG       AT       CA       CC       CG       CT
#4.945086 4.534408 8.181527 4.506151 6.518094 6.606635 2.505510 7.970537
#GA       GC       GG       GT       TA       TC       TG       TT
#6.503024 6.418251 9.624550 6.260008 3.867528 5.973664 9.230827 6.354200

#probabilites of first bases
#A        C        G        T
#27.92863 21.83863 13.84794 36.38479

#code adapted from http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/src/chapter10.html
#generate sequence using Markov chain

#define the alphabet of nucleotides
nucleotides <- c('A', 'C', 'G', 'T')

#set the values of the probabilities, where the previous nucleotide was "A"
afterAprobs <- c(4.945086, 4.534408, 8.181527, 4.506151)

# Set the values of the probabilities, where the previous nucleotide was "C"
afterCprobs <- c(6.518094, 6.606635, 2.505510, 7.970537)

# Set the values of the probabilities, where the previous nucleotide was "G"
afterGprobs <- c(6.503024, 6.418251, 9.624550, 6.260008)

# Set the values of the probabilities, where the previous nucleotide was "T"
afterTprobs <- c(3.867528, 5.973664, 9.230827, 6.354200)

# Create a 4 x 4 matrix
mytransitionmatrix <- matrix(c(afterAprobs, afterCprobs, afterGprobs, afterTprobs), 4, 4, byrow = TRUE)
rownames(mytransitionmatrix) <- nucleotides
colnames(mytransitionmatrix) <- nucleotides
# Print out the transition matrix
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

myinitialprobs <- c(27.92863, 21.83863, 13.84794, 36.38479)

my_number <- 1000000

my_outfile <- paste('my_random_seq_', my_size, '.fa', sep='')

my_random_seq <- sapply(1:my_number, generatemarkovseq, transitionmatrix=mytransitionmatrix, initialprobs=myinitialprobs, seqlength=my_size)

my_random_seq <- DNAStringSet(my_random_seq)
names(my_random_seq) <- 1:my_number
writeXStringSet(my_random_seq, my_outfile, format="fasta")
