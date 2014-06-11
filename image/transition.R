#!/bin/env Rscript

#store arguments in object
args <- commandArgs(TRUE)
#check for 1 argument
args_length <- length(args)
if (args_length != 1){
   stop("Requires transition matrix object")
}
my_object <- args[1]

load(my_object)

#install package if missing
required_package <- 'diagram'
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check){
   install.packages("diagram", repos='http://cran.us.r-project.org')
}

#load package
library(diagram)

nucleotides <- c('A', 'C', 'G', 'T')

#for plotmat the rows specific to and columns from
#so we need to transpose our matrix
mytransitionmatrix <- t(round(mytransitionmatrix,2))

my_outfile <- my_object
my_outfile <- gsub("_.*", '.eps', my_outfile, perl=T)
my_outfile <- paste('image/', my_outfile, sep='', collapse='')

postscript(my_outfile)
plotmat(mytransitionmatrix, #transition matrix
        relsize=0.7,
        name = nucleotides, #names of the states
        box.lwd=1, #outline of state
        cex.txt=1.2, #size of probabilities
        box.prop=1, #size of box
        box.type = 'circle',
        self.cex = 1.2, #size of self probability
        lwd = 1, #outline of probabilities
        box.cex=2 #size of text in box
		  )
dev.off()

quit()
