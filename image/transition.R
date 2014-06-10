#install package if missing
required_package <- 'diagram'
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check){
  install.packages("diagram")
}

#load package
library(diagram)

nucleotides <- c('A', 'C', 'G', 'T')

# Set the values of the probabilities, where the previous nucleotide was "A"
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

#for plotmat the rows specific to and columns from
#so we need to transpose our matrix
mytransitionmatrix <- t(round(mytransitionmatrix,2))

plotmat(mytransitionmatrix, #transition matrix
        relsize=0.7,
        name = nucleotides, #names of the states
        box.lwd=1, #outline of state
        cex.txt=0.5, #size of probabilities
        box.prop=1, #size of box
        box.type = 'circle',
        self.cex = 0.5, #size of self probability
        lwd = 1, #outline of probabilities
        box.cex=2 #size of text in box
)