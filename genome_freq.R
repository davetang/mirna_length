#frequency in hg38
a <- 898285419
c <- 623727342
g <- 626335137
t <- 900967885
acgt <- sum(a, c, g, t)
#set probabilities
my_prob <- c(a*100/acgt, c*100/acgt, g*100/acgt, t*100/acgt)
save(my_prob, file='data/hg38_nuc_freq.Robject')

#frequency in c. elegans
a <- 32371753
c <- 17782012
g <- 17759040
t <- 32373265
acgt <- sum(a, c, g, t)
#set probabilities
my_prob <- c(a*100/acgt, c*100/acgt, g*100/acgt, t*100/acgt)
save(my_prob, file='data/ce10_nuc_freq.Robject')

#frequency in zebrafish
a <- 446195643
c <- 258711649
g <- 258833640
t <- 446029177
acgt <- sum(a, c, g, t)
#set probabilities
my_prob <- c(a*100/acgt, c*100/acgt, g*100/acgt, t*100/acgt)
save(my_prob, file='data/danRer7_nuc_freq.Robject')

#frequency in mouse
a <- 773280124
c <- 552647817
g <- 552690118
t <- 774165441
acgt <- sum(a, c, g, t)
#set probabilities
my_prob <- c(a*100/acgt, c*100/acgt, g*100/acgt, t*100/acgt)
save(my_prob, file='data/mm10_nuc_freq.Robject')
