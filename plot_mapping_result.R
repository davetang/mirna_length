#install package if missing
required_package <- c('ggplot2', 'reshape2')
my_check <- required_package %in% installed.packages()[,"Package"]
if(!my_check[1]){
   install.packages("ggplot2", repos='http://cran.us.r-project.org')
}
if(!my_check[2]){
   install.packages("reshape2", repos='http://cran.us.r-project.org')
}

#load package
library(ggplot2)
library(reshape2)

#equal frequency perfect match
equal_perfect <- c(622903, 371803, 161709, 53697, 15071, 4020, 978, 272, 76, 14, 8, 1, 1)

#genome frequency perfect match
genome_perfect <- c(736400, 508275, 260889, 100013, 32196, 9423, 2589, 716, 193, 56, 12, 3, 0)

#markov chain perfect match
markov_perfect <- c(853182, 623820, 310149, 110331, 32937, 9104, 2601, 689, 170, 43, 8, 7, 0)

count_perfect <- data.frame(length=factor(15:27),
                            markov=markov_perfect,
                            genome=genome_perfect,
                            equal=equal_perfect
                            )

count_perfect_melt <- melt(count_perfect, id.vars='length')

#larger font
theme_set(theme_gray(base_size = 20))

ggplot(count_perfect_melt, aes(x=length, y=value, fill=variable)) +
geom_bar(position="dodge", stat="identity") +
xlab('Length of random sequence') +
ylab('Number of perfect matches')

#equal frequency mapped
equal_mapped <- c(999190, 999999, 999982, 998992, 984139, 908715, 723867, 456454, 216818, 80016, 25271, 7290, 2098, 604, 144, 46)

#genome frequency mapped
genome_mapped <- c(999575, 999997, 999989, 999608, 992752, 954001, 837775, 619178, 358830, 162785, 61092, 20422, 6339, 1883, 564, 176)

#markov chain mapped
markov_mapped <- c(999965, 999996, 999996, 999979, 999094, 988158, 924670, 727549, 415376, 171260, 58484, 18194, 5537, 1626, 450, 146)

count_mapped <- data.frame(length=factor(15:30),
                            markov=markov_mapped,
                            genome=genome_mapped,
                            equal=equal_mapped
                           )

count_mapped_melt <- melt(count_mapped, id.vars='length')

ggplot(count_mapped_melt, aes(x=length, y=value, fill=variable)) +
  geom_bar(position="dodge", stat="identity") +
  xlab('Length of random sequence') +
  ylab('Number mapped')
