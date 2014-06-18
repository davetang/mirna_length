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
library(scales)

my_org <- c('hg38', 'mm10', 'danRer7', 'ce10')
for (org in my_org){
   my_file <- paste(org, '/', 'equal_perfect_mapped_1000000.txt', sep='')
   equal_perfect <- scan(my_file)
   my_file <- paste(org, '/', 'gen_freq_perfect_mapped_1000000.txt', sep='')
   genome_perfect <- scan(my_file)
   my_file <- paste(org, '/', 'markov_perfect_mapped_1000000.txt', sep='')
   markov_perfect <- scan(my_file)
   count_perfect <- data.frame(length=factor(15:30),
                               markov=markov_perfect+1,
                               genome=genome_perfect+1,
                               equal=equal_perfect+1
                              )
   count_perfect_melt <- melt(count_perfect, id.vars='length')
   #larger font
   theme_set(theme_gray(base_size = 20))

   my_outfile <- paste('image/', org, '_perfect_mapped_1000000.eps', sep='')

   p <- ggplot(count_perfect_melt, aes(x=length, y=value, fill=variable)) +
   geom_bar(position="dodge", stat="identity") +
   xlab('Length of random sequence') +
   ylab('Number of perfect matches') +
   scale_y_continuous(trans = log2_trans(),
                      labels=trans_format("log2", math_format(2^.x)),
                      breaks=mapply(function(x) 2^x, seq(0,20,2)))

   ggsave(p, file=my_outfile)

   my_file <- paste(org, '/', 'equal_mapped_1000000.txt', sep='')
   equal_mapped <- scan(my_file)
   my_file <- paste(org, '/', 'gen_freq_mapped_1000000.txt', sep='')
   genome_mapped <- scan(my_file)
   my_file <- paste(org, '/', 'markov_mapped_1000000.txt', sep='')
   markov_mapped <- scan(my_file)
   count_mapped <- data.frame(length=factor(15:30),
                              markov=markov_mapped+1,
                              genome=genome_mapped+1,
                              equal=equal_mapped+1
                             )
   count_mapped_melt <- melt(count_mapped, id.vars='length')

   my_outfile <- paste('image/', org, '_mapped_1000000.eps', sep='')

   p <- ggplot(count_mapped_melt, aes(x=length, y=value, fill=variable)) +
   geom_bar(position="dodge", stat="identity") +
   xlab('Length of random sequence') +
   ylab('Number mapped') +
   scale_y_continuous(trans = log2_trans(),
                      labels=trans_format("log2", math_format(2^.x)),
                      breaks=mapply(function(x) 2^x, seq(0,20,2)))

   ggsave(p, file=my_outfile)
}

file.remove('Rplots.pdf')
