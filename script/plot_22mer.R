exon <- c(4.15,4.17,4.21,4.15,4.11,3.59,3.57,3.60,3.57,3.61,4.45,4.50,4.48,4.51,4.45)
intron <- c(53.01,52.97,53.12,53.22,53.20,51.97,52.06,52.04,51.96,51.96,53.89,53.80,53.91,53.77,53.83)
intergenic <- c(42.83,42.86,42.67,42.63,42.69,44.43,44.36,44.36,44.47,44.43,41.66,41.69,41.62,41.72,41.72)
my_type <- c('Equal', 'Genome', 'Markov')
type <- factor(rep(my_type, each=5), levels=my_type)

df <- data.frame(exon=exon, intron=intron, intergenic=intergenic, type=type)

sem <- function(x){
   sd(x)/sqrt(length(x))
}

library(ggplot2)

#larger font
theme_set(theme_gray(base_size = 20))

for (i in c('exon', 'intron', 'intergenic')){
   my_mean <- aggregate(subset(df, select=colnames(df)==i), by=list(type), mean)
   my_sem  <- aggregate(subset(df, select=colnames(df)==i), by=list(type), sem)
   mean_sem <- data.frame(mean=my_mean[,2], sem=my_sem[,2], group=my_type)

   my_outfile <- paste('image/', '22mer_', i, '.eps', sep='', collapse='')
   #plot using ggplot
   p <- ggplot(mean_sem, aes(x=group, y=mean)) +
     geom_bar(stat='identity') +
     geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                   width=.2) +
     ylab('Percent mapped')
   ggsave(p, file=my_outfile)
}

file.remove('Rplots.pdf')
