args <- commandArgs(T)
data <- read.table(args[1],header=F)
p <- hist(data$V1,breaks = 100,col='blue',main='Subreads Length Distribution',xlab='Read Length',ylab='Read Number')
pdf('read_distribution_hist.pdf')
png('read_distribution_hist.png')
dev.off()

