args <- commandArgs(T)
data <- read.table(args[1],header=F)
prefix <- args[2]
pdf(paste(prefix,'.read_distribution_hist.pdf',sep=""))
a<-dev.cur()
png(paste(prefix,'.read_distribution_hist.png',sep=""))
dev.control("enable")
p <- hist(data$V1,breaks = 100,col='blue',main='Subreads Length Distribution',xlab='Read Length',ylab='Read Number')
p
dev.copy(which=a)
dev.off()
dev.off()

