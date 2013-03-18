#! /usr/bin/Rscript


options <- commandArgs(trailingOnly = TRUE) 
                                            
tab.1 = read.table(options[1], header=T)   

pdf(paste(options[1],".pdf",sep=""))

plot(log10(sort(tab.1[,"alpha"])), main="alpha", xlab="index", ylab="log10(value)")

tmp = "ac"
plot(log10(sort(tab.1[,tmp])), main=tmp, xlab="index", ylab="log10(value)")

tmp = "ag"
plot(log10(sort(tab.1[,tmp])), main=tmp, xlab="index", ylab="log10(value)")

tmp = "at"
plot(log10(sort(tab.1[,tmp])), main=tmp, xlab="index", ylab="log10(value)")

tmp = "cg"
plot(log10(sort(tab.1[,tmp])), main=tmp, xlab="index", ylab="log10(value)")

tmp = "ct"
plot(log10(sort(tab.1[,tmp])), main=tmp, xlab="index", ylab="log10(value)")

bla = dev.off()

