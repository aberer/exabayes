#! /usr/bin/Rscript

args<-commandArgs(TRUE)

tab = read.table(args[1], header=T)
name= args[2]
opt=args[3]

pdf(name)

if(opt == "line"){  
  plot(tab[,1], type="l", main=names(tab)[1])
} else{
  hist(tab[,1], breaks="FD", main=names(tab)[1])
}
bla = dev.off()
