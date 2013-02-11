#! /usr/bin/Rscript

args<-commandArgs(TRUE)

tab = read.table(args[1], header=F)
name= args[2]
opt=args[3]

pdf(name)

if(opt == "line"){  
  plot(tab[,1], type="l")
} else{
  hist(tab[,1], breaks="FD")
}
bla = dev.off()
