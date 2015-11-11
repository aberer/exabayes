#! /usr/bin/Rscript


args = commandArgs(trailingOnly = T )

if(length(args) != 1 )
  {
    cat("./script file\n")
    quit()
  }

tab = read.table(args[1] , skip=1 , header=T)

start = dim(tab)[1] / 25
end =  dim(tab)[1] 


pdf(paste(args[1], ".pdf", sep=""))
for( i in 2:length(names(tab)))
{
  plot(sort(tab[start:end,i]), xlab="gen", ylab=names(tab)[i])
  boxplot(tab[start:end,i], ylab=names(tab)[i])
}
bla = dev.off()
