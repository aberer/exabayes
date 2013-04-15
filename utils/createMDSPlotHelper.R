#! /usr/bin/Rscript


args = commandArgs(trailingOnly = TRUE)


fn = args[1]                          #rf distances
numTreeFile = args[2]                      # number of trees per run 
name = args[3]                          # an id ofr the run

tab = read.table(fn, header=F)

numTrees = read.table(numTreeFile, header=F)

maximum = max(tab[,2])

x = matrix(0,nrow=maximum+1, ncol=maximum+1 )
x[upper.tri(x)] = tab[,3]
x[lower.tri(x)] = t(x[upper.tri(x)])

rfDists = x

## rfDists
scaled = cmdscale(rfDists)
## scaled
## stopifnot(dim(scaled)[1] == sum(numTrees)-1 )

pdf(paste("mdsplot-", name, ".pdf", sep=""))
matplot(NA,NA, xlim=range(scaled[,1]),ylim=range(scaled[,2]), xlab="dim A", ylab="dim B")
end = 0 
for (i in 1:dim(numTrees)[1])
  {
    start = end + 1 
    end = start + numTrees[i,1] - 1 

    ## if(i == dim(numTrees)[1])
    ##   {
    ##     end = end - 1 
    ##   }

    matlines(scaled[start:end,1],scaled[start:end,2], col=i)
  }

bla = dev.off()

