#! /usr/bin/Rscript

library(lattice)


myVioBoxPlot = function(form,dta)
{
  bwplot(form , data=dta, panel = function(..., box.ratio) {
    panel.violin(..., col = "lightblue", varwidth = FALSE, box.ratio = box.ratio)
    panel.bwplot(..., col='black', cex=0.8, pch='|', fill='gray', box.ratio = .1)
  })
}


args =  commandArgs(trailingOnly=T)

supertab = c()
for (arg in args)
  {
    tab = read.table(arg, skip=1, header=T)
    tab = tail(tab, dim(tab)[1] * 3 / 4  ) 
    tab$name = arg
    supertab = rbind(supertab, tab)
  }

pdf("exa.pdf")

myVioBoxPlot(LnL ~ name , supertab)
myVioBoxPlot(TL ~ name , supertab)
myVioBoxPlot(r0AC ~ name , supertab)
myVioBoxPlot(r0AG ~ name , supertab)
myVioBoxPlot(r0AT ~ name , supertab)
myVioBoxPlot(r0CG ~ name , supertab)
myVioBoxPlot(r0CT ~ name , supertab)
myVioBoxPlot(r0GT ~ name , supertab)
myVioBoxPlot(pi0A ~ name , supertab)
myVioBoxPlot(pi0C ~ name , supertab)
myVioBoxPlot(pi0G ~ name , supertab)
myVioBoxPlot(pi0T ~ name , supertab)
myVioBoxPlot(alpha0 ~ name , supertab)

bla = dev.off()

