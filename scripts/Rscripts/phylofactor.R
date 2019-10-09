
# install

# conda create phylofactor r rstudio r-devtools bioconductor-ggtree bioconductor-biostrings

for (package in c("ggtree","Biostrings"))#,"zCompositions"))
{
  if (!(package %in% rownames(installed.packages())))
  {
    BiocManager::install(package)
  }
}
# devtools::install_github('reptalex/phylofactor')
#setwd('~/Desktop/CRCEmicrobiota/CE/16S_CE/Analysis/')
setwd('~/Desktop/CRCEmicrobiota/CR/16S_CR/Analysis/')



# setwd('~/Desktop/WarmMicrobiota/Sequencing_C/Analysis/Analysis/')

library(Biostrings)
library(ggtree)
library(phylofactor)
library(zCompositions)
#library(phangorn)

tree= read.tree("../taxonomy/otu_tree.nwk")

#plot(tree)

D= read.delim("./data.tsv",row.names = 1)
metadata = read.delim('metadata.tsv',row.names = 1)
#taxonomy= read.delim("../taxonomy/Silva.tsv",row.names = 1)
taxonomy= read.delim("taxonomy_gg.tsv")
row.names(taxonomy)= taxonomy$X



subset= row.names(metadata)#[metadata$Source =='Cecum',])  # [metadata$Time==6,]) #[(metadata$Group=='SH-RT')|(metadata$Group=='SH-H') ,])#
data= D[subset,]
metadata= metadata[subset,]

# it is important to first filter to remove rows that are exclusively 0 values
data <- data[rowSums(data) > 0,]

# we are using the Count Zero Multiplicative approach
data <- cmultRepl(t(data), method="CZM", label=0)

# generate the centered log-ratio transformed data
# samples by row

d.clr= apply(data, 2, function(x) log(x) - mean(log(x)))


pf= PhyloFactor(as.matrix(data),tree,X= metadata,frmla = Data~ Group, ncores = 3,nfactors = 10)


barplot(pf$ExplainedVar[1:50])
summary(pf)
#pf$factors


pf.tree(pf,layout='rectangular',factors = 3)$ggplot

species.groups <- pf.groupsTospecies(pf)
pf.plot(pf,factors = 3)

fofi=3

pf.taxa(pf,taxonomy,factor = fofi)$group1
cluster= species.groups[[fofi]][[1]]
#taxonomy[cluster, 2]

cluster

FactorSummary=pf.summary(pf,taxonomy,factor=fofi)

library(beeswarm)



plot_data=data.frame( ILR=FactorSummary$ilr, group=metadata$Group)
beeswarm(ILR~  group, data =plot_data)
bxplot(ILR~  group, data =plot_data, add = TRUE)



### plot with summary
par(mfrow=c(1,2))
plot(FactorSummary$ilr,ylab='ILR coordinate',main='ILR coordinate of factor',
     xlab='sample no.',pch=16)
lines(FactorSummary$fitted.values,lwd=2,col='blue')
legend(x=1,y=-5,list('data','prediction'),pch=c(16,NA),lty=c(NA,1),
       col=c('black','blue'),lwd=c(NA,2))

plot(FactorSummary$MeanRatio,ylab='ILR coordinate',main='Mean Ratio of Grp1/Grp2',
     xlab='sample no.',pch=16)
lines(FactorSummary$fittedMeanRatio,lwd=2,col='blue')
legend(x=1,y=-5,list('data','prediction'),pch=c(16,NA),lty=c(NA,1),
       col=c('black','blue'),lwd=c(NA,2))



###


FactorSummary$ilr


plot(colSums(data[cluster,]))
plot(colMeans(d.clr[cluster,]))

#Gastranaerophilales decreased


pf.heatmap(tree=pf$tree,Data=  apply(pf$Data, 2, function(x) log(x) - mean(log(x))))



