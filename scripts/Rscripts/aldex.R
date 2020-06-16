
library(ALDEx2)

sink(snakemake@log[[1]])



metadata= read.delim(snakemake@input[['metadata']], stringsAsFactors = F,row.names = 1)
data= read.delim(file=snakemake@input[['data']],row.names = 1)

grouping_variable= metadata[,snakemake@config[['grouping_variable']]]
Comparisons= snakemake@config[['Comparisons']]


for (comparison in names(Comparisons))
{
  print(comparison)

  selected= grouping_variable %in% Comparisons[[comparison]]

  conditions = grouping_variable[selected]

  x <- aldex.clr(reads= t(data[selected,]), conds= conditions, mc.samples=128,denom='all',useMC=TRUE)

  d.eff <- aldex.effect(x,useMC = T)



  d.tt <- aldex.ttest(x)
  res.all <- data.frame(d.eff,d.tt)
  # get 'significant' set

  res.all$label = rep(0,nrow(res.all))
  res.all[res.all$wi.eBH < 0.1,"label"] = "Significant"
  res.all[which(res.all$label == 0), "label"] = "Not significant"


  write.table(res.all,paste0(snakemake@params$output_folder,'/',comparison,'/stats_aldex.tsv')   ,sep='\t')

  write.table(res.all[res.all$wi.eBH < 0.1,],paste0(snakemake@params$output_folder,'/',comparison,'/sig.tsv')   ,sep='\t')


}



  d.tt <- aldex.ttest(x, conditions)
  res.all <- data.frame(d.eff,d.tt)
  # get 'significant' set

  res.all$label = rep(0,nrow(res.all))
  res.all[res.all$wi.eBH < 0.1,"label"] = "Significant"
  res.all[which(res.all$label == 0), "label"] = "Not significant"


  write.table(res.all,paste0(snakemake@params$output_folder,'/',comparison,'/stats_aldex.tsv')   ,sep='\t')

  write.table(res.all[res.all$wi.eBH < 0.1,],paste0(snakemake@params$output_folder,'/',comparison,'/sig.tsv')   ,sep='\t')


}


#
