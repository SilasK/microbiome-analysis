library(propr)


output_folder= snakemake@params[['output_folder']]

#metadata= read.delim(snakemake@input[['metadata']], stringsAsFactors = F,row.names = 1)
data= read.delim(file=snakemake@input[['data']],row.names = 1)



rho = perb(data)

write.table(rho@matrix,paste0(output_folder,"/rho.tsv"),sep = '\t')

pdf(file = paste0(output_folder,"/rho_dist.pdf"))
hist(rho@matrix)
dev.off()


rho= updateCutoffs(rho,cutoff= c(0.65,0.75,0.85,0.95)) ; 



write.table(rho@fdr,paste0(output_folder,"/rho_fdr.tsv"),sep = '\t')


#choose lowest cutoff, which has a fdr below choosen value
cutoff <- min(rho@fdr[rho@fdr$FDR<  snakemake@config[['fdr_correlations']] ,'cutoff'])

# get the subset of OTUs that are joined by one or more low phi connections

rho=rho[">",cutoff]
rho=rho["<",- cutoff]

print(paste0("Found ",length(rho@pairs)," interactions with cutoff ", cutoff))

pdf(file = paste0(output_folder,"/Graph.pdf"))
graph= cytescape(rho)
dev.off()



write.table(graph,paste0(output_folder,"/Graph.tsv"),sep = '\t')

