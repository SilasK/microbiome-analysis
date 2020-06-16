import pandas as pd
import os


# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


#configfile: "config.yaml"
Comparisons= config['Comparisons'].keys()

rule all:
    input:
        'data.tsv',
        'metadata.tsv',
        expand("Comparisons/{comparison}/{file}",
               comparison=Comparisons,
               file=['stats_aldex.tsv','aldex_plot.pdf'] ), #'stats_relab.tsv'
        "Comparisons/cobined_stats_aldex.tsv",
        #"Comparisons/cobined_stats_relab.tsv"
        #"Correlations/rho.tsv"


localrules: init
rule init:
    input:
        config['data'],
        config['metadata']
    output:
        'data.tsv',
        'metadata.tsv'
    run:
        D= pd.read_table(input[0],index_col=0)

        metadata= pd.read_table(input[1],index_col=0)

        subset= D.index.intersection(metadata.index)

        metadata.loc[subset].to_csv(output[1],sep='\t')
        D=D.loc[subset]

        D= D.loc[:,D.mean()> config['min_count']]

        D.to_csv(output[0],sep='\t')


rule aldex:
    input:
        data='data.tsv',
        metadata='metadata.tsv'
    output:
        expand("Comparisons/{comparison}/stats_aldex.tsv",comparison=Comparisons )
    log:
        "logs/aldex2.txt"
    params:
        output_folder= "Comparisons"
    script:
        "scripts/Rscripts/aldex.R"


rule aldex_plot:
    input:
        "Comparisons/{comparison}/stats_aldex.tsv"
    output:
        "Comparisons/{comparison}/aldex_plot.pdf"
    run:
        import matplotlib
        import matplotlib.pylab as plt
        matplotlib.rcParams['pdf.fonttype']=42

        import sys
        sys.path.append(os.path.join(os.path.dirname(workflow.snakefile),'scripts'))
        from common import effect_plot as EP


        S= pd.read_table(input[0])
        EP.aldex_plot(S)
        plt.suptitle(wildcards.comparison)

        plt.savefig(output[0])


rule combine_stats:
    input:
        expand("Comparisons/{comparison}/stats_aldex.tsv",comparison=config['Comparisons'].keys() )
    output:
        "Comparisons/cobined_stats_aldex.tsv"
    params:
        comparisons = config['Comparisons'].keys()
    run:
        S={}
        for i, comparison in enumerate(params.comparisons):
            S[comparison] = pd.read_table(input[i])


        S=pd.concat(S,axis=1,sort=True)
        S.columns= S.columns.swaplevel()
        S.sort_index(axis=1,inplace=True)

        S.to_csv(output[0],sep='\t')






rule propr_rho:
    input:
        data='data.tsv',
        #metadata='metadata.tsv'
    output:
        "Correlations/rho.tsv"
    log:
        "logs/propr_rho.txt"
    params:
        output_folder= lambda wc,output: os.path.dirname(output[0])
    script:
        "scripts/Rscripts/correlations.R"



rule relab_analysis:
    input:
        data='data.tsv',
        metadata='metadata.tsv'
    output:
        expand("Comparisons/{comparison}/stats_relab.tsv",comparison=Comparisons )
    params:
        Comparisons= Comparisons,
        output_folder= "Comparisons"
    run:
        import helper_scripts as hs
        import numpy as np


        def calculate_pseudo_count(data,a=0.65):
            M= data.values
            return M[np.nonzero(M)].min()*a


        metadata= pd.read_table(input.metadata,index_col=0)
        data= pd.read_table(input.data,index_col=0)
        rel_data= (data.T/data.sum(1)).T
        pseudo_count= calculate_pseudo_count(rel_data)

        grouping_variable= metadata[config['grouping_variable']]

        for comparison_name in params.Comparisons:
            compared_groups = config['Comparisons'][comparison_name]


            V= hs.Viewpoint(rel_data,grouping_variable,Pairwise_Sig='kruskal',
                            Selection_rows=grouping_variable.isin(compared_groups),
                            order=compared_groups )


            S= V.Sumarize_Table(correction="Benjamini-Hochberg")
            S[('log2FC','log2FC')]= np.log2(S[(compared_groups[1]  ,'median')] +pseudo_count) - np.log2(S[(compared_groups[0]  ,'median')] +pseudo_count)

            S.to_csv(os.path.join(params.output_folder,comparison_name,'stats_relab.tsv'),sep='\t')



rule combine_relb_stats:
    input:
        expand("Comparisons/{comparison}/stats_relab.tsv",comparison=config['Comparisons'].keys() )
    output:
        "Comparisons/cobined_stats_relab.tsv"
    params:
        comparisons = config['Comparisons'].keys()
    run:
        S={}
        for i, comparison in enumerate(params.comparisons):
            d= pd.read_table(input[i],index_col=0,header=[0,1])
            d= d.loc[:,['median_diff','P_values','P Benjamini-Hochberg','log2FC']]
            d.columns= d.columns.droplevel(level=1)
            S[comparison] = d

        S=pd.concat(S,axis=1,sort=True)
        S.columns= S.columns.swaplevel()
        S.sort_index(axis=1,inplace=True)

        S.to_csv(output[0],sep='\t')
