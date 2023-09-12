library(tidyverse)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(data.table)
new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(total_rna_sig = case_when(padj < 0.1 & log2FoldChange > 0.75 ~ " Upregulated",
                                     padj < 0.1 & log2FoldChange < -0.75 ~ "Downregulated",
                                     T ~ "not_significant"))

new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(log2Fold_newRNA = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    mutate(gene = gsub("\\..*", "", gene))

uni = new_ratio_bayesian_p_de |>  pull(gene) |> unique()


sig_genes_df = new_ratio_bayesian_p_de |> 
    filter(new_rna_sig != 'unclear') |> 
    filter(total_rna_sig != 'not_significant') |> 
    dplyr::select(gene,time,new_rna_sig,total_rna_sig) |> 
    group_by(time,         new_rna_sig, total_rna_sig) |> 
    add_count() |> 
    ungroup() |> 
    filter(n >= 10) |> 
    dplyr::select(-n) |> 
    unique()


formula_res_all <- compareCluster(gene~time+total_rna_sig+new_rna_sig, 
                               data=sig_genes_df, 
                               fun="enrichGO",
                               OrgDb='org.Hs.eg.db',
                               keyType = "ENSEMBL",
                               universe = uni)

dotplot(formula_res_all) +
    facet_wrap(~time,scales = 'free_x')

k = pairwise_termsim(formula_res_all)

enrichplot::treeplot(k)

dotplot(formula_res_all)


formula_res_all@compareClusterResult |> 
    as.data.table() |> 
    separate_rows(geneID,sep = '/') |> 
    left_join(annotables::grch38 |> select(ensgene,symbol),by = c('geneID' = 'ensgene')) |> 
    as.data.table() |> 
    filter(grepl("MAP",Description))


matrix_genes = formula_res_all@compareClusterResult |> 
    as.data.table() |> 
    separate_rows(geneID,sep = '/') |> 
    left_join(annotables::grch38 |> select(ensgene,symbol),by = c('geneID' = 'ensgene')) |> 
    as.data.table() |> 
    filter(grepl("matrix",Description)) |> filter(time == 2) |> 
    pull(geneID)

formula_res_all@compareClusterResult |> 
    as.data.table() |> 
    separate_rows(geneID,sep = '/') |> 
    left_join(annotables::grch38 |> select(ensgene,symbol),by = c('geneID' = 'ensgene')) |> 
    as.data.table() |> 
    filter(grepl("ribo",Description)) |> filter(time == 2) 
         
new_ratio_bayesian_p_de |> 
    filter(time == 6) |> 
    filter(!is.na(mean_control_ntr)) |> 
    mutate(matrix_gene = gene %in% matrix_genes) |> 
    ggplot(aes(color = matrix_gene, x = mean_control_ntr)) + 
    stat_ecdf()
