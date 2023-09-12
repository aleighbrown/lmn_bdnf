library(tidyverse)
library(ggrepel)
library(clusterProfiler)
library(data.table)
new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(total_rna_sig = case_when(padj < 0.1 & log2FoldChange > 0.75 ~ " Upregulated",
                                     padj < 0.1 & log2FoldChange < -0.75 ~ "Downregulated",
                                     T ~ "not_significant"))

new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(log2Fold_newRNA = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    mutate(gene = gsub("\\..*", "", gene))



# 1 hour total rna volcano ------------------------------------------------


mapping = gprofiler2::gconvert(new_ratio_bayesian_p_de$gene,target = 'ENTREZGENE_ACC')
entrezed_table = new_ratio_bayesian_p_de |> 
    left_join(mapping,by = c("gene" = 'input')) |> 
    filter(!is.na(target)) |> 
    select(target,time,total_rna_sig) |> 
    unique()

sig_genes_df = new_ratio_bayesian_p_de |> 
    filter(total_rna_sig != "not_significant") |> 
    select(gene,time,total_rna_sig)

sig_genes_df_enz = entrezed_table |> 
    filter(total_rna_sig != "not_significant") |> 
    select(target,time,total_rna_sig)
# expressed_universe = new_ratio_bayesian_p_de |> 
#     filter(baseMean >=10) |> 
#     pull(gene) |> unique()
uni = new_ratio_bayesian_p_de |>  
    mutate(ens = gsub("\\..*", "",gene)) |> 
    pull(ens) |> unique()


formula_res <- compareCluster(gene~total_rna_sig+time, 
                              data=sig_genes_df, 
                              fun="enrichGO",
                              OrgDb='org.Hs.eg.db',
                              keyType = "ENSEMBL",
                              universe = uni)

# formula_res_enz <- compareCluster(target~total_rna_sig+time, 
#                               data=sig_genes_df_enz, 
#                               fun="enrichGO",
#                               OrgDb='org.Hs.eg.db',
#                               ont = 'ALL')
# 
# formula_res_nobg <- compareCluster(gene~total_rna_sig+time, 
#                               data=sig_genes_df, 
#                               fun="enrichGO",
#                               OrgDb='org.Hs.eg.db',
#                               keyType = "ENSEMBL")

 dotplot(formula_res, x="total_rna_sig") + 
    facet_grid(~time) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab(element_blank()) +
    theme(strip.text = element_text(size = 18)) +
    theme(axis.text.y  = element_text(size = 9)) 

# heatplot(formula_res, x="total_rna_sig") 


# 
# with_entrz2 = new_ratio_bayesian_p_de |> 
#     filter(time == 2) |> 
#     mutate(both_rna = glue::glue("{new_rna_sig}-{total_rna_sig}")) |> 
#     filter(both_rna != "unclear-not_significant") |> 
#     left_join(conversion_name, by = c("gene" = "input")) |> 
#     select(target,time,both_rna,new_rna_sig,total_rna_sig) |> 
#     unique() |> 
#     filter(!is.na(target))
# 
# formula_res2 <- compareCluster(target~total_rna_sig+new_rna_sig, data=with_entrz2, fun="enrichGO",OrgDb='org.Hs.eg.db')
# dotplot(formula_res2) +
#     facet_wrap(~total_rna_sig+new_rna_sig,scales = 'free') + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#     xlab(element_blank())
