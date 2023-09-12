library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)
library(data.table)
library(msigdbr)
#doi: https://doi.org/10.1101/2022.03.25.485808  from this
motor_neuron_markers = c("MT1X", "STMN2", "MAP1B", "NEFH", "NEFM", "PRPH", "KIF21A", 
                         "PRUNE2", "UTS2", "ATP1B1", "TMSB10", "RTN1", "RTN4", "RTN3", 
                         "TUBA4A", "ACLY", "SLC5A7", "SOD1", "GAPDH", "ACTB", "NEFL", 
                         "UCHL1 ", "KLC1", "S100B", "PVALB", "PLP1", "SNCG", "TUBB4B", 
                         "DYNC1H1", "YWHAG", "TUBB2A", "TUBA1B", "TMSB4X", "PKM", "HSP90AB1", 
                         "ACTG1", "TUBA1A", "UBB", "CALM1", "ANXA2", "HSP90AA1", "HSPA8", 
                         "AHNAK2", "CLU", "SPARCL1", "SPP1", "HSPB1", "S100A10", "LGALS1", 
                         "MCAM")

new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(total_rna_sig = case_when(padj < 0.1 & log2FoldChange > 0.75 ~ " upregulated",
                                     padj < 0.1 & log2FoldChange < -0.75 ~ "downregulated",
                                     T ~ "not_significant"))

new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(log2Fold_newRNA = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    mutate(gene = gsub("\\..*", "", gene))

# Category plot log2Fold 0.75 free ----------------------------------------

new_ratio_bayesian_p_de %>% 
    dplyr::select(gene,time,new_rna_sig,total_rna_sig) %>% 
    unique() %>% 
    filter(new_rna_sig != "unclear") |>
    filter(total_rna_sig != "not_significant") |>
    janitor::tabyl(total_rna_sig,new_rna_sig,time) %>% 
    melt() |> 
    ggplot(aes(x = total_rna_sig,
               y = value,
               fill = variable)) + 
    geom_col() +
    coord_flip() +
    ggpubr::theme_pubr() +
    facet_wrap(~(L1),nrow = 3,scales = 'free') +
    # theme(
    #     strip.background = element_blank(),
    #     strip.text.x = element_blank()
    # ) +
    # ggtitle('1hr BDNF Treatment') + 
    theme(legend.position = 'top') + 
    ylab("N genes") + 
    # scale_fill_manual(name="",
    #                     breaks=c("bdnf_equals_control", "bdnf_higher_new_rna", "bdnf_lower_new_rna"),
    #                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A"),
    #                     labels=c("Equivalent new RNA", "Higher new RNA", "Lower new RNA")) + 
    scale_x_discrete(labels=c("downregulated" = "Downregulated", "upregulated" = "Upregulated")) + 
    xlab(element_blank())




# Category plot log2Fold 0.75 fixed ----------------------------------------

new_ratio_bayesian_p_de %>% 
    dplyr::select(gene,time,new_rna_sig,total_rna_sig) %>% 
    unique() %>% 
    # filter(new_rna_sig != "unclear") |>
    filter(total_rna_sig != "not_significant") |>
    group_by(total_rna_sig,new_rna_sig,time) %>% 
    summarise(n = n_distinct(gene)) 
    ggplot(aes(x = total_rna_sig,
               y = n,
               fill = new_rna_sig)) + 
    geom_col() +
    coord_flip() +
    ggpubr::theme_pubr() +
    facet_wrap(~time,scales = 'free',nrow = 3) +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
    # ggtitle('1hr BDNF Treatment') + 
    theme(legend.position = 'top') + 
    ylab("N genes")  +
    scale_fill_manual(name="",
                      breaks=c("bdnf_equals_control", "bdnf_higher_new_rna", "bdnf_lower_new_rna",'unclear'),
                      values = c("#A6CEE3", "#1F78B4", "#B2DF8A","#666666"),
                      labels=c("Equivalent new RNA", "Higher new RNA", "Lower new RNA", "Not assessed")) 
    # scale_x_discrete(labels=c("downregulated" = "Downregulated", "upregulated" = "Upregulated")) + 
    # xlab(element_blank())
    # 



new_ratio_bayesian_p_de %>% 
    filter(time == 1) %>% 
    dplyr::select(gene,gene_name,time,new_rna_sig,total_rna_sig) %>% 
    unique() %>% 
    mutate(is_mn_marker = (gene_name %in% motor_neuron_markers))

new_ratio_bayesian_p_de %>% 
    filter(time == 2) %>% 
    filter(new_rna_sig != 'unclear') %>%
    select(gene,time,new_rna_sig,total_rna_sig) %>% 
    unique() %>% 
    janitor::tabyl(total_rna_sig,new_rna_sig) %>% 
    melt() %>% 
    ggplot(aes(x = 1,
               y = value,
               fill = variable)) + 
    geom_col() +
    coord_flip() +
    theme_minimal() +
    # theme(
    #     strip.background = element_blank(),
    #     strip.text.x = element_blank()
    # ) +
    scale_fill_brewer(palette = "Paired") + 
    ggtitle('6hr BDNF Treatment') + 
    theme(legend.position = 'top') + 
    ylab("N genes")

new_ratio_bayesian_p_de %>% 
    filter(time == 2) %>% 
    filter(new_rna_sig != 'unclear') %>%
    select(gene,time,new_rna_sig,total_rna_sig) %>% 
    unique() %>% 
    janitor::tabyl(total_rna_sig,new_rna_sig) %>% 
    melt() %>% 
    ggplot(aes(x = 1,
               y = value,
               fill = variable)) + 
    geom_col() +
    facet_wrap(~total_rna_sig,scales = 'free',nrow = 3,) +
    coord_flip() +
    theme_minimal() +
    # theme(
    #     strip.background = element_blank(),
    #     strip.text.x = element_blank()
    # ) +
    scale_fill_brewer(palette = "Paired") + 
    ggtitle('2hr BDNF Treatment') + 
    theme(legend.position = 'top') + 
    ylab("N genes")



new_ratio_bayesian_p_de %>% 
    filter(gene_name %in% c("ARC","EGR1","EGR2",
                            "EGR3","FOS","IER2","JUNB")) %>% 
    ggplot(aes(x = time, y = log2Fold_newRNA,
               group = gene_name,fill = gene_name)) + 
    geom_line() + 
    geom_point(size = 4,pch = 21) + 
    geom_hline(yintercept = 0) +
    theme_minimal() +
    ylab("Log2Fold Change") +
    scale_x_continuous(breaks = c(0,1,6)) 
    



new_ratio_bayesian_p_de |> 
    filter(!is.na(mean_bdnf_ntr) & !is.na(mean_control_ntr)) |> 
    mutate(confidentally_assessable = n_samp_passing_bdnf >=2 & n_samp_passing_control >=2) |> 
    mutate(confidentally_assessable = ifelse(is.na(confidentally_assessable),FALSE,confidentally_assessable)) |> 
    dplyr::select(gene,time,new_rna_sig,confidentally_assessable) |> 
    unique() |> 
    group_by(time,new_rna_sig,confidentally_assessable) |> 
    summarize(n = n_distinct(gene)) |> 
    ggplot(aes(x = as.character(time),
               y = n,
               fill = confidentally_assessable)) + 
    geom_col() + 
    labs(fill = 'Did 2 / 3 samples meet width criteria?')





new_ratio_bayesian_p_de |> 
    # mutate(new_rna_sig = case_when(bayesian_p > 0.9   ~ "bdnf_lower_new_rna",
    #                                bayesian_p < 0.1  ~ "bdnf_higher_new_rna",
    #                                bayesian_p > 0.4 & bayesian_p < 0.6  ~ "bdnf_equal_controls"
    # )) %>% 
    # filter( n_samp_passing_bdnf >=2 & n_samp_passing_control >=2) |> 
    filter(new_rna_sig != 'unclear') |> 
    filter(!is.na(new_rna_sig)) |> 
    dplyr::select(gene,time,new_rna_sig) |> 
    unique() |> 
    group_by(time,new_rna_sig) |> 
    summarize(n = n_distinct(gene)) |> 
    ggplot(aes(x = as.character(time),
               y = n,
               fill = new_rna_sig)) + 
    geom_col()


new_ratio_bayesian_p_de |> 
    ggplot(aes(x = log2FoldChange, y = log2Fold_newRNA)) + 
    geom_hex() + 
    facet_wrap(~time) +
    ggpubr::theme_pubr()

new_ratio_bayesian_p_de |> 
    ggplot(aes(x = log2FoldChange)) + 
    geom_density(size = 1.4) + 
    facet_wrap(~time,ncol = 1) + 
    ggpubr::theme_pubr() +
    geom_vline(linetype = 'dashed',xintercept = 0) 


new_ratio_bayesian_p_de |> 
    dplyr::select(gene,time, mean_bdnf_ntr,mean_control_ntr) |> 
    filter(!is.na(mean_bdnf_ntr),!is.na(mean_control_ntr)) |> 
    melt(id.vars = c("gene","time")) |> 
    ggplot(aes(x = value, color = variable)) + 
    geom_density(size = 1.4) + 
    facet_wrap(~time,ncol = 1) + 
    ggpubr::theme_pubr() +
    scale_color_manual(values = c("#0E7FFE", "#FE02FF")) + 
    theme(legend.position = 'none')


# 
# 
# 
# up_bois = hour_one$results_table %>% filter(padj < 0.1) %>% slice_max(log2FoldChange,n = 12) %>% pull(gene_name)
# 
# new_ratio_bayesian_p_de %>% 
#     filter(time == 1) %>% 
#     filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) %>% 
#     filter(gene_name %in% up_bois) %>% 
#     select(mean_bdnf_ntr,gene_name,mean_control_ntr,log2FoldChange)  %>% 
#      melt(id.vars = c("gene_name","log2FoldChange")) %>% 
#     mutate(gene_name =fct_reorder(gene_name,log2FoldChange)) %>% 
#     ggplot(aes(x = gene_name, y = value,
#                fill = variable)) + 
#     geom_col( position = "dodge") + 
#     coord_flip() +
#     theme_minimal() 
# 
# 
# estimate_list_full[time == 1 & symbol %in% 'ARC'] %>% 
#     ggplot(aes(x = symbol,y = map,fill = condition)) + 
#     geom_errorbar(aes(ymin = `005quantile`, ymax = `095quantile`)) + 
#     geom_point(pch = 21) + 
#     coord_flip()
# 
# library(clusterProfiler)
# 
# search_kegg_organism('hsa', by='kegg_code')
# 
# hr1 = new_ratio_bayesian_p_de %>% 
#     filter(time == 1) %>% 
#     filter(total_rna_sig != "not_significant") %>% 
#     select(gene,time,new_rna_sig,total_rna_sig,log2FoldChange) %>% 
#     unique() %>% 
#     mutate(gene = gsub("\\..*","",gene))
# 
# entrz_lookup = gprofiler2::gconvert(hr1$gene,target = "ENTREZGENE_ACC")
# 
# gene1hr = hr1 %>% 
#     left_join(entrz_lookup,by = c("gene" = "input")) %>% 
#     filter(total_rna_sig != "not_significant") %>% 
#     pull(log2FoldChange)
# 
# names(gene1hr) = hr1 %>% 
#     left_join(entrz_lookup,by = c("gene" = "input")) %>% 
#     filter(total_rna_sig != "not_significant") %>% 
#     pull(target)
# 
# enrichme = hr1 %>% 
#     left_join(entrz_lookup,by = c("gene" = "input")) %>% 
#     filter(total_rna_sig != "not_significant") %>% 
#     pull(target)
# 
# kk <- enrichKEGG(gene         = enrichme,
#                  organism     = 'hsa',
#                  pvalueCutoff = 0.05)
# 
# kk@result %>% slice_min(p.adjust,n = 10)
# library("pathview")
# hsa04110 <- pathview(gene.data  = gene1hr,
#                      pathway.id = "hsa04370",
#                      species    = "hsa")
# 
# 
# lmn = readxl::read_excel('/Users/annaleigh/Downloads/first_pass_half_life_stability_measurements.xlsx') |> 
#     as.data.table() 
# lmn = lmn |> 
#     mutate(gene = gsub("\\..*", "", Geneid)) 
# 
# cortical = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/cort_half_lifes_decay_and_synthesis_rates.csv")
# cortical = cortical |> 
#     mutate(gene = gsub("\\..*", "", gene))
# 
# 
# new_ratio_bayesian_p_de %>% 
#     mutate(gene = gsub("\\..*", "", gene)) %>%
#     filter(new_rna_sig != 'unclear') |> 
#     dplyr::select(gene,gene_name,time,new_rna_sig,total_rna_sig,log2FoldChange) |> 
#     # left_join(cortical) |> 
#     left_join(lmn[,.(gene,half_life_control)]) |> 
#     unique() |> 
#     ggplot(aes(x = new_rna_sig,
#                y = half_life_control,fill = total_rna_sig)) +
#     geom_boxplot() +
#     facet_wrap(~time,nrow = 3) + ggpubr::stat_compare_means(label = 'p.signif') +
#     scale_y_continuous(trans = scales::pseudo_log_trans()) +
#     ggtitle('Half-life in control i3 LMN')
# 
# 
# new_ratio_bayesian_p_de %>% 
#     mutate(gene = gsub("\\..*", "", gene)) %>%
#     filter(new_rna_sig != 'unclear') |> 
#     dplyr::select(gene,gene_name,time,new_rna_sig,total_rna_sig,log2FoldChange) |> 
#     left_join(cortical) |>
#     # left_join(lmn[,.(gene,half_life_control)]) |> 
#     unique() |> 
#     ggplot(aes(x = new_rna_sig,
#                y = half_life_ctrl,fill = total_rna_sig)) +
#     geom_boxplot() +
#     facet_wrap(~time,nrow = 3) + ggpubr::stat_compare_means(label = 'p.signif') +
#     ggtitle('Half-life in control i3 cortical')
