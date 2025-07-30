library(data.table)
library(tidyverse)
source("scripts/map_phosphorylation_sites.R")

bdnf1hr_v_control = fread("data/bdnf1hr_v_control_phos.csv")

bdnf6hr_v_control = fread('data/bdnf6hr_v_control_phos.csv')


bdnf1hr_v_control2 = map_phosphorylation_sites(bdnf1hr_v_control)
# OLD uniprot id
bdnf1hr_v_control2 = bdnf1hr_v_control2 %>% filter(!is.na(phospho_site))

bdnf6hr_v_control2 = map_phosphorylation_sites(bdnf6hr_v_control)
bdnf6hr_v_control2 = bdnf6hr_v_control2 %>% filter(!is.na(phospho_site))

bdnf1hr_v_control2 |> 
    mutate(remapped_amino = str_sub(clean_peptide,phospho_site,phospho_site)) %>% 
    distinct(uniprot_id,gene_names,remapped_amino,phospho_position,log_fc,p_value,adj_p_val,phospho_sty_probabilities,seq) %>% 
    mutate(protein_length = nchar(seq)) %>% 
    mutate(needs_padding = case_when(phospho_position - 7 <= 0 ~ "front",
                                     phospho_position + 7 > protein_length ~ 'back',
                                     TRUE ~ 'no')) %>% 
    mutate(padding_length = case_when(needs_padding == "front" ~ 7 - (phospho_position - 1),
                                      needs_padding == "back" ~ 7 - (protein_length - phospho_position),
                                      TRUE ~ 0)) %>% 
    mutate(padding = ifelse(
        padding_length == 0,
        '',
        mapply(function(n) paste0(rep("_", n), collapse = ""), padding_length)
    )) %>% 
    mutate(uniprot_flank_larger = case_when(needs_padding == 'no' ~ stringr::str_sub(seq,as.numeric(phospho_position) - 7, as.numeric(phospho_position) + 7),
                                            needs_padding == 'front' ~ glue::glue("{padding}{stringr::str_sub(seq,1, as.numeric(phospho_position) + 7)}"),
                                            needs_padding == 'back' ~ glue::glue("{stringr::str_sub(seq,as.numeric(phospho_position) - 7,protein_length)}{padding}"))) %>% 
    select(gene_names,phospho_position,uniprot_flank_larger,remapped_amino,needs_padding,padding,phospho_sty_probabilities,log_fc,p_value,adj_p_val) %>% 
    fwrite('data/phospho_bdnf_1hr_for_enrichment_psp_version2.tsv',sep = '\t')


bdnf6hr_v_control2 |> 
    mutate(remapped_amino = str_sub(clean_peptide,phospho_site,phospho_site)) %>% 
    distinct(uniprot_id,gene_names,remapped_amino,phospho_position,log_fc,p_value,adj_p_val,phospho_sty_probabilities,seq) %>% 
    mutate(protein_length = nchar(seq)) %>% 
    mutate(needs_padding = case_when(phospho_position - 7 <= 0 ~ "front",
                                     phospho_position + 7 > protein_length ~ 'back',
                                     TRUE ~ 'no')) %>% 
    mutate(padding_length = case_when(needs_padding == "front" ~ 7 - (phospho_position - 1),
                                      needs_padding == "back" ~ 7 - (protein_length - phospho_position),
                                      TRUE ~ 0)) %>% 
    mutate(padding = ifelse(
        padding_length == 0,
        '',
        mapply(function(n) paste0(rep("_", n), collapse = ""), padding_length)
    )) %>% 
    mutate(uniprot_flank_larger = case_when(needs_padding == 'no' ~ stringr::str_sub(seq,as.numeric(phospho_position) - 7, as.numeric(phospho_position) + 7),
                                            needs_padding == 'front' ~ glue::glue("{padding}{stringr::str_sub(seq,1, as.numeric(phospho_position) + 7)}"),
                                            needs_padding == 'back' ~ glue::glue("{stringr::str_sub(seq,as.numeric(phospho_position) - 7,protein_length)}{padding}"))) %>% 
    fwrite('data/phospho_bdnf_6hr_for_enrichment_psp_version2.tsv',sep = '\t')



# bdnf1hr_kinase_enrich = fread('data/enrichment-analysis-result-table_onehour_pvalue.txt')
# bdnf6hr_kinase_enrich = fread('data/enrichment-analysis-result-table_sixhour_pvalue.txt')

bdnf1hr_kinase_enrich = fread('/Users/annaleigh/Documents/GitHub/lmn_bdnf/data/kinase_library_prediction_1hr.csv')
bdnf6hr_kinase_enrich = fread('/Users/annaleigh/Documents/GitHub/lmn_bdnf/data/kinase_library_prediction_6hr.csv')

# one_kin = GOfuncR::get_anno_categories(unique(bdnf1hr_kinase_enrich$kinase))
# these_mapk = one_kin |> filter(grepl("MAPK",name)) |> pull(gene) |> unique()
# top_p = bdnf1hr_kinase_enrich |> 
#     slice_max(n = 5,dominant_enrichment_value_log2) |> 
#     pull(kinase) |> unique()
#re-make that plot

one_plot = bdnf1hr_kinase_enrich |> 
    mutate(plotAlpha = most_sig_fisher_pval < 0.01) |> 
    mutate(plotKinase = ifelse(most_sig_fisher_pval < 0.01,V1,NA_character_)) |> 
    # mutate(plotKinase = ifelse(dominant_p_value < 0.01 & (kinase %in% these_mapk| kinase %in% top_p),kinase,NA_character_)) |> 
    # filter(!(dominant_direction == 'upregulated set' & dominant_p_value < 0.1 & dominant_enrichment_value_log2 < 0)) |> 
    # mutate(plotFill = case_when(dominant_p_value >0.01 ~ 'grey',
    #                             dominant_direction == 'upregulated set' & dominant_p_value < 0.1 ~ 'red',
    #                             dominant_direction == 'downregulated set'~ 'blue')) |> 
    mutate(plotKinase = case_when(plotKinase == 'P90RSK' ~ "RSK1",
                                  plotKinase == 'P70S6K' ~ "S6K1",
                                  plotKinase == 'P70S6KB' ~ "S6K2",
                                  TRUE ~ plotKinase)) %>% 
    ggplot(aes(x = most_sig_log2_freq_factor,
               y = -log10(most_sig_fisher_pval),
               alpha = plotAlpha)) + 
    geom_point(pch = 21,fill = 'black') +
    theme_minimal() + 
    geom_hline(yintercept = 2,linetype = 2) + 
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = c("#446dea",'grey','red')) +
    scale_y_continuous(trans = scales::pseudo_log_trans()) + 
    ggrepel::geom_label_repel(aes(label = plotKinase)) + 
    xlab(expression(paste("Lo", g[2], " kinase enrichment"))) +
    ylab(expression(paste("-Lo", g[10], " enrichment p-value"))) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 18)) +
    ggtitle('1 h BDNF treatment')

six_plot = bdnf6hr_kinase_enrich |> 
    mutate(plotAlpha = most_sig_fisher_pval < 0.01) |> 
    mutate(plotKinase = ifelse(most_sig_fisher_pval < 0.01,V1,NA_character_)) |> 
    mutate(plotKinase = case_when(plotKinase == 'P90RSK' ~ "RSK1",
                                  plotKinase == 'P70S6K' ~ "S6K1",
                                  plotKinase == 'P70S6KB' ~ "S6K2",
                                  TRUE ~ plotKinase)) %>% 
    ggplot(aes(x = most_sig_log2_freq_factor,
               y = -log10(most_sig_fisher_pval),
               alpha = plotAlpha)) + 
    geom_point(pch = 21,fill = 'black') +
    theme_minimal() + 
    geom_hline(yintercept = 2,linetype = 2) + 
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = c("#446dea",'grey','red')) +
    scale_y_continuous(trans = scales::pseudo_log_trans()) + 
    ggrepel::geom_label_repel(aes(label = plotKinase)) + 
    xlab(expression(paste("Lo", g[2], " kinase enrichment"))) +
    ylab(expression(paste("-Lo", g[10], " enrichment p-value"))) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 18)) +
    ggtitle('6 h BDNF treatment') 

ggpubr::ggarrange(one_plot,six_plot,nrow = 1)
