
bdnf1hr_input = fread("data/phospho_bdnf_1hr_for_enrichment_psp_version2.tsv")
bdnf6hr_input = fread("data/phospho_bdnf_6hr_for_enrichment_psp_version2.tsv")
one_hour_prediction = fread("/Users/annaleigh/Documents/GitHub/lmn_bdnf/data/phospho_bdnf_1hr_kinome_library_prediction.tsv")
 

these_kianses = c("P70S6K", "P70S6KB", "P90RSK", "RSK2", "RSK3", "RSK4","ERK1","ERK2")

figure_s3a_proteins = c("MACF1","MAP1A","MICAL3","MAP1B","MLLT4",
                        "MAP2","MPRIP","PALLD",
                        "MAP7D1","DMTN","MAPRE1","STMN1","CAMSAP3",
                        "IQGAP2","MYO5A","CLASP1","CLASP2")


tmp = one_hour_prediction %>% 
    filter(gene_names %in% figure_s3a_proteins) %>% 
    arrange(`Percentile Rank`) %>% 
    filter(kinase %in% these_kianses) %>% 
    left_join(bdnf1hr_input %>% distinct(gene_names,remapped_amino,phospho_position,adj_p_val)) %>% 
    filter(adj_p_val < 0.1) %>% 
    mutate(plot_site = glue::glue('{gene_names}-{remapped_amino}{phospho_position}')) %>% 
    distinct(plot_site,kinase,`Percentile Rank`) %>% 
    arrange(kinase) %>% 
    pivot_wider(names_from = 'kinase',values_from = 'Percentile Rank') %>% 
    as.data.frame()

rownames(tmp) = tmp$plot_site
tmp$plot_site = NULL

annotation_df = bdnf1hr_input %>% distinct(gene_names,remapped_amino,phospho_position,p_value,adj_p_val,log_fc) %>% 
    filter(adj_p_val < 0.1) %>% 
    mutate(plot_site = glue::glue('{gene_names}-{remapped_amino}{phospho_position}')) %>% 
    left_join(tmp %>% rownames_to_column('plot_site')) %>% 
    filter(!is.na(ERK1)) %>% 
    distinct(plot_site,log_fc)
    
rownames(annotation_df) = annotation_df$plot_site    
annotation_df$plot_site = NULL
annotation_colors = list(
    log_fc = colorRampPalette(c("blue", "white", "red"))(100)
)
default_colors = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

pheatmap::pheatmap(tmp,cutree_rows = 3,cluster_cols = FALSE,
                   annotation_row = annotation_df,
                   annotation_colors = annotation_colors,
                   color = default_colors)    
    
new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv") 
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de %>% 
    dplyr::select(-gene_name) %>% 
    mutate(gene = gsub("\\..*", "", gene)) %>%  
    left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol), by = c("gene" = "ensgene")) %>% 
    unique()

one_hour_prediction %>% 
    left_join(bdnf1hr_input %>% distinct(gene_names,remapped_amino,phospho_position,log_fc,adj_p_val)) %>% 
    filter(`Percentile Rank` <=5) %>% 
    filter(kinase %in% these_kianses) %>% 
    filter(adj_p_val < 0.1) %>% 
    filter()
    distinct(gene_names) %>% clipr::write_clip()
    # filter(kinase %in% c("AKT1","AKT2","AKT3"))

up_reg_phos_1hr = clusterProfiler::enrichGO(sig_up_one,
                                            universe = unique(new_ratio_bayesian_p_de$gene),
                                            keyType = 'ENSEMBL',
                                            OrgDb = org.Hs.eg.db,
                                            ont = 'ALL',
                                            readable = TRUE)