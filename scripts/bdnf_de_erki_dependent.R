library(clusterProfiler)
library(tidyverse)
library(ggpubr)
library(ggrepel)
new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv") 

new_ratio_bayesian_p_de = new_ratio_bayesian_p_de %>% 
    dplyr::select(-gene_name) %>% 
    mutate(gene = gsub("\\..*", "", gene)) %>%  
    left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol), by = c("gene" = "ensgene")) %>% 
    unique()

bg_expressed_in_LMN = new_ratio_bayesian_p_de[baseMean > 15,gene]
# read in data, load helper function --------------------------------------


add_the_annotation <- function(dt, column) {
    
    
    # Perform gconvert and select input/target columns
    converted_data <- gprofiler2::gconvert(dt |> pull(column),numeric_ns = 'ENTREZGENE_ACC')
    
    converted_data = converted_data |> 
        select(input, target)
    
    # Perform left join on uniprot or x1 column (depending on input)
    xname = 'input'
    dt <- dt %>%
        left_join(converted_data, by=setNames(xname, column))  # Removed extra '='
    
    # Return the modified data table
    return(dt)
}

bdnf_v_control_axon_total = fread("/Users/annaleigh/Desktop/bdnf_v_control_axon_total.csv")
bdnf_v_control_soma_total = fread("/Users/annaleigh/Desktop/bdnf_v_control_soma_total.csv")

erki_v_bdnf_soma = fread("/Users/annaleigh/Downloads/DE results/erki_v_bdnf_soma_total.csv")
erki_v_bdnf_axon = fread("/Users/annaleigh/Downloads/DE results/erki_v_bdnf_axon_total.csv")

erki_v_control_soma = fread("/Users/annaleigh/Downloads/DE results/erki_v_control_soma_total.csv")
erki_v_control_axon = fread("/Users/annaleigh/Downloads/DE results/erki_v_control_axon_total.csv")



bdnf_effect_both = bdnf_v_control_soma_total %>% 
    full_join(bdnf_v_control_axon_total,by = c("un"),
              suffix = c("_soma","_axon"))
    
bdnf_effect_both = bdnf_effect_both %>% 
    mutate(sig = case_when(adj_p_val_soma < 0.1 & abs(log_fc_soma) > 0.5 & 
                               adj_p_val_axon < 0.1 & abs(log_fc_axon) > 0.5 ~ "Sig. Both",
                           adj_p_val_soma < 0.1 & abs(log_fc_soma) > 0.5 ~ ' Sig. Soma',
                           adj_p_val_axon < 0.1 & abs(log_fc_axon) > 0.5 ~ " Sig. Axon")) %>% 
    mutate(axon = case_when(log_fc_axon > 0.5 & adj_p_val_axon < 0.1 ~ "Upregulated",
                            log_fc_axon < -0.5 & adj_p_val_axon < 0.1 ~ "Downregulated",
                                           TRUE ~ "NS")) %>% 
    mutate(soma = case_when(log_fc_soma > 0.5 & adj_p_val_soma < 0.1 ~ "Upregulated",
                            log_fc_soma < -0.5 & adj_p_val_soma < 0.1 ~ "Downregulated",
                            TRUE ~ "NS")) 

# Step 1: Calculate the number of unique values in each 'sig' category
counts <- bdnf_effect_both %>% 
    group_by(sig) %>% 
    summarize(count = n_distinct(un))

# Create a text label containing this information
label_text <- paste(counts$sig, counts$count, sep = ": ", collapse = "\n")


ggplot() + 
    geom_point(data = bdnf_effect_both,
               aes(x = log_fc_soma,
                   y = log_fc_axon)) + 
    geom_point(data = subset(bdnf_effect_both,!is.na(sig)),
               aes(x = log_fc_soma,
                   y = log_fc_axon,color = sig)) + 
    # Add stat_cor for data where sig is NA
    stat_cor(data = subset(bdnf_effect_both, is.na(sig)),
             aes(x = log_fc_soma, y = log_fc_axon),
             label.y.npc = 1) +  # Adjust label.y to position this label slightly higher
    # Add stat_cor for data where sig is "Sig. Both"
    stat_cor(data = subset(bdnf_effect_both, sig == 'Sig. Both'),
             aes(x = log_fc_soma, y = log_fc_axon),
             label.y.npc = 0.95,color = "#ab8aceff")  +  # Adjust label.y to p
    theme_classic() +
    geom_hline(yintercept = 0,lintype = 'dotted') +
    geom_vline(xintercept = 0,lintype = 'dotted') + 
    ylab("Log2FC BDNF vs Control - Axon") +
    xlab("Log2FC BDNF vs Control - Soma") +
    annotate("text", 
             x = Inf, y = -Inf, 
             label = label_text, 
             hjust = 1.1, vjust = -1.1, 
             size = 5, 
             color = "black") +
    scale_color_manual(values = c("#abd7e6ff","#e08e8eff","#ab8aceff"))
    

# GO BDNF effect ---------------------------------------------------
library(org.Hs.eg.db)
amigoingmad()

both_go_dt = bdnf_effect_both %>% 
    filter(axon != "NS" | soma != 'NS' ) %>% 
    mutate(direction = as.character(glue::glue("Soma: {soma}\n Axon: {axon}"))) %>% 
    mutate(conver = case_when(axon == soma ~ "Convergent",
                            axon == 'NS' ~ "Soma only",
                            soma == 'NS' ~ 'Axon only',
                            TRUE ~ 'Divergent')) %>% 
    mutate(go_fc = case_when(conver == 'Axon only' ~ axon,
                             conver == 'Soma only' ~ soma,
                             conver == 'Convergent' ~ soma,
                             conver == 'Divergent' ~ direction)) %>% 
    separate(un,into = c("uniprot"),remove = FALSE) %>%
    add_the_annotation('uniprot') %>%
    filter(conver != 'Divergent') %>% 
    distinct(target,conver,go_fc)

both_go3 = clusterProfiler::compareCluster(target~conver+go_fc,
                                           data = both_go_dt,
                                           universe = bg_expressed_in_LMN,
                                           keyType = 'ENSEMBL',
                                           OrgDb = org.Hs.eg.db,
                                           ont = 'ALL',
                                           readable = TRUE)

both_go3@compareClusterResult = both_go3@compareClusterResult %>% 
    mutate(conver = fct_relevel(conver,"Soma only")) %>% 
    separate(GeneRatio,remove = FALSE,into = c("top","bottom"),convert = TRUE) %>% 
    filter(top >1)

library(clusterProfiler)
dotplot(both_go3, x="go_fc") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10)) +
    xlab(element_blank()) +
    facet_grid(rows = vars(conver), cols = vars(ONTOLOGY), scales = 'free') + 
    coord_flip()

cnetplot(both_go3) + scale_fill_manual(values = c("#AAAFD7",
                                                  "#DE71A8",
                                                  "#B2E26D",
                                                  "#B15CDC",
                                                  "#DBA879"))


# soma Volcano plot BDNF vs Control --------------------------------------------

bdnf_v_control_soma_total %>% 
    mutate(short_name = str_replace(un, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = ifelse(adj_p_val < 0.1 & abs(log_fc) > 0.5,short_name,NA_character_)) %>% 
    ggplot( aes(x = log_fc, y = -log10(adj_p_val))) +
    geom_point(aes(alpha = adj_p_val < 0.1 & abs(log_fc) > 0.5), size = 2) +
    scale_color_manual(values = c("red", "blue"), name = "Significance") +
    geom_hline(yintercept = -log10(0.1), linetype = "dotted") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    labs( x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() +
    xlim(c(-5.5, 5)) +
    ylim(c(0, 6.2)) +
    theme(legend.position = 'none') 
# soma Volcano plot ERKI vs Control --------------------------------------------

erki_v_control_soma %>% 
    mutate(short_name = str_replace(un, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = ifelse(adj_p_val < 0.1 & abs(log_fc) > 0.5,short_name,NA_character_)) %>% 
    ggplot( aes(x = log_fc, y = -log10(adj_p_val))) +
    geom_point(aes(alpha = adj_p_val < 0.1 & abs(log_fc) > 0.5), size = 2) +
    scale_color_manual(values = c("red", "blue"), name = "Significance") +
    geom_hline(yintercept = -log10(0.1), linetype = "dotted") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    labs( x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() +
    xlim(c(-5.5, 5)) +
    ylim(c(0, 6.2)) +
    theme(legend.position = 'none') 

# soma Volcano plot ERKI vs BDNF --------------------------------------------

erki_v_bdnf_soma %>% 
    mutate(short_name = str_replace(un, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = ifelse(adj_p_val < 0.1 & abs(log_fc) > 0.5,short_name,NA_character_)) %>% 
    ggplot( aes(x = log_fc, y = -log10(adj_p_val))) +
    geom_point(aes(alpha = adj_p_val < 0.1 & abs(log_fc) > 0.5), size = 2) +
    scale_color_manual(values = c("red", "blue"), name = "Significance") +
    geom_hline(yintercept = -log10(0.1), linetype = "dotted") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    labs( x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() +
    xlim(c(-5.5, 5)) +
    ylim(c(0, 6.2)) +
    theme(legend.position = 'none') 

# axon Volcano plot BDNF vs Control --------------------------------------------

bdnf_v_control_axon_total %>% 
    mutate(short_name = str_replace(un, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = ifelse(adj_p_val < 0.1 & abs(log_fc) > 0.5,short_name,NA_character_)) %>% 
    ggplot( aes(x = log_fc, y = -log10(adj_p_val))) +
    geom_point(aes(alpha = adj_p_val < 0.1 & abs(log_fc) > 0.5), size = 2) +
    scale_color_manual(values = c("red", "blue"), name = "Significance") +
    geom_hline(yintercept = -log10(0.1), linetype = "dotted") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    labs( x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() +
    xlim(c(-5.5, 5)) +
    ylim(c(0, 6.2)) +
    theme(legend.position = 'none')

# axon Volcano plot ERKI vs Control --------------------------------------------

erki_v_control_axon %>% 
    mutate(short_name = str_replace(un, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = ifelse(adj_p_val < 0.1 & abs(log_fc) > 0.5,short_name,NA_character_)) %>% 
    ggplot( aes(x = log_fc, y = -log10(adj_p_val))) +
    geom_point(aes(alpha = adj_p_val < 0.1 & abs(log_fc) > 0.5), size = 2) +
    scale_color_manual(values = c("red", "blue"), name = "Significance") +
    geom_hline(yintercept = -log10(0.1), linetype = "dotted") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    labs( x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() +
    xlim(c(-5.5, 5)) +
    ylim(c(0, 6.2)) +
    theme(legend.position = 'none')
# axon Volcano plot ERKI vs BDNF --------------------------------------------

erki_v_bdnf_axon %>% 
    mutate(short_name = str_replace(un, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = ifelse(adj_p_val < 0.1 & abs(log_fc) > 0.5,short_name,NA_character_)) %>% 
    ggplot( aes(x = log_fc, y = -log10(adj_p_val))) +
    geom_point(aes(alpha = adj_p_val < 0.1 & abs(log_fc) > 0.5), size = 2) +
    scale_color_manual(values = c("red", "blue"), name = "Significance") +
    geom_hline(yintercept = -log10(0.1), linetype = "dotted") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
    labs( x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_classic() +
    xlim(c(-5.5, 5)) +
    ylim(c(0, 6.2)) +
    theme(legend.position = 'none') 




# define the proteins that are DE ERK vs BDNF -----------------------------

de_erki_axon = erki_v_bdnf_axon %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) >0.5)
# define the proteins that are BDNF v control -----------------------------
de_bdnf_de_erki = bdnf_v_control_axon_total %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) >0.5) %>% 
    filter(un %in% de_erki_axon$un) %>% 
    pull(un)



# axon fc compare with ERKI depend highligh ----------------------------
bdnf_v_control_soma_total %>% 
    filter(un %in% up_bdnf_soma | un %in% down_bdnf_soma) %>% 
    filter(un %in% de_erki_soma$un) %>% 
    full_join(erki_v_control_soma,by = 'un', suffix = c("_bdnf","_ierk")) %>% 
    ggplot(aes(x = log_fc_bdnf,
               y = log_fc_ierk)) + 
    geom_point(alpha = 0.3) + 
    ylim(c(-3,6)) + 
    xlim(c(-3,6)) + 
    geom_abline() + 
    # geom_text_repel(aes(label = un)) + 
    theme_classic() + 
    geom_vline(xintercept = -0.5,linetype = 'dotted') +
    geom_vline(xintercept = 0.5,linetype = 'dotted') +
    geom_hline(yintercept = 0.5,linetype = 'dotted') + 
    geom_hline(yintercept = -0.5,linetype = 'dotted')
    
    


# axon spaghetti -----------------------------------------------------------
axon_erki_dependent = erki_v_control_axon %>% 
    full_join(bdnf_v_control_axon_total,by = 'un', suffix = c("_erki","_bdnf")) %>% 
    filter(un %in% de_bdnf_de_erki) %>% 
    select(un,log_fc_erki,log_fc_bdnf,adj_p_val_erki,adj_p_val_bdnf) %>% 
    mutate(up_down = ifelse(log_fc_bdnf > 0, "up","down")) 

axon_erki_dependent %>% 
    melt(id.vars = c("un",'up_down')) %>% 
    mutate(variable = fct_relevel(variable,"log_fc_bdnf")) %>% 
    ggplot(aes(x = variable,
               y = value,
               group = un,color = up_down)) + 
    geom_line(alpha = 0.3,show.legend = FALSE) + 
    theme_classic() + 
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted") +
    geom_hline(yintercept = 0) + 
    ylab("Log2 FC") + 
    xlab(element_blank()) +
    scale_x_discrete(labels=c("BDNF\nvs\ncontrol", "iERK + BDNF\nvs\nControl")) + 
    scale_color_manual(values = c("#35666bff","#BC2C1A"))


# axon GO erki dependent -------------------------------------------------------

axon_go_dt_erki = axon_erki_dependent %>% 
    separate(un,into = c("uniprot"),remove = FALSE) %>%
    add_the_annotation('uniprot') %>%
    distinct(target,up_down)

go_erki_axon = clusterProfiler::compareCluster(target~up_down,
                                           data = axon_go_dt_erki,
                                           universe = bg_expressed_in_LMN,
                                           keyType = 'ENSEMBL',
                                           OrgDb = org.Hs.eg.db,
                                           ont = 'ALL',
                                           readable = TRUE)

cnetplot(go_erki_axon,shadowtext = 'none') + scale_fill_manual(values = c("#7FB7BE","#BC2C1A"))
dotplot(go_erki_axon) + 
    facet_wrap(~up_down,scales = 'free')
# soma define the proteins that are DE ERK vs BDNF -----------------------------

de_erki_soma = erki_v_bdnf_soma %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) >0.5)
# soma define the proteins that are BDNF v control -----------------------------
soma_de_bdnf_de_erki = bdnf_v_control_soma_total %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) >0.5) %>% 
    filter(un %in% de_erki_soma$un) %>% 
    pull(un)
# soma spaghetti -----------------------------------------------------------
erki_v_control_soma %>%
    full_join(bdnf_v_control_soma_total,by = 'un', suffix = c("_erki","_bdnf")) %>% 
    filter(un %in% soma_de_bdnf_de_erki) %>% 
    mutate(sig_erki = adj_p_val_erki < 0.1 & abs(log_fc_erki) > 0.5) %>% 
    mutate(up_down = ifelse(log_fc_bdnf > 0, "up","down")) %>% 
    mutate(erki_effect = case_when(sig_erki == FALSE & log_fc_bdnf > 0 ~ 'BDNF upregulated - ERKi negated',
                                   sig_erki == FALSE & log_fc_bdnf < 0 ~ 'BDNF downregulated - ERKi negated',
                                   sig_erki & log_fc_erki * log_fc_bdnf < 0 & log_fc_bdnf > 0  ~ 'BDNF upregulated - ERKi downregulated',
                                   sig_erki & log_fc_erki * log_fc_bdnf < 0 & log_fc_bdnf < 0  ~ 'BDNF downregulated - ERKi upregulated',
                                   sig_erki & log_fc_erki < log_fc_bdnf & log_fc_bdnf > 0 ~ 'BDNF upregulated - mediated ERKi',
                                   sig_erki & log_fc_erki > log_fc_bdnf & log_fc_bdnf < 0 ~ 'BDNF downregulated - mediated ERKi',
                                   sig_erki & log_fc_erki > log_fc_bdnf & log_fc_bdnf > 0 ~ 'BDNF upregulated - potentiated ERKi',
                                   sig_erki & log_fc_erki < log_fc_bdnf & log_fc_bdnf < 0 ~ 'BDNF downregulated - potentiated ERKi',
    )) %>% 
    select(un,log_fc_erki,log_fc_bdnf,up_down,erki_effect) %>% 
    melt(id.vars = c("un",'up_down','erki_effect')) %>% 
    mutate(variable = fct_relevel(variable,"log_fc_bdnf")) %>% 
    mutate(negated = (grepl("negated",erki_effect))) %>% 
    mutate(up_down = paste0(up_down,negated)) %>% 
    ggplot(aes(x = variable,
               y = value,
               group = un,
               color = up_down)) + 
    geom_line(alpha = 0.3,show.legend = FALSE) + 
    theme_classic() + 
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted") +
    geom_hline(yintercept = 0) + 
    ylab("Log2 FC") + 
    xlab(element_blank()) +
    scale_x_discrete(labels=c("BDNF\nvs\ncontrol", "iERK + BDNF\nvs\nControl")) + 
    scale_color_manual(values = c(downFALSE = "#7FB7BE",
                                  upFALSE = "#BC2C1A",
                                  downTRUE = "#35666bff",
                                  upTRUE = "#470f0aff")) 


t1 = erki_v_control_soma %>% 
    full_join(bdnf_v_control_soma_total,by = 'un', suffix = c("_erki","_bdnf")) %>% 
    filter(un %in% soma_de_bdnf_de_erki) %>% 
    mutate(sig_erki = adj_p_val_erki < 0.1 & abs(log_fc_erki) > 0.5) %>% 
    select(un,log_fc_erki,log_fc_bdnf,sig_erki) %>% 
    mutate(up_down = ifelse(log_fc_bdnf > 0, "up","down")) %>% 
    filter(sig_erki == TRUE) %>% 
    mutate(erki_effect = case_when(log_fc_erki > log_fc_bdnf ~ "iERK enhanced",
        log_fc_erki * log_fc_bdnf < 0 ~ "Divergent",
                                   log_fc_erki * log_fc_bdnf > 0 ~ "Convergent")) %>% 
    filter(erki_effect == 'iERK enhanced')


t1 %>% 
    separate(un,into = c("uniprot"),remove = FALSE) %>%
    add_the_annotation('uniprot') %>%
    distinct(target,conver,go_fc)


new_soma_imputed_f %>% 
    rownames_to_column('un') %>% 
    filter(un %in% t1$un) %>% 
    column_to_rownames('un') %>% 
    pheatmap::pheatmap(scale = 'row',
                       cluster_cols = FALSE)



erki_v_control_soma %>% 
    ggplot(aes(x = log_fc, y = -log10(adj_p_val))) +
    geom_point(aes(alpha = adj_p_val < 0.1 & abs(log_fc) > 0.5), size = 2) +
    scale_alpha_manual(values = c(0.4, 1), name = "Significance") +
    geom_hline(yintercept = -log10(0.1), linetype = "dotted", color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    labs(title = "ERKi + BDNF vs Control - soma", x = "Log2 Fold Change", y = "-log10(p-value)") +
    theme_minimal()


erki_v_control_soma %>% 
    left_join(erki_v_control_swapped,by = c("un"),
              suffix = c("_soma","_axon")) %>% 
    ggplot(aes(x = log_fc_soma,
               y = log_fc_axon)) + 
    geom_point() + 
    theme_classic() + 
    geom_hline(yintercept = 0,lintype = 'dotted') +
    geom_vline(xintercept = 0,lintype = 'dotted') + 
    ylab("Log2FC ERKI vs Control - Axon - Control & BDNF swapped") +
    xlab("Log2FC ERKI vs Control - Soma") + 
    ggtitle('All proteins')

erki_v_control_soma %>% 
    left_join(erki_v_control_swapped,by = c("un"),
              suffix = c("_soma","_axon")) %>% 
    filter(p_value_soma < 0.05 & p_value_axon < 0.05) %>% 
    
    ggplot(aes(x = log_fc_soma,
               y = log_fc_axon)) + 
    geom_point() + 
    theme_classic() + 
    geom_hline(yintercept = 0,lintype = 'dotted') +
    geom_vline(xintercept = 0,lintype = 'dotted') + 
    ylab("Log2FC ERKI vs Control - Axon - Control & BDNF swapped") +
    xlab("Log2FC ERKI vs Control - Soma") + 
    ggtitle('All proteins')



# axon - barely change erk -----------------------------------------------

bdnf_axon_bare_change_erki = axon_erki_dependent %>% 
    filter(abs(log_fc_erki) < 0.5) %>% 
    arrange(-log_fc_bdnf) %>% 
    separate(un,into = c("uniprot"),remove = FALSE) %>%
    add_the_annotation('uniprot') %>% 
    distinct(target,up_down)

go_bdnf_axon_bare_change_erki = clusterProfiler::compareCluster(target~up_down,
                                           data = bdnf_axon_bare_change_erki,
                                           universe = bg_expressed_in_LMN,
                                           keyType = 'ENSEMBL',
                                           OrgDb = org.Hs.eg.db,
                                           ont = 'ALL',
                                           readable = TRUE)
dotplot(go_bdnf_axon_bare_change_erki) +
    facet_wrap(~up_down) 
cnetplot(go_bdnf_axon_bare_change_erki) + scale_fill_manual(values = c("#38676dff","#49110b4c"))
# soma - barely change erk -----------------------------------------------
soma_erki_effects = erki_v_control_soma %>% 
    full_join(bdnf_v_control_soma_total,by = 'un', suffix = c("_erki","_bdnf")) %>% 
    filter(un %in% soma_de_bdnf_de_erki) %>% 
    mutate(sig_erki = adj_p_val_erki < 0.1 & abs(log_fc_erki) > 0.5) %>% 
    mutate(erki_effect = case_when(sig_erki == FALSE & log_fc_bdnf > 0 ~ 'BDNF upregulated - ERKi negated',
                                   sig_erki == FALSE & log_fc_bdnf < 0 ~ 'BDNF downregulated - ERKi negated',
                                   sig_erki & log_fc_erki * log_fc_bdnf < 0 & log_fc_bdnf > 0  ~ 'BDNF upregulated - ERKi downregulated',
                                   sig_erki & log_fc_erki * log_fc_bdnf < 0 & log_fc_bdnf < 0  ~ 'BDNF downregulated - ERKi upregulated',
                                   sig_erki & log_fc_erki < log_fc_bdnf & log_fc_bdnf > 0 ~ 'BDNF upregulated - mediated ERKi',
                                   sig_erki & log_fc_erki > log_fc_bdnf & log_fc_bdnf < 0 ~ 'BDNF downregulated - mediated ERKi',
                                   sig_erki & log_fc_erki > log_fc_bdnf & log_fc_bdnf > 0 ~ 'BDNF upregulated - potentiated ERKi',
                                   sig_erki & log_fc_erki < log_fc_bdnf & log_fc_bdnf < 0 ~ 'BDNF downregulated - potentiated ERKi',
                                   )) %>% 
    select(un,log_fc_erki,log_fc_bdnf,sig_erki,adj_p_val_erki,adj_p_val_bdnf,erki_effect) 

soma_erki_effects_go_dt = soma_erki_effects %>% 
    separate(un,into = c("uniprot"),remove = FALSE) %>%
    add_the_annotation('uniprot') %>% 
    distinct(target,erki_effect)

go_bdnf_soma_change_erki = clusterProfiler::compareCluster(target~erki_effect,
                                                                data = soma_erki_effects_go_dt,
                                                                universe = bg_expressed_in_LMN,
                                                                keyType = 'ENSEMBL',
                                                                OrgDb = org.Hs.eg.db,
                                                                ont = 'ALL',
                                                                readable = TRUE)

go_bdnf_soma_change_erki@compareClusterResult = go_bdnf_soma_change_erki@compareClusterResult %>% 
    filter(Count >=4)
dotplot(go_bdnf_soma_change_erki) +
    facet_wrap(~erki_effect,scales = 'free')
cnetplot(go_bdnf_soma_change_erki,shadowtext = 'none') + scale_fill_manual(values = c("#7FB7BE","#BC2C1A","#49110b4c"))


# clusterPlotBoth Axon + soma ---------------------------------------------
soma_erki_effects_go_dt = soma_erki_effects %>% 
    separate(un,into = c("uniprot"),remove = FALSE) %>%
    add_the_annotation('uniprot') %>% 
    distinct(target,erki_effect)

full_erki_go_dt = axon_go_dt_erki %>% 
    mutate(erki_effect = ifelse(up_down == "up",
                                "BDNF upregulated - ERKi negated",
                                "BDNF downregulated - ERKi negated")) %>% 
    mutate(compartment = 'axon') %>% 
    select(-up_down) %>% 
    rbind(soma_erki_effects_go_dt %>% 
              mutate(compartment = 'soma')) 

go_bdnf_full_erki = clusterProfiler::compareCluster(target~erki_effect+compartment,
                                                           data = full_erki_go_dt,
                                                           universe = bg_expressed_in_LMN,
                                                           keyType = 'ENSEMBL',
                                                           OrgDb = org.Hs.eg.db,
                                                           ont = 'ALL',
                                                           readable = TRUE)

go_bdnf_full_erki@compareClusterResult = go_bdnf_full_erki@compareClusterResult %>% 
    filter(Count >=4)
dotplot(go_bdnf_full_erki,x = 'compartment') +
    facet_wrap(~erki_effect,scales = 'free')
cnetplot(go_bdnf_full_erki,shadowtext = 'none')

k = go_bdnf_full_erki 

k@compareClusterResult = k@compareClusterResult %>% 
    filter(erki_effect == 'BDNF upregulated - ERKi negated') 

dotplot(k)
# eulerr plot erki depedent -----------------------------------------------


axon_erki_dependent = erki_v_control_axon %>% 
    full_join(bdnf_v_control_axon_total,by = 'un', suffix = c("_erki","_bdnf")) %>% 
    filter(un %in% de_bdnf_de_erki) %>% 
    select(un,log_fc_erki,log_fc_bdnf,adj_p_val_erki,adj_p_val_bdnf) %>% 
    mutate(up_down = ifelse(log_fc_bdnf > 0, "up","down")) 


up_bdnf_soma = bdnf_effect_both %>% filter(soma == 'Upregulated') %>% pull(un)
up_bdnf_axon = bdnf_effect_both %>% filter(axon == 'Upregulated') %>% pull(un)
down_bdnf_soma = bdnf_effect_both %>% filter(soma == 'Downregulated') %>% pull(un)
down_bdnf_axon = bdnf_effect_both %>% filter(axon == 'Downregulated') %>% pull(un)


up_erki_axon = erki_v_bdnf_axon %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) > 0.5) %>% 
    filter(log_fc > 0.5) %>% pull(un)

down_erki_axon = erki_v_bdnf_axon %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) > 0.5) %>% 
    filter(log_fc < -0.5) %>% pull(un)


up_erki_soma = erki_v_bdnf_soma %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) > 0.5) %>% 
    filter(log_fc > 0.5) %>% pull(un)

down_erki_soma = erki_v_bdnf_soma %>% 
    filter(adj_p_val < 0.1 & abs(log_fc) > 0.5) %>% 
    filter(log_fc < -0.5) %>% pull(un)

eulr_tbl = data.table(un = c(up_bdnf_axon,
                             down_bdnf_axon)) %>% 
    unique()

eulr_tbl %>% 
    mutate(up_bdnf_axon = un %in% up_bdnf_axon) %>% 
    mutate(down_bdnf_axon = un %in% down_bdnf_axon) %>% 
    mutate(down_erki_axon = un %in% down_erki_axon) %>% 
    mutate(up_erki_axon = un %in% up_erki_axon) %>% 
    filter(up_bdnf_axon == TRUE | down_bdnf_axon == TRUE) %>% 
    tibble::column_to_rownames("un") %>% 
    eulerr::euler() %>% 
    plot(quantities = TRUE)



eulr_tbl = data.table(un = c(up_bdnf_soma,
                             down_bdnf_soma)) %>% 
    unique()

eulr_tbl %>% 
    mutate(up_bdnf_soma = un %in% up_bdnf_soma) %>% 
    mutate(down_bdnf_soma = un %in% down_bdnf_soma) %>% 
    mutate(down_erki_soma = un %in% down_erki_soma) %>% 
    mutate(up_erki_soma = un %in% up_erki_soma) %>% 
    filter(up_bdnf_soma == TRUE | down_bdnf_soma == TRUE) %>% 
    tibble::column_to_rownames("un") %>% 
    eulerr::euler() %>% 
    plot(quantities = TRUE)

# axon fc compare with ERKI depend highligh ----------------------------
bdnf_v_control_axon_total %>% 
    filter(un %in% up_bdnf_axon | un %in% down_bdnf_axon) %>% 
    filter(un %in% de_erki_axon$un) %>% 
    full_join(erki_v_control_axon,by = 'un', suffix = c("_bdnf","_erki")) %>% 
    mutate(sig_erki = abs(log_fc_erki) > 0.5 & adj_p_val_erki < 0.1) %>% 
    mutate(short_name = str_replace(un, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = ifelse(adj_p_val_bdnf < 0.1 & abs(log_fc_bdnf) > 2,short_name,NA_character_)) %>% 
    mutate(erki_effect = case_when(sig_erki == FALSE & log_fc_bdnf > 0 ~ 'BDNF upregulated - ERKi negated',
                                   sig_erki == FALSE & log_fc_bdnf < 0 ~ 'BDNF downregulated - ERKi negated',
                                   sig_erki & log_fc_erki * log_fc_bdnf < 0 & log_fc_bdnf > 0  ~ 'BDNF upregulated - ERKi downregulated',
                                   sig_erki & log_fc_erki * log_fc_bdnf < 0 & log_fc_bdnf < 0  ~ 'BDNF downregulated - ERKi upregulated',
                                   sig_erki & log_fc_erki < log_fc_bdnf & log_fc_bdnf > 0 ~ 'BDNF upregulated - mediated ERKi',
                                   sig_erki & log_fc_erki > log_fc_bdnf & log_fc_bdnf < 0 ~ 'BDNF downregulated - mediated ERKi',
                                   sig_erki & log_fc_erki > log_fc_bdnf & log_fc_bdnf > 0 ~ 'BDNF upregulated - potentiated ERKi',
                                   sig_erki & log_fc_erki < log_fc_bdnf & log_fc_bdnf < 0 ~ 'BDNF downregulated - potentiated ERKi',
    )) %>% 
    filter(!is.na(log_fc_bdnf)) %>% 
    ggplot(aes(x = log_fc_bdnf,
               y = log_fc_erki,color = erki_effect)) + 
    geom_point() + 
    ylim(c(-5,5)) + 
    xlim(c(-5,5)) + 
    geom_abline() + 
    # geom_text_repel(aes(label = un)) + 
    theme_classic() + 
    geom_vline(xintercept = -0.5,linetype = 'dotted') +
    geom_vline(xintercept = 0.5,linetype = 'dotted') +
    geom_hline(yintercept = 0.5,linetype = 'dotted') + 
    geom_hline(yintercept = -0.5,linetype = 'dotted') + 
    scale_color_manual(values = c("#7FB7BE","#BC2C1A","#49110b4c")) +
    ylab("BDNF + ERKi vs control Log2FC") +
    xlab("BDNF vs control Log2FC") +
    geom_text_repel(aes(label = short_name))

soma_erki_effects %>% 
    ggplot(aes(x = log_fc_bdnf,
               y = log_fc_erki,
               color = erki_effect)) + 
    geom_point() + 
    geom_abline() + 
    # geom_text_repel(aes(label = un)) + 
    theme_classic() + 
    geom_vline(xintercept = -0.5,linetype = 'dotted') +
    geom_vline(xintercept = 0.5,linetype = 'dotted') +
    geom_hline(yintercept = 0.5,linetype = 'dotted') + 
    geom_hline(yintercept = -0.5,linetype = 'dotted') +
    ylim(c(-5,5)) +
    xlim(c(-5,5))

# Heatmap Top soma/axon erki dependent -------------------------------------------------------
imputed_total_soma = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/SOMA_impseq_Original_Result_1721384471.75033.csv")
imputed_total_axon = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/All_Normed_first_imseqrob_AxonOriginal_Result_1721383562.17019.csv")

axon_top_erki = axon_erki_dependent %>% arrange(-log_fc_bdnf)  %>% slice_max(log_fc_bdnf,n = 25) %>% pull(un)
axon_bottom_erki = axon_erki_dependent %>%  arrange(log_fc_bdnf) %>% slice_min(log_fc_bdnf,n = 25) %>% pull(un)

top_bottom_ax = c(axon_top_erki,axon_bottom_erki)

hot_ax = imputed_total_axon |> 
    # left_join(imputed_total_axon,by = c("V1")) |> 
    na.omit()  |>
    select(-matches("b1|d4|d1")) %>%
    select(-matches("c5|c4")) |>
    # mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
    filter(V1 %in% c(top_bottom_ax)) %>%
    unique() %>% 
    relocate(matches('control')) %>% 
    tibble::column_to_rownames('V1')

hot_ax = hot_ax[top_bottom_ax, ] %>% 
    tibble::rownames_to_column('V1') %>% 
    mutate(short_name = str_replace(V1, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = make.unique(short_name)) %>% 
    tibble::column_to_rownames("short_name") %>% 
    select(-V1) 


my_sample_col = tibble(samples = colnames(hot_ax)) |> 
    separate(samples,into = c("condition","compartment","replicate"),remove = FALSE) |> 
    tibble::column_to_rownames('samples') |> 
    select(-replicate)

my_colour = list(
    condition = c(control = "#F93FFF", BDNF = "#1B7DFD", iERKBDNF = "#72f54aff"),
    compartment = c(axon = "black", soma = "grey")
)


pheatmap::pheatmap(hot_ax,
                   annotation_colors = my_colour,
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE,
                   annotation_col = my_sample_col,
                   scale = 'row',
                   show_colnames = FALSE)


# heatmap soma top erki ---------------------------------------------------
soma_top_erki = soma_erki_effects %>% arrange(-log_fc_bdnf) %>% slice_max(log_fc_bdnf,n = 25) %>% pull(un)
soma_bottom_erki = soma_erki_effects %>% arrange(log_fc_bdnf) %>% slice_min(log_fc_bdnf,n = 25) %>% pull(un)
top_bottom_soma = c(soma_top_erki,soma_bottom_erki)

hot_soma = imputed_total_soma |> 
    # left_join(imputed_total_axon,by = c("V1")) |> 
    na.omit()  |>
    select(-matches("b1|d4|d1")) %>%
    select(-matches("c5|c4")) |>
    # mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
    filter(V1 %in% c(top_bottom_soma)) %>%
    unique() %>% 
    relocate(matches('control')) %>% 
    tibble::column_to_rownames('V1')

hot_soma = hot_soma[top_bottom_soma, ] %>% 
    tibble::rownames_to_column('V1') %>% 
    mutate(short_name = str_replace(V1, "^[^_]*_", "")) %>% 
    mutate(short_name = str_replace_all(short_name, "_HUMAN", "")) %>% 
    mutate(short_name = make.unique(short_name)) %>% 
    tibble::column_to_rownames("short_name") %>% 
    select(-V1) 


my_sample_col = tibble(samples = colnames(hot_soma)) |> 
    separate(samples,into = c("condition","compartment","replicate"),remove = FALSE) |> 
    tibble::column_to_rownames('samples') |> 
    select(-replicate)

my_colour = list(
    condition = c(control = "#F93FFF", BDNF = "#1B7DFD", iERKBDNF = "#72f54aff"),
    compartment = c(axon = "black", soma = "grey")
)


pheatmap::pheatmap(hot_soma,
                   annotation_colors = my_colour,
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE,
                   annotation_col = my_sample_col,
                   scale = 'row',
                   show_colnames = FALSE)





soma_erki_effects %>% 
    count(erki_effect)
axon_go_dt_erki %>% 
    mutate(erki_effect = ifelse(up_down == "up",
                                "BDNF upregulated - ERKi negated",
                                "BDNF downregulated - ERKi negated")) %>% 
    mutate(compartment = 'axon')