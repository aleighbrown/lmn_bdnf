library(data.table)
library(limma)
library(ggrepel)
#read in the data
phos <- as.data.table(janitor::clean_names(data.table::fread("data/Phospho_STYSites_interactive_processed4.csv")))
#excel added rows to the bottom as it does, get rid of those
phos = phos[!is.na(t)]

# one hour phos differential ----------------------------------------------

onehour_phos = phos %>% 
    dplyr::select(gene_names,amino_acid,phospho_sty_probabilities,protein,position,a1:b3) %>% #just take the control v bdnf1hr which are A1>B3,
    mutate(un = glue::glue("{gene_names}|{amino_acid}|{position}|{phospho_sty_probabilities}|{protein}")) %>% #turn the row names into something useful
    tibble::column_to_rownames("un") %>% #and into a rowname
    dplyr::select(-(gene_names:position)) %>% #bye-bye to those unneeded columns
    mutate(log(across(a1:b3),base = 2))

#Make a factor of conditions, set control to the first level
group <- factor(c("control","control","control",
                  "bdnf_1hr","bdnf_1hr","bdnf_1hr"),
                levels = c("control","bdnf_1hr"))

#????? this is from limma vignette, should we use a 0 here?
design <-model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
#????? this is from limma vignette, should it really be this "minus"?
contr.matrix <- makeContrasts(
    ControlvBDNF1hr = control- bdnf_1hr,
    levels = colnames(design))

#the rest just your pseudocode filled in
vfit <- lmFit(onehour_phos, design)
fit2 <- contrasts.fit(vfit, contr.matrix)
efit2 <- eBayes(fit2)

bdnf1hr_v_control <- topTreat(efit2, coef=1, adjust="fdr",sort.by="logFC",number=50000)

bdnf1hr_v_control = bdnf1hr_v_control %>% 
    tibble::rownames_to_column('un') %>% 
    separate(un, into =  c("gene_names","amino_acid","site","phospho_sty_probabilities","uniprot_id"),sep = "\\|") %>% 
    janitor::clean_names() %>% 
    as.data.table() %>% 
    mutate(log_fc = -log_fc) %>% 
    separate_rows(gene_names) %>% as.data.table() 


#
# six hour phos differential ----------------------------------------------

six_phos = phos %>% 
    dplyr::select(gene_names,amino_acid,phospho_sty_probabilities,protein,position,a1:a3,c1:c3) %>% #just take the control v bdnf1hr which are A1>B3,
    mutate(un = glue::glue("{gene_names}|{amino_acid}|{position}|{phospho_sty_probabilities}|{protein}")) %>% #turn the row names into something useful
    tibble::column_to_rownames("un") %>% #and into a rowname
    dplyr::select(-(gene_names:position)) %>% #bye-bye to those unneeded columns
    mutate(log(across(a1:c3),base = 2))


#Make a factor of conditions, set control to the first level
group <- factor(c("control","control","control",
                  "bdnf_6hr","bdnf_6hr","bdnf_6hr"),
                levels = c("control","bdnf_6hr"))

#????? this is from limma vignette, should we use a 0 here?
design <-model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
#????? this is from limma vignette, should it really be this "minus"?
contr.matrix <- makeContrasts(
    ControlvBDNF6hr = control- bdnf_6hr,
    levels = colnames(design))

#the rest just your pseudocode filled in
vfit <- lmFit(six_phos, design)
fit2 <- contrasts.fit(vfit, contr.matrix)
efit2 <- eBayes(fit2)

bdnf6hr_v_control <- topTreat(efit2, coef=1, adjust="fdr",sort.by="logFC",number=50000)


bdnf6hr_v_control = bdnf6hr_v_control %>% 
    tibble::rownames_to_column('un') %>% 
    separate(un, into =  c("gene_names","amino_acid","site","phospho_sty_probabilities","uniprot_id"),sep = "\\|") %>% 
    janitor::clean_names() %>% 
    as.data.table() %>% 
    mutate(log_fc = -log_fc) %>% 
    separate_rows(gene_names) %>% as.data.table()



bdnf1hr_v_control = fread('bdnf1hr_v_control_phos.csv')
bdnf6hr_v_control = fread('bdnf6hr_v_control_phos.csv')
# outputting files for phosphosite -------------------------------------
print("Seqs hour 1")
seqs = UniprotR::GetSequences(bdnf1hr_v_control$uniprot_id)$Sequence
bdnf1hr_v_control$seq = seqs
bdnf1hr_v_control = bdnf1hr_v_control |> 
    mutate(uniprot_flank = stringr::str_sub(seq,as.numeric(site) - 2, as.numeric(site) + 1)) |> 
    mutate(contains_pxstp_domain = grepl("P[A-Z][ST]P",uniprot_flank))

print("Seqs hour 6")
seqs6 = UniprotR::GetSequences(bdnf6hr_v_control$uniprot_id)$Sequence

bdnf6hr_v_control$seq = seqs6
bdnf6hr_v_control = bdnf6hr_v_control |> 
    mutate(uniprot_flank = stringr::str_sub(seq,as.numeric(site) - 2, as.numeric(site) + 1)) |> 
    mutate(contains_pxstp_domain = grepl("P[A-Z][ST]P",uniprot_flank))



bdnf1hr_v_control %>% 
    fwrite("bdnf1hr_v_control_phos.csv")

bdnf6hr_v_control %>% 
    fwrite("bdnf6hr_v_control_phos.csv")


# write out the sites for DE analysis -------------------------------------
window_length = 7
bdnf1hr_v_control |> 
    mutate(uniprot_flank = stringr::str_sub(seq,as.numeric(site) - window_length, as.numeric(site) + window_length)) |> 
    select(uniprot_flank,log_fc,p_value) |> 
    unique() |> 
    fwrite('for_de_test_phosphosite_one_pval.tsv',sep = '\t')
    
bdnf6hr_v_control |> 
    mutate(uniprot_flank = stringr::str_sub(seq,as.numeric(site) - window_length, as.numeric(site) + window_length)) |> 
    select(uniprot_flank,log_fc,p_value) |> 
    unique() |> 
    fwrite('for_de_test_phosphosite_six_pval.tsv',sep = '\t')
library(ggrepel)

one_hour_pps_de = fread('enrichment-analysis-result-table_onehour_pvalue.txt')

one_hour_pps_de |> 
    mutate(plot_label = ifelse(dominant_p_value < 0.05,kinase,NA_character_)) |> 
    mutate(plot_alpha = ifelse(dominant_p_value < 0.05,'sig','ns')) |> 
    filter(dominant_adjusted_p_value != 1) |> 
    ggplot(aes(x = dominant_enrichment_value_log2,
               y = -log10(dominant_p_value))) + 
    geom_point(aes(alpha = plot_alpha),show.legend = FALSE) + 
    geom_hline(linetype = 2, yintercept = -log10(0.1)) +
    geom_vline(xintercept = 0) +
    geom_label_repel(aes(label = plot_label),
                     max.overlaps = 30) +
    ggpubr::theme_pubr() + 
    scale_y_continuous(trans = scales::pseudo_log_trans()) +
    ylab(bquote('-Log'[10]~ 'enrichment p-value')) +
    xlab(bquote('Log'[2]~ 'kinase enrichment')) 


six_hour_pps_de = fread('enrichment-analysis-result-table_sixhour_pvalue.txt')

six_hour_pps_de |> 
    mutate(plot_label = ifelse(dominant_p_value < 0.05,kinase,NA_character_)) |> 
    mutate(plot_alpha = ifelse(dominant_p_value < 0.05,'sig','ns')) |> 
    filter(dominant_adjusted_p_value != 1) |> 
    ggplot(aes(x = dominant_enrichment_value_log2,
               y = -log10(dominant_p_value))) + 
    geom_point(aes(alpha = plot_alpha),show.legend = FALSE) + 
    geom_hline(linetype = 2, yintercept = -log10(0.1)) +
    geom_vline(xintercept = 0) +
    geom_label_repel(aes(label = plot_label),
                     max.overlaps = 30) +
    ggpubr::theme_pubr() + 
    scale_y_continuous(trans = scales::pseudo_log_trans()) +
    ylab(bquote('-Log'[10]~ 'enrichment p-value')) +
    xlab(bquote('Log'[2]~ 'kinase enrichment')) 

