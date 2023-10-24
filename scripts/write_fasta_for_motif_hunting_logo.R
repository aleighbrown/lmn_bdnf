library(data.table)
library(ggseqlogo)
library(tidyverse)

conflicted::conflict_prefer('select','dplyr')
conflicted::conflicts_prefer(dplyr::filter)
bdnf1hr_v_control = 
    fread("data/bdnf1hr_v_control_phos.csv")

bdnf6hr_v_control = fread('data/bdnf6hr_v_control_phos.csv')

#write the log2Fold and p-value for the kinase DE analysis https://www.phosphosite.org/kinaseLibraryAction
bdnf1hr_v_control |> 
    mutate(uniprot_flank_larger = stringr::str_sub(seq,as.numeric(site) - 7, as.numeric(site) + 7)) |> 
    select(uniprot_flank_larger,log_fc,p_value)  |> 
    unique() |> 
    fwrite('data/phospho_bdnf_1hr_for_enrichment_psp.tsv',sep = '\t',col.names = FALSE)
#GO categories
bdnf1hr_kinase_enrich = fread('data/enrichment-analysis-result-table_onehour_pvalue.txt')
bdnf6hr_kinase_enrich = fread('data/enrichment-analysis-result-table_sixhour_pvalue.txt')

one_kin = GOfuncR::get_anno_categories(unique(bdnf1hr_kinase_enrich$kinase))
these_mapk = one_kin |> filter(grepl("MAPK",name)) |> pull(gene) |> unique()
top_p = bdnf1hr_kinase_enrich |> 
    slice_max(n = 5,dominant_enrichment_value_log2) |> 
    pull(kinase) |> unique()
#re-make that plot

bdnf1hr_kinase_enrich |> 
    mutate(plotAlpha = dominant_p_value < 0.01) |> 
    mutate(plotKinase = ifelse(dominant_p_value < 0.01,kinase,NA_character_)) |> 
    # mutate(plotKinase = ifelse(dominant_p_value < 0.01 & (kinase %in% these_mapk| kinase %in% top_p),kinase,NA_character_)) |> 
    filter(!(dominant_direction == 'upregulated set' & dominant_p_value < 0.1 & dominant_enrichment_value_log2 < 0)) |> 
    mutate(plotFill = case_when(dominant_p_value >0.01 ~ 'grey',
                                dominant_direction == 'upregulated set' & dominant_p_value < 0.1 ~ 'red',
                                dominant_direction == 'downregulated set'~ 'blue')) |> 
    ggplot(aes(x = dominant_enrichment_value_log2,
               y = dominant_p_value_log10_abs,
               alpha = plotAlpha)) + 
    geom_point(pch = 21,fill = 'black') +
    theme_minimal() + 
    geom_hline(yintercept = 2,linetype = 2) + 
    geom_vline(xintercept = 0) + 
    scale_fill_manual(values = c("#446dea",'grey','red')) +
    scale_y_continuous(trans = scales::pseudo_log_trans()) + 
    ggrepel::geom_label_repel(aes(label = plotKinase),min.segment.length = 0,max.overlaps = 100) + 
    xlab(expression(paste("Lo", g[2], " kinase enrichment"))) +
    ylab(expression(paste("-Lo", g[10], " enrichment p-value"))) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none') 

bdnf6hr_kinase_enrich |> 
    mutate(plotAlpha = dominant_p_value < 0.01) |> 
    mutate(plotKinase = ifelse(dominant_p_value < 0.01,kinase,NA_character_)) |> 
    # mutate(plotKinase = ifelse(dominant_p_value < 0.01 & (kinase %in% these_mapk| kinase %in% top_p),kinase,NA_character_)) |> 
    filter(!(dominant_direction == 'upregulated set' & dominant_p_value < 0.1 & dominant_enrichment_value_log2 < 0)) |> 
    mutate(plotFill = case_when(dominant_p_value >0.01 ~ 'grey',
                                dominant_direction == 'upregulated set' & dominant_p_value < 0.1 ~ 'red',
                                dominant_direction == 'downregulated set'~ 'blue')) |> 
    ggplot(aes(x = dominant_enrichment_value_log2,
               y = dominant_p_value_log10_abs,
               alpha = plotAlpha)) + 
    geom_point(pch = 21,fill = 'black') +
    theme_minimal() + 
    geom_hline(yintercept = 2,linetype = 2) + 
    geom_vline(xintercept = 0) + 
    scale_fill_manual(values = c("#446dea",'grey','red')) +
    scale_y_continuous(trans = scales::pseudo_log_trans()) + 
    ggrepel::geom_label_repel(aes(label = plotKinase),min.segment.length = 0,max.overlaps = 100) + 
    xlab(expression(paste("Lo", g[2], " kinase enrichment"))) +
    ylab(expression(paste("-Lo", g[10], " enrichment p-value"))) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none') 

sig_to_write_1hr = bdnf1hr_v_control |> 
    mutate(uniprot_flank_larger = stringr::str_sub(seq,as.numeric(site) - 5, as.numeric(site) + 5)) |> 
    # mutate(phos_flank = glue::glue("{stringr::str_sub(seq,as.numeric(site) - 5,as.numeric(site))}{tolower(stringr::str_sub(seq,as.numeric(site),as.numeric(site)))}*{stringr::str_sub(seq,as.numeric(site)+1,as.numeric(site) +5)}")) |> 
    filter(p_value < 0.05) |>
    unique() |> 
    mutate(fasta_name = glue::glue(">{gene_names}|{amino_acid}{site}|{uniprot_id}")) |> 
    select(gene_names,fasta_name,uniprot_flank_larger,p_value,log_fc,amino_acid)

sig_to_write_6hr = bdnf6hr_v_control |> 
    mutate(uniprot_flank_larger = stringr::str_sub(seq,as.numeric(site) - 5, as.numeric(site) + 5)) |> 
    # mutate(phos_flank = glue::glue("{stringr::str_sub(seq,as.numeric(site) - 5,as.numeric(site))}{tolower(stringr::str_sub(seq,as.numeric(site),as.numeric(site)))}*{stringr::str_sub(seq,as.numeric(site)+1,as.numeric(site) +5)}")) |> 
    filter(p_value < 0.05) |>
    unique() |> 
    mutate(fasta_name = glue::glue(">{gene_names}|{amino_acid}{site}|{uniprot_id}")) |> 
    select(gene_names,fasta_name,uniprot_flank_larger,p_value,log_fc,amino_acid)

letts_found = str_split(sig_to_write_1hr$uniprot_flank_larger,"") |> 
    purrr::simplify() |> unique()

letts_found = letts_found[!is.na(letts_found)]

cs2 = make_col_scheme(chars=letts_found, values=1:20)
adj_u_1hr = sig_to_write_1hr |> 
    filter(nchar(uniprot_flank_larger) == 11) |> 
    # filter(gene_names %in% make_motif) |> 
    pull(uniprot_flank_larger) |> unique()

adj_u_6hr = sig_to_write_6hr |> 
    filter(nchar(uniprot_flank_larger) == 11) |> 
    # filter(gene_names %in% make_motif) |> 
    pull(uniprot_flank_larger) |> unique()

#MAKE sequence LOGO 1 hour
ggplot() + geom_logo( adj_u_1hr, method = 'prob',col_scheme = cs2, high_col = "black",low_col = 'black') + theme_logo()
#MAKE sequence LOGO 6 hour
ggplot() + geom_logo( adj_u_6hr, method = 'prob',col_scheme = cs2, high_col = "black",low_col = 'black') + theme_logo()


f <- file(description = 'onehour_pvalue_sig_5aa_updown.fasta', open = "a")

written = c()
for (i in 1 : nrow(sig_to_write)){
    header = sig_to_write[i,fasta_name]
    tail = sig_to_write[i,uniprot_flank_larger]
    if(!(tail %in% written)){
        if(!is.na(tail)){
            if(nchar(tail) > 6){
                writeLines(header, f)
                writeLines(tail, f)
            }
        }
        
    }
    
    written = c(tail,written)
    
    
}
close(f)


make_motif = c("PALLD", "KIF21A", "DMTN", "MAPRE1", "MAP2", "MAP7D1", "MAP1A", 
                "MAP1B", "MACF1", "IQGAP2", "STMN1", "CAMSAP3", "MYO5A", "MLLT4", 
                "MICAL3", "MPRIP", "CLASP2", "TSC2", "NTRK2", "MAP2K1", "RANBP2", 
                "MAPK1","DYNC1LI1", "DYNC1I1", "DYNC1LI1", "DYNC1I2")


f <- file(description = 'force_motif_adj.txt', open = "a")
written = c()
sig_to_write 
for (i in 1 : nrow(sig_to_write)){
    header = sig_to_write[i,fasta_name]
    tail = sig_to_write[i,uniprot_flank_larger]
    if(!(tail %in% written)){
        if(!is.na(tail)){
            if(nchar(tail) ==11){
                # writeLines(header, f)
                writeLines(tail, f)
            }
        }
        
    }
    
    written = c(tail,written)
    
    
}
close(f)




# find the akt motif ------------------------------------------------------

akt_pat = "R.{1}R.{2}[ST]" #RXRXXS/T
px_pat = "P.{1}[ST]P" #RXRXXS/T

bdnf1hr_v_control |> 
    mutate(uniprot_flank_akt = stringr::str_sub(seq,as.numeric(site) - 5, as.numeric(site))) |> 
    mutate(uniprot_flank_pxsp = stringr::str_sub(seq,as.numeric(site) - 2, as.numeric(site) + 1)) |>
    select(gene_names, uniprot_flank_akt,uniprot_flank_pxsp,site,amino_acid,p_value,log_fc) |> 
    mutate(rxrxxST = grepl(akt_pat,uniprot_flank_akt)) |> 
    mutate(pxSTp = grepl(px_pat,uniprot_flank_pxsp)) |>
    mutate(STp = grepl("[ST]P",str_sub(uniprot_flank_pxsp, start= -2))) |> 
    fwrite("~/Desktop/bdnf1hr_v_control_with_motif.csv")
    # filter(gene_names == "STAU2")
bdnf6hr_v_control |> 
    mutate(uniprot_flank_akt = stringr::str_sub(seq,as.numeric(site) - 5, as.numeric(site))) |> 
    mutate(uniprot_flank_pxsp = stringr::str_sub(seq,as.numeric(site) - 2, as.numeric(site) + 1)) |>
    select(gene_names, uniprot_flank_akt,uniprot_flank_pxsp,site,amino_acid,p_value,log_fc) |> 
    mutate(rxrxxST = grepl(akt_pat,uniprot_flank_akt)) |> 
    mutate(pxSTp = grepl(px_pat,uniprot_flank_pxsp)) |>
    mutate(STp = grepl("[ST]P",str_sub(uniprot_flank_pxsp, start= -2))) |> 
    fwrite("~/Desktop/bdnf6hr_v_control_with_motif.csv")
