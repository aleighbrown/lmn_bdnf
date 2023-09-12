



splicing_dots_tables_function <- function(input_splicing, 
                                          list_gene = "",
                                          psi_cutoff = 0,
                                          probability_changing_cutoff = 0,
                                          upper_bound = 4.5,
                                          prob_cut = 0.9) {
    
    if(psi_cutoff != 0 ){
        psi_genes <- input_splicing %>% 
            filter(mean_dpsi_per_lsv_junction >= psi_cutoff) %>% 
            filter(probability_changing >= 0.9)
    }
    
    if(probability_changing_cutoff != 0){
        psi_genes <- input_splicing %>% 
            filter(mean_dpsi_per_lsv_junction >= psi_cutoff) %>% 
            filter(probability_changing >= probability_changing_cutoff)
    }
    
    if(exists('psi_genes')){
        
        psi_genes = psi_genes %>% pull(gene_name)
        list_gene = unique(c(list_gene,psi_genes))
    }
    
    splicing_dots_tables <- input_splicing %>% 
        mutate(junction_name = case_when(gene_name %in% list_gene &
                                             mean_dpsi_per_lsv_junction > 0.2 & 
                                             probability_changing > prob_cut ~ gene_name,
                                         T ~ "")) %>%
        mutate(log10_test_stat = -log10(1 - probability_changing)) %>%
        mutate(log10_test_stat = ifelse(is.infinite(log10_test_stat), upper_bound, log10_test_stat)) %>%
        mutate(graph_alpha = ifelse(probability_changing > 0.9, 1, 0.2)) %>%
        mutate(label_junction = case_when(gene_name %in% list_gene &
                                              mean_dpsi_per_lsv_junction > 0.2 &
                                              probability_changing > prob_cut ~ junction_name,
                                          T ~ NA_character_))
    

    fig1 <- ggplot() +
        ggrastr::geom_point_rast(data = splicing_dots_tables %>% dplyr::filter(splicing_location == "CDS"),
                   aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat,
                       alpha = graph_alpha), fill = "#4E3822", pch = 21, size = 2) +
        ggrastr::geom_point_rast(data = splicing_dots_tables %>% dplyr::filter(splicing_location == "3'UTR"),
                   aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat,
                       alpha = graph_alpha), fill = "#01FDF6", pch = 21, size = 2) +
        ggrastr::geom_point_rast(data = splicing_dots_tables %>% dplyr::filter(in_five == TRUE),
                   aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat,
                       alpha = graph_alpha), fill = "#E3D26F", pch = 21, size = 2) +
        geom_text_repel(data = splicing_dots_tables[!is.na(label_junction)],
                        aes(x = mean_dpsi_per_lsv_junction, y =log10_test_stat,
                            label = label_junction), 
                        point.padding = 0.3,
                        nudge_y = 0.2,
                        min.segment.length = 0,
                        box.padding  = 2,
                        max.overlaps = 10000,
                        size=4, 
                        show.legend = F) +
        geom_hline(yintercept = -log10(1 - .9)) +
        geom_vline(xintercept = -0,linetype="dotted") +
        guides(alpha = "none", size = "none") +
        theme(legend.position = 'top') +
        ggpubr::theme_pubr() +
        xlab("delta PSI") +
        ylab(expression(paste("-Lo", g[10], " Test Statistic"))) +
        scale_x_continuous(labels = scales::percent, limits = c(-1,1))
    return(fig1)
}

return_the_location_of_junction = function(my_grange){
    
    this_table = my_grange |> 
        annotate_regions(annotations = annotations) |> 
        as.data.table() |> 
        mutate(tx = gsub("\\..*", "", annot.tx_id)) 
    
    this_table = this_table |> 
        filter(tx %in% mane_transcripts) |> 
        filter(!(annot.type %in% c("hg38_genes_promoters","hg38_genes_1to5kb"))) |> 
        filter(gene_name == annot.symbol) |> 
        select(paste_into_igv_junction,gene_name,mean_dpsi_per_lsv_junction:lsv_type,annot.type,tx,annot.tx_id) |> 
        select(-lsv_type) |> 
        unique() |> 
        filter(annot.type %in% c("hg38_genes_5UTRs","hg38_genes_3UTRs"))  
        
    
    my_df = my_grange |> 
        mutate(in_five = (paste_into_igv_junction %in% this_table[annot.type == "hg38_genes_5UTRs"]$paste_into_igv_junction)) |>
        mutate(in_three = (paste_into_igv_junction %in% this_table[annot.type == "hg38_genes_3UTRs"]$paste_into_igv_junction)) |>
        as.data.table() |> 
        select(paste_into_igv_junction,seqnames:lsv_type,in_five,in_three) |> 
        select(-lsv_type,-lsv_id)
    # 
    return(my_df)
}

if (!require("pacman")) install.packages("pacman")
library(pacman)
library(data.table)
library(GenomicFeatures)
library(tidyverse)
pacman::p_load(annotatr,janitor,tidyverse,AnnotationDbi,org.Hs.eg.db,data.table,ggrepel)
library(plyranges)
z <- makeTxDbFromGFF("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.1.ensembl_genomic.gff.gz")
mane_transcripts = transcripts(z)
mane_transcripts = mane_transcripts |> 
    mutate(tx = gsub("\\..*", "", tx_name)) |> 
    as.data.table() |> 
    pull(tx) |> 
    unique()

#build the annotations for the hg38 basic genes set using the annotatr package
annotations = build_annotations(genome = "hg38", annotations = c("hg38_basicgenes"))
annotations = keepStandardChromosomes(annotations, pruning.mode = "coarse")

#find the significant junctions
hour_one = fread("/Users/annaleigh/Desktop/bdnf-4su_splicing/control1-bdnf1_annotated_junctions.csv")
hour_two = fread("/Users/annaleigh/Desktop/bdnf-4su_splicing/control2-bdnf2_annotated_junctions.csv")
hour_six = fread("/Users/annaleigh/Desktop/bdnf-4su_splicing/control6-bdnf6_annotated_junctions.csv")


junctions_hour_one = hour_one|> 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) 
junctions_hour_two = hour_two|> 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) 
junctions_hour_six = hour_six|> 
    makeGRangesFromDataFrame(,keep.extra.columns = TRUE) 


junctions_hour_one = return_the_location_of_junction(junctions_hour_one)
junctions_hour_one = junctions_hour_one |> 
    filter(!(in_five == TRUE & in_three == TRUE)) |> 
    mutate(splicing_location = case_when(in_five == TRUE ~ "5'UTR",
                                         in_three == TRUE ~ "3'UTR",
                                         TRUE ~ "CDS")) 
junctions_hour_two = return_the_location_of_junction(junctions_hour_two)
junctions_hour_two = junctions_hour_two |> 
    filter(!(in_five == TRUE & in_three == TRUE)) |> 
    mutate(splicing_location = case_when(in_five == TRUE ~ "5'UTR",
                                         in_three == TRUE ~ "3'UTR",
                                         TRUE ~ "CDS")) 
junctions_hour_six = return_the_location_of_junction(junctions_hour_six)
junctions_hour_six = junctions_hour_six |> 
    filter(!(in_five == TRUE & in_three == TRUE)) |> 
    mutate(splicing_location = case_when(in_five == TRUE ~ "5'UTR",
                                         in_three == TRUE ~ "3'UTR",
                                         TRUE ~ "CDS")) 



both_time_points_sig = junctions_hour_two[paste_into_igv_junction %in% 
                                              intersect(junctions_hour_two[probability_changing > 0.9,paste_into_igv_junction],
          junctions_hour_six[probability_changing > 0.9,paste_into_igv_junction])] |> 
    select(gene_name,splicing_location,paste_into_igv_junction) |> unique()
          


# identify transcription factors with 5'utr splicing ----------------------
#29425488
tfs = fread("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt") |> janitor::clean_names()
tfs = tfs[is_tf == 'Yes']
total_tf |> 
    separate(pathway, into = c("tf","dir"),convert = TRUE,remove = FALSE) |> 
    dplyr::filter(confidence_levels == 'A') |> 
    select(tf,padj,NES,timepoint) |> 
    filter(padj < 0.1) |> 
    arrange(tf)

one_gnes = junctions_hour_one |> 
    mutate(gene_id = gsub("\\..*", "", gene_id)) |> 
    filter(gene_id %in% tfs$ensembl_id) |> 
    filter(!grepl("^CH[0-9]",gene_name)) |> 
    # filter(splicing_location == "5'UTR") |> 
    filter(!grepl("^AC0",gene_name)) |> 
    filter(abs(mean_dpsi_per_lsv_junction) > 0.2) |> 
    select(splicing_location,mean_dpsi_per_lsv_junction,bdnf1_mean_psi,control1_mean_psi,probability_changing,gene_name,paste_into_igv_junction) 

two_gnes = junctions_hour_two |> 
    mutate(gene_id = gsub("\\..*", "", gene_id)) |> 
    filter(gene_id %in% tfs$ensembl_id) |> 
    filter(!grepl("^CH[0-9]",gene_name)) |> 
    # filter(splicing_location == "5'UTR") |> 
    filter(!grepl("^AC0",gene_name)) |> 
    filter(abs(mean_dpsi_per_lsv_junction) > 0.2) |> 
    select(splicing_location,mean_dpsi_per_lsv_junction,bdnf2_mean_psi,control2_mean_psi,probability_changing,gene_name,paste_into_igv_junction)


six_gnes = junctions_hour_six |> 
    mutate(gene_id = gsub("\\..*", "", gene_id)) |> 
    filter(gene_id %in% tfs$ensembl_id) |> 
    filter(!grepl("^CH[0-9]",gene_name)) |> 
    # filter(splicing_location == "5'UTR") |> 
    filter(!grepl("^AC0",gene_name)) |> 
    filter(abs(mean_dpsi_per_lsv_junction) > 0.2) |> 
    select(splicing_location,mean_dpsi_per_lsv_junction,bdnf6_mean_psi,control6_mean_psi,probability_changing,gene_name,paste_into_igv_junction)

 
six_gnes[splicing_location == "5'UTR"] |> 
    pull(gene_name) |> 
    unique()

splicing_dots_tables_function(junctions_hour_one,
                              c("NTRK2","ATF3"))

splicing_dots_tables_function(junctions_hour_two,
                              c("AKT1","ATF3"))

splicing_dots_tables_function(junctions_hour_six,
                              c("ATF3"))

ggpubr::ggarrange(plotlist = list(hour1,hour2,hour6),ncol = 3)
