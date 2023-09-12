library(org.Hs.eg.db)
library(tidyverse)
library(tidyr)
library(dplyr)
library(data.table)
library(janitor)
library(limma)
library(ggrepel)

new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv") 
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de %>% 
    dplyr::select(-gene_name) %>% 
    mutate(gene = gsub("\\..*", "", gene)) %>%  
    left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol), by = c("gene" = "ensgene")) %>% 
    unique()

# Reading in and cleaning 1 hour
onehour_total = as.data.table(clean_names(fread("data/proteome/total_proteome_bdnf1hr_control.csv")))
onehour_total = onehour_total %>% 
    tidyr::separate_rows(v1)

conversion = gprofiler2::gconvert(onehour_total$v1)

onehour_total = onehour_total %>% 
    left_join(conversion %>% dplyr::select(input,name,target),
              by = c("v1" = "input"))


onehour_total = onehour_total %>% 
    mutate(de = case_when(p_value > 0.05 ~ 'not sig',
                          p_value <= 0.05 & log_fc > 0 ~ "upregulated",
                          p_value <= 0.05 & log_fc < 0 ~ "downregulated"))

# Readign in and cleanign 6 hour
sixhour_total = as.data.table(clean_names(fread("data/proteome/total_proteome_bdnf6hr_control.csv")))
sixhour_total = sixhour_total %>% 
    tidyr::separate_rows(v1)

conversion = gprofiler2::gconvert(sixhour_total$v1)

sixhour_total = sixhour_total %>% 
    left_join(conversion %>% dplyr::select(input,name,target),
              by = c("v1" = "input"))

sixhour_total = sixhour_total %>% 
    mutate(de = case_when(p_value > 0.05 ~ 'not sig',
                          p_value <= 0.05 & log_fc > 0 ~ "upregulated",
                          p_value <= 0.05 & log_fc < 0 ~ "downregulated"))

# GO analysis on one hours phososites -------------------------------------------
#1h up
bdnf1hr_v_control = fread('bdnf1hr_v_control_phos.csv')

sig_up_one <- bdnf1hr_v_control %>% filter(adj_p_val < 0.1 &
                                              log_fc > 0) %>% 
    pull(gene_names) %>% unique() %>% 
    gprofiler2::gconvert() %>% 
    pull(target)

sig_down_one <- bdnf1hr_v_control %>% filter(adj_p_val < 0.1 &
                                               log_fc < 0) %>% 
    pull(gene_names) %>% unique() %>% 
    gprofiler2::gconvert() %>% 
    pull(target)


up_reg_phos_1hr = clusterProfiler::enrichGO(sig_up_one,
                                            universe = unique(new_ratio_bayesian_p_de$gene),
                                            keyType = 'ENSEMBL',
                                            OrgDb = org.Hs.eg.db,
                                            ont = 'ALL',
                                            readable = TRUE)

down_reg_phos_1hr = clusterProfiler::enrichGO(sig_down_one,
                                              universe = unique(new_ratio_bayesian_p_de$gene),
                                              keyType = 'ENSEMBL',
                                              OrgDb = org.Hs.eg.db,
                                              ont = 'ALL',
                                              readable = TRUE)

sig_phos_1hr = clusterProfiler::enrichGO(unique(c(sig_down_one,sig_up_one)),
                                         # universe = unique(unique(new_ratio_bayesian_p_de$gene_name)),
                                         keyType = 'ENSEMBL',
                                         OrgDb = org.Hs.eg.db,
                                         ont = 'ALL',
                                         readable = TRUE)

clusterProfiler::dotplot(up_reg_phos_1hr) + ggtitle("BDNF 1hr\nUp-Phosoproteome")
clusterProfiler::cnetplot(up_reg_phos_1hr) + ggtitle("BDNF 1hr\nUp-Phosoproteome")

clusterProfiler::dotplot(down_reg_phos_1hr) + ggtitle("BDNF 1hr\nDown-Phosoproteome")
clusterProfiler::cnetplot(down_reg_phos_1hr) + ggtitle("BDNF 1hr\nDown-Phosoproteome")

clusterProfiler::dotplot(sig_phos_1hr) + ggtitle("BDNF 1hr\nSig DifferentPhosoproteome")
clusterProfiler::cnetplot(sig_phos_1hr) + ggtitle("BDNF 1hr\nSig DifferentPhosoproteome")


# GO analyais on six hours phos -------------------------------------------
bdnf6hr_v_control = fread('bdnf6hr_v_control_phos.csv')

sig_up_six <- bdnf6hr_v_control %>% filter(adj_p_val < 0.1 &
                                               log_fc > 0) %>% 
    pull(gene_names) %>% unique() %>% 
    gprofiler2::gconvert() %>% 
    pull(target)

sig_down_six <- bdnf6hr_v_control %>% filter(adj_p_val < 0.1 &
                                                 log_fc < 0) %>% 
    pull(gene_names) %>% unique() %>% 
    gprofiler2::gconvert() %>% 
    pull(target)

up_reg_phos_6hr = clusterProfiler::enrichGO(sig_up_six,
                                            universe = unique(new_ratio_bayesian_p_de$gene),
                                            keyType = 'ENSEMBL',
                                            OrgDb = org.Hs.eg.db,
                                            ont = 'ALL',
                                            readable = TRUE)

down_reg_phos_6hr = clusterProfiler::enrichGO(sig_down_six,
                                              universe = unique(new_ratio_bayesian_p_de$gene),
                                              keyType = 'ENSEMBL',
                                              OrgDb = org.Hs.eg.db,
                                              ont = 'ALL',
                                              readable = TRUE)

sig_phos_6hr = clusterProfiler::enrichGO(unique(c(sig_up_six,sig_down_six)),
                                         #universe = unique(unique(new_ratio_bayesian_p_de$gene_name)),
                                         keyType = 'ENSEMBL',
                                         OrgDb = org.Hs.eg.db,
                                         ont = 'ALL',
                                         readable = TRUE)

clusterProfiler::dotplot(up_reg_phos_6hr) + ggtitle("BDNF 6hr\nUp-Phosoproteome")
clusterProfiler::cnetplot(up_reg_phos_6hr) + ggtitle("BDNF 6hr\nUp-Phosoproteome")

clusterProfiler::dotplot(down_reg_phos_6hr) + ggtitle("BDNF 6hr\nDown-Phosoproteome")
clusterProfiler::cnetplot(down_reg_phos_6hr) + ggtitle("BDNF 6hr\nDown-Phosoproteome")

clusterProfiler::dotplot(sig_phos_6hr) + ggtitle("BDNF 6hr\nSig DifferentPhosoproteome")
clusterProfiler::cnetplot(sig_phos_6hr) + ggtitle("BDNF 6hr\nSig DifferentPhosoproteome")




# bdnf1hr_v_control = fread(here::here('bdnf1hr_v_control_phos.csv'))
# bdnf6hr_v_control = fread("bdnf6hr_v_control_phos.csv")
# volcano plot hilighting microtubles 1 hour -------------------------------------
# microtubles = janitor::clean_names(readxl::read_excel('/Users/annaleigh/Documents/GitHub/bdnf_4su/data/Microtubule list_Mass spec.xlsx'))

these_genes = c("PALLD", "KIF21A", "DMTN", "MAPRE1", "MAP2", "MAP7D1", "MAP1A", 
  "MAP1B", "MACF1", "IQGAP2", "STMN1", "CAMSAP3", "MYO5A", "MLLT4", 
  "MICAL3", "MPRIP", "CLASP1", "CLASP2", "TSC2", "NTRK2", "MAP2K1", "RANBP2", 
  "MAPK1","DCX","LIN28A")

results_edit = bdnf1hr_v_control %>% 
    mutate(site = glue::glue("{gene_names}-{amino_acid}{site}")) %>% 
    mutate(plotAlpha = ifelse( adj_p_val < 0.1,0.4,0.2)) %>% 
    mutate(plotAlpha = ifelse( gene_names %in% these_genes,0.8,plotAlpha)) %>% 
    mutate(plotColor = case_when(gene_names %in% these_genes  ~ "#8baac4",
                              TRUE ~ '#ffffff')) %>% 
    mutate(plotSize = ifelse(gene_names %in% these_genes
                          ,3,
                          1.5)) %>% 
    group_by(gene_names) %>%
    mutate(max_change = max(abs(log_fc))) %>%
    ungroup() %>%
    mutate(plotlabel = case_when(gene_names =='LIN28A' & abs(log_fc) == max_change ~ site,
                                 adj_p_val < 0.1 & gene_names %in% these_genes ~ site,
                                 adj_p_val >= 0.1 ~ NA_character_,
                                 gene_names %in% these_genes & abs(log_fc) == max_change ~ site,
                                 TRUE ~ NA_character_))
    

ggplot(results_edit,aes(x = log_fc, y = -log10(adj_p_val))) + 
    geom_point(aes(alpha = plotAlpha),show_guide  = F,pch = 21,
               fill = results_edit$plotColor, size = results_edit$plotSize,stroke = 0.6) + 
    geom_hline(yintercept = -log10(0.1),linetype="dotted") + 
    geom_vline(xintercept = -0.5,linetype="dotted") + 
    geom_vline(xintercept = 0.5,linetype="dotted") + 
    ggpubr::theme_classic2() + 
    geom_label_repel(aes(label = plotlabel),size = 5,min.segment.length = 0) +
    ylab(bquote('-Log'[10]~ 'Adjusted p-value')) + 
    xlab(bquote('Log'[2]~ 'Fold Change')) +
    theme(text = element_text(size = 25))



# volcano plot higlighting top genes 1 hr ---------------------------------
top_n = bdnf1hr_v_control %>% slice_max(abs(log_fc),n = 20) %>% pull(gene_names)
top_n = bdnf1hr_v_control %>% filter(log_fc > 1.8) %>% pull(gene_names)

re_top_1 = bdnf1hr_v_control %>% 
    mutate(site = glue::glue("{gene_names}-{amino_acid}{site}")) %>% 
    mutate(plotAlpha = ifelse( adj_p_val < 0.1,1,0.2)) %>% 
    mutate(plotColor = case_when(gene_names %in% top_n & adj_p_val < 0.1 ~ "#cf793e",
                                 TRUE ~ '#808080')) %>% 
    mutate(plotSize = ifelse(gene_names %in% top_n
                             ,3,
                             1.5)) %>% 
    group_by(gene_names) %>%
    mutate(max_change = max(abs(log_fc))) %>%
    ungroup() %>%
    mutate(plotlabel = case_when(adj_p_val < 0.1 & gene_names %in% top_n ~ site,
                                 gene_names %in% top_n & abs(log_fc) == max_change ~ site,
                                 TRUE ~ NA_character_))


ggplot(re_top_1,aes(x = log_fc, y = -log10(adj_p_val))) + 
    geom_point(aes(alpha = plotAlpha),show_guide  = F,pch = 21,
               fill = re_top_1$plotColor, size = re_top_1$plotSize ) + 
    geom_hline(yintercept = -log10(0.1),linetype="dotted") + 
    geom_vline(xintercept = -0.5,linetype="dotted") + 
    geom_vline(xintercept = 0.5,linetype="dotted") + 
    ggpubr::theme_classic2() + 
    geom_text_repel(min.segment.length = 0.1,aes(label = plotlabel),size = 2) +
    ylab(bquote('-Log'[10]~ 'Adjusted p-value')) + 
    xlab(bquote('Log'[2]~ 'Fold Change')) +
    ggtitle('1 hr BDNF Phosphoproteome') + 
    coord_cartesian(xlim=c(1,5))


# volcano plot higlighting DYNC 1 hr ---------------------------------
matr3 = c("MATR3", "YLPM1", "PPP1CA", "PPP1CC", "SLC2A1", "COPA", "COPE", "RPL14", "PTBP1", "RBM26", "RPL21", "PPP1CB", "PTBP3", "ZNF738", "RBM14", "SLC25A5", "SLC2A3", "HNRNPH1", "DDX5", "SLC25A6", "RPL15", "RPL7A", "DNMT3B", "RPS11", "RPL23", "RPL4", "DDX3X", "NKRF", "COPB2", "NAT10", "EIF4A1", "RPL11", "RPS8", "RPS9", "RPL27", "MCM7", "DDX21", "RPL18", "RPL32", "DDX17", "DDX39B", "OAT", "RPL3", "RPL10", "ADAR", "XRN2", "DRG1", "RPL10A", "RPS3", "RPL22", "RALY", "TCP1", "L1TD1", "ILF3", "ESPL1", "DHX9", "DHCR7", "NUP205", "IGF2BP1", "TUBA1A", "RPL28", "EIF4A3", "STT3A", "TRAP1", "HNRNPA0", "NONO", "RACK1", "DNAJB6", "RPS4X", "CCT3", "RPL5", "RPL7", "PCBP1", "RPL36", "RPS23", "GARS", "PCNA", "RPL6", "RPS3A", "HNRNPU", "PHB2", "FADS2", "TUBB3", "SEC61A1", "PSMD3", "EIF3CL", "RSL1D1", "ATP5F1A", "CCT6A", "RPL9", "RAN", "RPL23A", "KPNA2", "RPS2", "TUBAL3", "DNAJA1", "UQCRC2", "SRSF9", "EIF3D", "XPO1", "TUBB", "AASS", "RPS16", "RPL30", "HNRNPM", "EIF3F", "HNRNPC", "IGF2BP3", "RPS12", "RBMX", "PHB", "ATP5F1C", "FBL", "PRKDC", "ARCN1", "CAND1", "VDAC2", "TOP2A", "SUGP2", "RPL17", "LARS", "RPL19", "YBX1", "RPS14", "DDX1", "HSP90AB1", "RPS7", "RPLP2", "ALDH18A1", "ZNF326", "DARS", "TRIM21", "PSIP1", "HNRNPR", "TUFM", "RPL18A", "SERPINH1", "RPS27L", "PCBP2", "PSMC5", "ATP5F1B", "ACOT7", "HIST1H1B", "SRSF1", "IARS", "NPM1", "SF3B1", "RPS18", "HYOU1", "PSMC2", "SYNCRIP", "EEF1A1", "COL18A1", "RPL38", "LAMB1")
matr3 = bdnf1hr_v_control$gene_names[grep("^DYNC",bdnf1hr_v_control$gene_names)]
ma_dync = bdnf1hr_v_control %>% 
    mutate(site = glue::glue("{gene_names}-{amino_acid}{site}")) %>% 
    mutate(plotAlpha = ifelse(p_value < 0.05,1,0.2)) %>% 
    mutate(plotColor = case_when(gene_names == 'DYNC1LI1' ~ "#f7cc31",
                                 TRUE ~ '#808080')) %>% 
    mutate(plotSize = ifelse(gene_names == 'DYNC1LI1',
                             3,
                             1.5)) %>% 
    mutate(plotlabel = case_when(gene_names == 'DYNC1LI1' ~ site,
                                 TRUE ~ NA_character_)) 
    # mutate(plotlabel = case_when(adj_p_val < 0.1 & grepl("^DYNC",gene_names) ~ site,
    #                              grepl("^DYNC",gene_names) & abs(log_fc) == max_change ~ site,
    #                              TRUE ~ NA_character_))


ggplot(ma_dync,aes(x = log_fc, y = -log10(p_value))) + 
    geom_point(aes(alpha = plotAlpha),show_guide  = F,
               pch = 21,
               fill = ma_dync$plotColor, 
               size = ma_dync$plotSize ) + 
    geom_hline(yintercept = -log10(0.05),linetype="dotted") + 
    geom_vline(xintercept = -0.5,linetype="dotted") + 
    geom_vline(xintercept = 0.5,linetype="dotted") + 
    ggpubr::theme_classic2() + 
    geom_label_repel(min.segment.length = 0.1, box.padding = 1,aes(label = plotlabel),size = 4) +
    ylab(bquote('-Log'[10]~ 'P-value')) + 
    xlab(bquote('Log'[2]~ 'Fold Change')) +
    ggtitle('1 hr BDNF Phosphoproteome')





# volcano plot hilighting microtubles 6 hour -------------------------------------

these_genes = c("PALLD", "KIF21A", "DMTN", "MAPRE1", "MAP2", "MAP7D1", "MAP1A", 
                "MAP1B", "MACF1", "IQGAP2", "STMN1", "CAMSAP3", "MYO5A", "MLLT4", 
                "MICAL3", "MPRIP", "CLASP1", "CLASP2", "TSC2", "NTRK2", "MAP2K1", "RANBP2", 
                "MAPK1","DCX","LIN28A")

results_edit = bdnf6hr_v_control %>% 
    mutate(site = glue::glue("{gene_names}-{amino_acid}{site}")) %>% 
    mutate(plotAlpha = ifelse( adj_p_val < 0.1,0.4,0.2)) %>% 
    mutate(plotAlpha = ifelse( gene_names %in% these_genes,0.8,plotAlpha)) %>% 
    mutate(plotColor = case_when(gene_names %in% these_genes  ~ "#8baac4",
                                 TRUE ~ '#ffffff')) %>% 
    mutate(plotSize = ifelse(gene_names %in% these_genes
                             ,3,
                             1.5)) %>% 
    group_by(gene_names) %>%
    mutate(max_change = max(abs(log_fc))) %>%
    ungroup() %>%
    mutate(plotlabel = case_when(gene_names =='LIN28A' & abs(log_fc) == max_change ~ site,
        adj_p_val < 0.1 & gene_names %in% these_genes ~ site,
        adj_p_val >= 0.1 ~ NA_character_,
        gene_names %in% these_genes & abs(log_fc) == max_change ~ site,
        TRUE ~ NA_character_))


ggplot(results_edit,aes(x = log_fc, y = -log10(adj_p_val))) + 
    geom_point(aes(alpha = plotAlpha),show_guide  = F,pch = 21,
               fill = results_edit$plotColor, size = results_edit$plotSize,stroke = 0.6) + 
    geom_hline(yintercept = -log10(0.1),linetype="dotted") + 
    geom_vline(xintercept = -0.5,linetype="dotted") + 
    geom_vline(xintercept = 0.5,linetype="dotted") + 
    ggpubr::theme_classic2() + 
    geom_label_repel(aes(label = plotlabel),size = 5,min.segment.length = 0) +
    ylab(bquote('-Log'[10]~ 'Adjusted p-value')) + 
    xlab(bquote('Log'[2]~ 'Fold Change')) +
    theme(text = element_text(size = 25))


# volcano plot higlighting top genes 6 hr ---------------------------------
top_n = bdnf6hr_v_control %>% slice_max(abs(log_fc),n = 20) %>% pull(gene_names)
re_top_6 = bdnf6hr_v_control %>% 
    mutate(site = glue::glue("{gene_names}-{amino_acid}{site}")) %>% 
    mutate(plotAlpha = ifelse( adj_p_val < 0.1,1,0.2)) %>% 
    mutate(plotColor = case_when(gene_names %in% top_n & adj_p_val < 0.1 ~ "#cf793e",
                                 TRUE ~ '#808080')) %>% 
    mutate(plotSize = ifelse(gene_names %in% top_n
                             ,3,
                             1.5)) %>% 
    group_by(gene_names) %>%
    mutate(max_change = max(abs(log_fc))) %>%
    ungroup() %>%
    mutate(plotlabel = case_when(adj_p_val < 0.1 & gene_names %in% top_n ~ site,
                                 gene_names %in% top_n & abs(log_fc) == max_change ~ site,
                                 TRUE ~ NA_character_))


ggplot(re_top_6,aes(x = log_fc, y = -log10(adj_p_val))) + 
    geom_point(aes(alpha = plotAlpha),show_guide  = F,pch = 21,
               fill = re_top_1$plotColor, size = re_top_1$plotSize ) + 
    geom_hline(yintercept = -log10(0.1),linetype="dotted") + 
    geom_vline(xintercept = -0.5,linetype="dotted") + 
    geom_vline(xintercept = 0.5,linetype="dotted") + 
    ggpubr::theme_classic2() + 
    geom_label_repel(aes(label = plotlabel),size = 4) +
    ylab(bquote('-Log'[10]~ 'Adjusted p-value')) + 
    xlab(bquote('Log'[2]~ 'Fold Change')) +
    ggtitle('6 hr BDNF Phosphoproteome')



# onehour_phos = as.data.table(clean_names(fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/proteome/phospho_proteome_bdnf1hr_control.csv")))
# onehour_phos = onehour_phos %>% 
#     separate_rows(id)

# conversion = gprofiler2::gconvert(onehour_phos$id)
# 
# onehour_phos = onehour_phos %>% 
#     left_join(conversion %>% dplyr::select(input,name,target),
#               by = c("id" = "input")) %>% dplyr::select(-id) %>% 
#     unique()
# 
# onehour_phos = onehour_phos %>% 
#     mutate(de = case_when(adj_p_val >= 0.1 ~ 'not sig',
#                           adj_p_val < 0.1 & log_fc > 0 ~ "upregulated",
#                           adj_p_val < 0.1 & log_fc < 0 ~ "downregulated"))








up_reg_total_1hr = clusterProfiler::enrichGO(onehour_total %>% filter(p_value < 0.05 &
                                                                        log_fc > 0) 
                                            %>% pull(target),
                                            universe = unique(new_ratio_bayesian_p_de$gene),
                                            keyType = 'ENSEMBL',
                                            OrgDb = org.Hs.eg.db,
                                            ont = 'ALL',
                                            readable = TRUE)

down_reg_total_1hr = clusterProfiler::enrichGO(onehour_total %>% filter(p_value < 0.05 &
                                                                          log_fc < 0) 
                                              %>% pull(target),
                                              universe = unique(new_ratio_bayesian_p_de$gene),
                                              keyType = 'ENSEMBL',
                                              OrgDb = org.Hs.eg.db,
                                              ont = 'ALL',
                                              readable = TRUE)

clusterProfiler::dotplot(up_reg_total_1hr) + ggtitle("BDNF 1hr\nUp-Total proteome")
clusterProfiler::dotplot(down_reg_total_1hr) + ggtitle("BDNF 1hr\nDown-Total proteome")
clusterProfiler::cnetplot(down_reg_total_1hr)




onehourphos = onehour_phos %>% filter(p_value < 0.05) %>% 
    pull(name) %>% unique()


onehourphos_gocateogries = as.data.table(GOfuncR::get_anno_categories(onehourphos))


test = list(axonal_gene = axonal_genes_either,
            sigonehourphos = onehourphos)

ggVennDiagram(test) + 
    scale_color_manual(values = c("black","black","black"))

new_ratio_bayesian_p_de %>% 
    mutate(gene = gsub("\\..*","",gene))

disbois = new_ratio_bayesian_p_de %>% filter(time ==2 & 
                                       total_rna_sig == "downregulated" & 
                                       new_rna_sig == "bdnf_higher_new_rna") %>% 
    pull(gene)

tmp = clusterProfiler::enrichGO(disbois,
                          universe = unique(new_ratio_bayesian_p_de$gene),
                          keyType = 'ENSEMBL',
                          OrgDb = org.Hs.eg.db,
                          ont = 'BP',
                          readable = TRUE)

tmp@result %>% filter(p.adjust < 0.1) %>% 
    dplyr::select(Description,p.adjust,GeneRatio,geneID) %>% 
    separate_rows(geneID) %>% 
    separate(GeneRatio,into = c("num","denom"),remove = FALSE) %>% 
    mutate(strength = as.numeric(num)/as.numeric(denom)) %>% 
    filter(grepl("splicing",Description)) %>% fwrite("~/Desktop/downregulated_twohours_plot_these_genes.csv")





bdnf1hr_v_control = bdnf1hr_v_control |> 
    left_join(mapped_bdnf1 |> select(input,target),
              by = c("gene_names" = "input")) |> unique()

ens_sig = bdnf1hr_v_control |> filter(adj_p_val < 0.1) |> 
    pull(target) |> unique()

total_sig_phos_go = clusterProfiler::enrichGO(ens_sig,
                                               keyType = 'ENSEMBL',
                                               OrgDb = org.Hs.eg.db,
                                               ont = 'ALL',
                                               readable = TRUE)

clusterProfiler::dotplot(total_sig_phos_go)
clusterProfiler::cnetplot(total_sig_phos_go)
