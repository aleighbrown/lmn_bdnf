
#data from https://doi.org/10.1073/pnas.1806660115
library(tidyverse)
library(data.table)
library(DESeq2)
library(clusterProfiler)
library(ggridges)
library(GeneOverlap)
# shrink lfc before input -------------------------------------------------


hour_one_featurecounts <- readRDS("data/hour_one_featurecounts.RDS")
hour_two_featurecounts <- readRDS("data/hour_two_featurecounts.RDS")
hour_six_featurecounts <- readRDS("data/hour_six_featurecounts.RDS")

shrunk_one <- lfcShrink((hour_one_featurecounts$deseq_obj),"comparison_BDNF_vs_control") |> 
    as.data.frame() |> 
    tibble::rownames_to_column('ensgene') |> 
    mutate(ensgene = gsub("\\..*", "", ensgene)) |> 
    left_join(annotables::grch38, by = c('ensgene')) |> as.data.table()
    
shrunk_two <- lfcShrink((hour_two_featurecounts$deseq_obj),"comparison_BDNF_vs_control") |> 
    as.data.frame() |> 
    tibble::rownames_to_column('ensgene') |> 
    mutate(ensgene = gsub("\\..*", "", ensgene)) |> 
    left_join(annotables::grch38, by = c('ensgene')) |> as.data.table()

shrunk_six <- lfcShrink((hour_six_featurecounts$deseq_obj),"comparison_BDNF_vs_control") |> 
    as.data.frame() |> 
    tibble::rownames_to_column('ensgene') |> 
    mutate(ensgene = gsub("\\..*", "", ensgene)) |> 
    left_join(annotables::grch38, by = c('ensgene')) |> as.data.table()

new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(total_rna_sig = case_when(padj < 0.1 & log2FoldChange > 0.75 ~ " Upregulated",
                                     padj < 0.1 & log2FoldChange < -0.75 ~ "Downregulated",
                                     T ~ "not_significant"))

new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(log2Fold_newRNA = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    mutate(gene = gsub("\\..*", "", gene))

# run de on the other BDNF paper ------------------------------------------

meta = fread("data/other_bdnf_paper/E-MTAB-6975.sdrf.txt")
meta = janitor::clean_names(meta)
meta = meta |> 
    select(source_name,extract_name,
           factor_value_sampling_time_point,factor_value_stimulus) |> 
    unique()

meta = meta |> 
    mutate(condition = case_when(factor_value_stimulus == 'brain-derived neurotrophic factor' ~ "BDNF",
                                 factor_value_stimulus == 'none' ~ "_Control",
                                 factor_value_stimulus == 'neurotrophin-4' ~ "NT4",
                                 TRUE ~ "TRKB_agonist")) 

raw_counts = fread("data/other_bdnf_paper/raw_counts.txt")


conv_counts = raw_counts |> 
    tibble::column_to_rownames("tracking_id")

conv_counts = conv_counts[meta$source_name]


# run differential expression on 2 hours -------------------------------------------

meta2bdnf = meta |> 
    filter(factor_value_sampling_time_point == 2) |> 
    filter(condition %in% c("BDNF","_Control")) |> 
    tibble::column_to_rownames('source_name')


raw_counts2bdnf = conv_counts[meta2bdnf$extract_name]
    

dds_bdnf2 = DESeqDataSetFromMatrix(raw_counts2bdnf,
                                        colData = meta2bdnf, 
                                        design = ~ condition)
#that just created the object
#to actually run the analysis we just call "DESeq"
dds_bdnf2 = DESeq(dds_bdnf2)

shrunk_other_two <- lfcShrink((dds_bdnf2),"condition_BDNF_vs__Control") |> 
    as.data.frame() |> 
    tibble::rownames_to_column('ensgene') |> 
    mutate(ensgene = gsub("\\..*", "", ensgene)) |> 
    left_join(annotables::grch38, by = c('ensgene')) |> as.data.table()

#we can quickly view the results with the function results()
hour2_res = results(dds_bdnf2) |> 
    as.data.frame() |> 
    tibble::rownames_to_column('ensgene') |> 
    left_join(annotables::grch38, by = c('ensgene')) |> as.data.table()
    

hour2_res |> 
    left_join(new_ratio_bayesian_p_de[time == 2],by = c('ensgene' = 'gene'),
              suffix = c("_other","_our")) |> 
    mutate(z_score_other = scale(log2FoldChange_other)) |> 
    mutate(z_score_ours = scale(log2FoldChange_our)) |> 
    filter(padj_other < 0.1 | padj_our < 0.1) |> 
    # filter(total_rna_sig_our == 'Downregulated') |> 
    # filter(padj_other < 0.1) |>
    # filter(padj_our < 0.1) |>
        # ggplot(aes(x = z_score_ours,
        #            y = z_score_other)) + 
        # 
    ggplot(aes(x = log2FoldChange_our,
               y = log2FoldChange_other)) +
    geom_point(alpha = 0.3) + 
    # geom_text_repel(aes(label = gene_name)) + 
    # geom_hline(linetype = 'dashed',yintercept = 0) +
    # geom_vline(linetype = 'dashed',xintercept = 0) +
    ggpubr::stat_cor(cor.coef.name = 'rho',method = 'spearman') +
    geom_smooth(method = 'lm') +
    ggpubr::theme_pubr() 




hour2_res |> 
    left_join(new_ratio_bayesian_p_de[time == 2],by = c('ensgene' = 'gene'),
              suffix = c("_other","_our")) |> 
    # filter(total_rna_sig == 'Downregulated') |> 
    # filter(new_rna_sig == 'bdnf_higher_new_rna') |> 
    filter(padj_other < 0.1 & padj_our < 0.1) |> 
    ggplot(aes(x = log2FoldChange_our,
               y = log2FoldChange_other)) + 
    geom_point(alpha = 0.1) + 
    # geom_text_repel(aes(label = gene_name)) + 
    geom_hline(linetype = 'dashed',yintercept = 0) +
    geom_vline(linetype = 'dashed',xintercept = 0) +
    ggpubr::stat_cor(cor.coef.name = 'rho',method = 'spearman') +
    geom_smooth(method = 'lm') +
    ggpubr::theme_pubr()

shrunk_other_two |> 
    left_join(shrunk_two,by = c('ensgene'),
              suffix = c("_other","_our")) |> 
    # filter(total_rna_sig == 'Downregulated') |> 
    # filter(new_rna_sig == 'bdnf_higher_new_rna') |> 
    filter(padj_other < 0.1 & padj_our < 0.1) |> 
    ggplot(aes(x = log2FoldChange_our,
               y = log2FoldChange_other)) + 
    geom_point(alpha = 0.1) + 
    # geom_text_repel(aes(label = gene_name)) + 
    geom_hline(linetype = 'dashed',yintercept = 0) +
    geom_vline(linetype = 'dashed',xintercept = 0) +
    ggpubr::stat_cor() +
    geom_smooth(method = 'lm') +
    ggpubr::theme_pubr()

# hypergeometic test -------------------------------------------
# union(hour2_res[baseMean > 5]$ensgene,uni) |> length() - n genes expressed in either

ours = new_ratio_bayesian_p_de[log2FoldChange < -0.75 & padj < 0.1 & time == 2,gene]
theirs = hour2_res[log2FoldChange < -0.75 & padj < 0.1,ensgene]

go.obj <- newGeneOverlap(ours,
                         theirs,
                         31039)

go.obj <- testGeneOverlap(go.obj)


# hour2_res |> 
#     left_join(new_ratio_bayesian_p_de[time == 2],by = c('ensgene' = 'gene'),
#               suffix = c("_other","_our")) |> 
#     filter(total_rna_sig_our == 'Downregulated') |> 
#     filter(new_rna_sig == 'bdnf_higher_new_rna')  
#     janitor::tabyl(total_rna_sig_other) |> 
#     dplyr::rename(n_gene = n) |> 
#     ggplot(aes(x = total_rna_sig_other,
#                y = n_gene)) + 
#     geom_col() 
#     coord_flip() +
#     ggpubr::theme_pubr() +
#     facet_wrap(~(L1),nrow = 3,scales = 'free') 
    


hour2_res = hour2_res |> 
    mutate(total_rna_sig = case_when(padj < 0.05 & log2FoldChange > 1 ~ "upregulated",
                                     padj < 0.05 & log2FoldChange < -1 ~ "downregulated",
                                     T ~ "not_significant")) 

other_uni = hour2_res[baseMean > 10,ensgene]

formula_res1_other <- compareCluster(ensgene~total_rna_sig, 
                               data=hour2_res[total_rna_sig != 'not_significant'], 
                               fun="enrichGO",
                               OrgDb='org.Hs.eg.db',
                               keyType = "ENSEMBL",
                               universe = other_uni)

dotplot(formula_res1_other, x="total_rna_sig") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab(element_blank()) +
    theme(axis.text.y  = element_text(size = 9)) 
    
    

    
    
    
dds_all = DESeqDataSetFromMatrix(conv_counts,
                                   colData = meta, 
                                   design = ~ condition)

dds_all = DESeq(dds_all)

# pcaExplorer::pcaExplorer(dds_all)

# chandran magenta genes --------------------------------------------------


#https://doi.org/10.1016/j.neuron.2016.01.034
all_genes = readxl::read_excel('data/chandra2016_modules.xlsx')
mapping = gprofiler2::gorth(all_genes$Symbol,source_organism = 'rnorvegicus','hsapiens')
all_genes = all_genes |> 
    left_join(mapping, by = c('Symbol'  = "input"))

magenta_genes = all_genes |> filter(Module == 'magenta')

# Jacobi single cell ------------------------------------------------------


up_C_P_ko = c("ENSG00000164326", "ENSG00000185129", "ENSG00000068383", "ENSG00000007402", 
  "ENSG00000113657", "ENSG00000069869", "ENSG00000168028", "ENSG00000288920", 
  "ENSG00000101210", "ENSG00000052802", "ENSG00000091664", "ENSG00000167522",
  "ENSG00000168036", "ENSG00000171246", "ENSG00000159363", "ENSG00000166165",
  "ENSG00000176887", "ENSG00000036257", "ENSG00000196136", "ENSG00000134369",
  "ENSG00000143862", "ENSG00000182253",  "ENSG00000170345", "ENSG00000166068", 
  "ENSG00000154478", "ENSG00000163399", "ENSG00000090863", "ENSG00000106772", 
  "ENSG00000171858", "ENSG00000206418", "ENSG00000176697", "ENSG00000172458", 
  "ENSG00000080371", "ENSG00000177519", "ENSG00000198727",  "ENSG00000188229", 
  "ENSG00000054793", "ENSG00000134046", "ENSG00000132824", "ENSG00000171208", 
  "ENSG00000069482", "ENSG00000125740", "ENSG00000168610", "ENSG00000130558", 
  "ENSG00000188483", "ENSG00000177733", "ENSG00000146676", "ENSG00000176396", 
  "ENSG00000198886", "ENSG00000166710", "ENSG00000105409", "ENSG00000196632", 
  "ENSG00000112972", "ENSG00000131773", "ENSG00000138083", "ENSG00000198840", 
  "ENSG00000111674", "ENSG00000130725", "ENSG00000102003", "ENSG00000258947", 
  "ENSG00000156642", "ENSG00000236279", "ENSG00000177606", "ENSG00000176788")

up_C_P_ko = data.table(ensgene = up_C_P_ko) |> 
    mutate(gs_name = "up_C_P_ko")

down_C_PS_ko = c("ENSG00000117632", "ENSG00000161970", "ENSG00000236782", 
               "ENSG00000145708", "ENSG00000147434", "ENSG00000115694",
               "ENSG00000126351")

down_C_PS_ko = data.table(ensgene = down_C_PS_ko) |> 
    mutate(gs_name = "down_C_PS_ko")

up_C_PS_ko = c("ENSG00000277586", "ENSG00000134369", "ENSG00000139289", 
               "ENSG00000176887", "ENSG00000091844", "ENSG00000155367",  
               "ENSG00000070961", "ENSG00000171246", "ENSG00000167552", 
               "ENSG00000130303", "ENSG00000141404", "ENSG00000119917", 
               "ENSG00000187608", "ENSG00000123095", "ENSG00000188517",
               "ENSG00000214063", "ENSG00000162852", "ENSG00000168610", 
               "ENSG00000145675", "ENSG00000163399", "ENSG00000175130", 
               "ENSG00000130429", "ENSG00000109501", "ENSG00000100362", 
               "ENSG00000130222", "ENSG00000177606", "ENSG00000160932", 
               "ENSG00000128245", "ENSG00000176697", "ENSG00000272398", 
               "ENSG00000164326",  "ENSG00000258947", "ENSG00000069482", 
               "ENSG00000204264", "ENSG00000177519", "ENSG00000204010", 
               "ENSG00000006128", "ENSG00000196136", "ENSG00000166710")

up_C_PS_ko = data.table(ensgene = up_C_PS_ko) |> 
    mutate(gs_name = "up_C_PS_ko")

the_genes_in306 = c("ENSG00000179820", "ENSG00000135631", "ENSG00000166963", 
                    "ENSG00000060762",  "ENSG00000141759",  "ENSG00000139641", 
                    "ENSG00000171488", "ENSG00000177519",  "ENSG00000069956", 
                    "ENSG00000177606", "ENSG00000165030", "ENSG00000175416", 
                    "ENSG00000015592", "ENSG00000168092", "ENSG00000183943", 
                    "ENSG00000242498", "ENSG00000095794", "ENSG00000142669", 
                    "ENSG00000117592", "ENSG00000165868", "ENSG00000087074", 
                    "ENSG00000137076", "ENSG00000144339", "ENSG00000100243", 
                    "ENSG00000165566", "ENSG00000109472", "ENSG00000068971",  
                    "ENSG00000168461", "ENSG00000102900", "ENSG00000115226", 
                    "ENSG00000103363", "ENSG00000141424", "ENSG00000178209", 
                    "ENSG00000197043", "ENSG00000164111", "ENSG00000213465", 
                    "ENSG00000196136", "ENSG00000129235",  "ENSG00000184076",
                    "ENSG00000162909", "ENSG00000105223",    "ENSG00000164733",
                    "ENSG00000154620", "ENSG00000131016",  "ENSG00000159713", 
                    "ENSG00000183150",  "ENSG00000013364", "ENSG00000170876", 
                    "ENSG00000117632", "ENSG00000173267", "ENSG00000153113",  
                    "ENSG00000100433", "ENSG00000156103", "ENSG00000144959", 
                    "ENSG00000162852", "ENSG00000070614", "ENSG00000184979", 
                    "ENSG00000198910", "ENSG00000120129", "ENSG00000110218", 
                    "ENSG00000248485",  "ENSG00000167535", "ENSG00000060982", 
                    "ENSG00000242173", "ENSG00000173801", "ENSG00000105835", 
                    "ENSG00000148120", "ENSG00000143179", "ENSG00000196136", 
                    "ENSG00000188523", "ENSG00000164054", "ENSG00000136274", 
                    "ENSG00000111639",  "ENSG00000106397", "ENSG00000143198",
                    "ENSG00000136153", "ENSG00000079616", "ENSG00000233276",
                    "ENSG00000196136", "ENSG00000110002", "ENSG00000112245", 
                    "ENSG00000112245", "ENSG00000162636",  "ENSG00000122694", 
                    "ENSG00000124535",  "ENSG00000148908", "ENSG00000112697", 
                    "ENSG00000139645", "ENSG00000103966", "ENSG00000160789", 
                    "ENSG00000048140", "ENSG00000075975",  "ENSG00000076706", 
                    "ENSG00000178980", "ENSG00000075223", "ENSG00000213977", 
                    "ENSG00000138741", "ENSG00000159228", "ENSG00000151572", 
                    "ENSG00000170915", "ENSG00000198948", "ENSG00000073792",  
                    "ENSG00000125995", "ENSG00000030582",   "ENSG00000196924",
                    "ENSG00000152503", "ENSG00000167972", "ENSG00000178567", 
                    "ENSG00000124570", "ENSG00000170502", "ENSG00000213719", 
                    "ENSG00000187608", "ENSG00000167460", "ENSG00000271303", 
                    "ENSG00000143549", "ENSG00000127314", "ENSG00000119242", 
                    "ENSG00000221946",  "ENSG00000173221", "ENSG00000165795", 
                    "ENSG00000175040", "ENSG00000136783", "ENSG00000165028", 
                    "ENSG00000136895", "ENSG00000184937", "ENSG00000113312",
                    "ENSG00000091436", "ENSG00000104435", "ENSG00000104219",
                    "ENSG00000160932",  "ENSG00000196136",  "ENSG00000163320", 
                    "ENSG00000227051",  "ENSG00000087157", "ENSG00000077238", 
                    "ENSG00000166920", "ENSG00000049130", "ENSG00000159200",
                    "ENSG00000135404", "ENSG00000197747", "ENSG00000175130", 
                    "ENSG00000213281", "ENSG00000173786",  "ENSG00000122012", 
                    "ENSG00000115267", "ENSG00000105939", "ENSG00000166342", 
                    "ENSG00000091664", "ENSG00000113594", "ENSG00000164687", 
                    "ENSG00000155366", "ENSG00000115758", "ENSG00000151929", 
                    "ENSG00000175535", "ENSG00000169750", "ENSG00000135636", 
                    "ENSG00000134352",   "ENSG00000089127", "ENSG00000172469", 
                    "ENSG00000075426", "ENSG00000129460", "ENSG00000120708", 
                    "ENSG00000118985", "ENSG00000130222", "ENSG00000096433", 
                    "ENSG00000131236", "ENSG00000157601", "ENSG00000102007", 
                    "ENSG00000102007", "ENSG00000203485", "ENSG00000159348", 
                    "ENSG00000137285", "ENSG00000107249", "ENSG00000135821",
                    "ENSG00000184144", "ENSG00000169627", "ENSG00000183336", 
                    "ENSG00000187091",   "ENSG00000131724",  "ENSG00000132256", 
                    "ENSG00000157379", "ENSG00000151651", "ENSG00000162704", 
                    "ENSG00000136238", "ENSG00000178531", "ENSG00000132196", 
                    "ENSG00000028137", "ENSG00000240065", "ENSG00000182718", 
                    "ENSG00000137959", "ENSG00000237515", "ENSG00000100097",  
                    "ENSG00000132530", "ENSG00000115884", "ENSG00000023191",  
                    "ENSG00000034510", "ENSG00000103522", "ENSG00000171621", 
                    "ENSG00000134207", "ENSG00000147065", "ENSG00000107201", 
                    "ENSG00000092820", "ENSG00000116661", "ENSG00000204264", 
                    "ENSG00000153048", "ENSG00000130429", "ENSG00000034510", 
                    "ENSG00000189159", "ENSG00000137965", "ENSG00000204010", 
                    "ENSG00000159176", "ENSG00000162772", "ENSG00000180354",
                    "ENSG00000133874", "ENSG00000141934", "ENSG00000165233", 
                    "ENSG00000130303", "ENSG00000185338", "ENSG00000134242", 
                    "ENSG00000171551", "ENSG00000134072", "ENSG00000182645",  
                    "ENSG00000106211", "ENSG00000110148", "ENSG00000126561", 
                    "ENSG00000221968", "ENSG00000105514", "ENSG00000164099", 
                    "ENSG00000131981", "ENSG00000067082", "ENSG00000162493", 
                    "ENSG00000177409", "ENSG00000169504", "ENSG00000163563", "ENSG00000291309", "ENSG00000119922", "ENSG00000184371", "ENSG00000026025", "ENSG00000196136",  "ENSG00000176014", "ENSG00000205116", "ENSG00000262406", "ENSG00000069482",  "ENSG00000183665", "ENSG00000168528", "ENSG00000124762", "ENSG00000026508", "ENSG00000196843", "ENSG00000172020",  "ENSG00000176170",  "ENSG00000101000", "ENSG00000100336", "ENSG00000100342", "ENSG00000128284", "ENSG00000128335", "ENSG00000092621", "ENSG00000172183", "ENSG00000198576",  "ENSG00000163191", "ENSG00000134321",  "ENSG00000088836", "ENSG00000135070", "ENSG00000291309", "ENSG00000158710",       "ENSG00000121797")


the_genes_in306 = data.table(ensgene = the_genes_in306) |> 
    mutate(gs_name = "the_genes_in306")


magenta_genes = data.table(ensgene = magenta_genes |> filter(!is.na(ortholog_ensg)) |> pull(ortholog_ensg)) |> 
    mutate(gs_name = "magenta_genes")


# DOI:https://doi.org/10.1016/j.neuron.2020.07.026
common_up_injury = readLines('data/common_up_injury.txt')
common_up_injury = data.table(ensgene = common_up_injury) |> 
    mutate(gs_name = "common_up_injury")

# gsea on injury genes ----------------------------------------------------


baseline_manipulations = rbind(up_C_P_ko,up_C_PS_ko,down_C_PS_ko,the_genes_in306,magenta_genes,common_up_injury) |> 
    dplyr::relocate(gs_name) |> 
    unique()

gsea1_total_shrunk = shrunk_one |> 
    # filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) |>
    filter(!is.na(log2FoldChange)) |>
    select(ensgene,log2FoldChange) |> 
    arrange(-log2FoldChange) %>% unique()

gsea_vec1_total_shrunk = gsea1_total_shrunk$log2FoldChange
names(gsea_vec1_total_shrunk) = gsea1_total_shrunk$ensgene


gsea2_total_shrunk = shrunk_two |> 
    # filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) |>
    filter(!is.na(log2FoldChange)) |>
    select(ensgene,log2FoldChange) |> 
    arrange(-log2FoldChange) %>% unique()

gsea_vec2_total_shrunk = gsea2_total_shrunk$log2FoldChange
names(gsea_vec2_total_shrunk) = gsea2_total_shrunk$ensgene

gsea6_total_shrunk = shrunk_six |> 
    # filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) |>
    filter(!is.na(log2FoldChange)) |>
    select(ensgene,log2FoldChange) |> 
    arrange(-log2FoldChange) %>% unique()

gsea_vec6_total_shrunk = gsea6_total_shrunk$log2FoldChange
names(gsea_vec6_total_shrunk) = gsea6_total_shrunk$ensgene


gsea1_new_cate = new_ratio_bayesian_p_de |> 
    filter(time == 1) |> 
    # filter(new_rna_sig == 'bdnf_higher_new_rna') |> 
    # filter(total_rna_sig == 'Downregulated') |> 
    filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) |>
    select(gene,log2Fold_newRNA) |> 
    arrange(-log2Fold_newRNA) %>% 
    unique() %>% 
    add_count(gene) %>% 
    filter(n == 1)

gsea_vec1_new = gsea1_new_cate$log2Fold_newRNA
names(gsea_vec1_new) = gsea1_new_cate$gene

gsea2_new_cate = new_ratio_bayesian_p_de |> 
    filter(time == 2) |> 
    # filter(new_rna_sig == 'bdnf_higher_new_rna') |> 
    # filter(total_rna_sig == 'Downregulated') |> 
    filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) |>
    select(gene,log2FoldChange,log2Fold_newRNA) |> 
    arrange(-log2Fold_newRNA) %>% unique() %>% 
    add_count(gene) %>% 
    filter(n == 1)

gsea_vec2_new = gsea2_new_cate$log2Fold_newRNA
names(gsea_vec2_new) = gsea2_new_cate$gene

gsea6_new_cate = new_ratio_bayesian_p_de |> 
    filter(time == 6) |> 
    # filter(new_rna_sig == 'bdnf_higher_new_rna') |> 
    # filter(total_rna_sig == 'Downregulated') |> 
    filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) |>
    select(gene,log2FoldChange,log2Fold_newRNA) |> 
    arrange(-log2Fold_newRNA) %>% unique() %>% 
    add_count(gene) %>% 
    filter(n == 1)

gsea_vec6_new = gsea6_new_cate$log2Fold_newRNA
names(gsea_vec6_new) = gsea6_new_cate$gene



em_shrunk_1 <- clusterProfiler::GSEA(gsea_vec1_total_shrunk, TERM2GENE = baseline_manipulations,pvalueCutoff = 0.9)
em_shrunk_2 <- clusterProfiler::GSEA(gsea_vec2_total_shrunk, TERM2GENE = baseline_manipulations,pvalueCutoff = 0.9)
em_shrunk_6 <- clusterProfiler::GSEA(gsea_vec6_total_shrunk, TERM2GENE = baseline_manipulations,pvalueCutoff = 0.9)

em1_new <- clusterProfiler::GSEA(gsea_vec1_new, TERM2GENE = baseline_manipulations,pvalueCutoff = 0.9)
em2_new <- clusterProfiler::GSEA(gsea_vec2_new, TERM2GENE = baseline_manipulations,pvalueCutoff = 0.9)
em6_new <- clusterProfiler::GSEA(gsea_vec6_new, TERM2GENE = baseline_manipulations,pvalueCutoff = 0.9)


clusterProfiler::ridgeplot(em_shrunk_1) + ggtitle('total 1 hr')
clusterProfiler::ridgeplot(em_shrunk_2) + ggtitle('total 2 hr')
clusterProfiler::ridgeplot(em_shrunk_6) + ggtitle('total 6 hr')



clusterProfiler::ridgeplot(em1_new) + ggtitle('new 1 hr')
clusterProfiler::ridgeplot(em2_new) + ggtitle('new 2 hr')
clusterProfiler::ridgeplot(em6_new) + ggtitle('new 6 hr')



new_rna_enrichment = rbind((em1_new@result |> as.data.table() |> 
    mutate(time = 1) |> 
    mutate(type = 'new_RNA')),
(em2_new@result |> as.data.table() |> 
    mutate(time = 2) |> 
    mutate(type = 'new_RNA')),
(em6_new@result |> as.data.table() |> 
    mutate(time = 6) |> 
    mutate(type = 'new_RNA')))

total_rna_enrichment = rbind((em_shrunk_1@result |> as.data.table() |> 
                                mutate(time = 1) |> 
                                mutate(type = 'total_RNA')),
                           (em_shrunk_2@result |> as.data.table() |> 
                                mutate(time = 2) |> 
                                mutate(type = 'total_RNA')),
                           (em_shrunk_6@result |> as.data.table() |> 
                                mutate(time = 6) |> 
                                mutate(type = 'total_RNA')))


full_enrichment = rbind(total_rna_enrichment,new_rna_enrichment)
shrunken_together = (rbind((shrunk_one |> 
                                mutate(time = 1)),
                           (shrunk_two |> 
                                mutate(time = 2)),
                           (shrunk_six |> 
                                mutate(time = 6)))) |> 
    unique() |> 
    select(ensgene,log2FoldChange,time) |> 
    # filter(ensgene %in% baseline_manipulations$ensgene) |> 
    dplyr::rename(gene = ensgene,
                  value = log2FoldChange) |> 
    mutate(variable = 'total_RNA')
    


new_ratio_bayesian_p_de |> 
    filter(n_samp_passing_bdnf >=2 & n_samp_passing_control >= 2) |>
    select(gene,time,log2Fold_newRNA) |>
    filter(gene %in% baseline_manipulations$ensgene) |> 
    dplyr::rename(value = log2Fold_newRNA) |> 
    mutate(variable = 'new_RNA') |> 
    rbind(shrunken_together) |> 
    unique() |> 
    left_join(baseline_manipulations,by = c('gene' = 'ensgene')) |> 
    unique() |> 
    filter(!is.na(gs_name))  |> 
    left_join(full_enrichment,by = c("gs_name" = 'ID','variable' = 'type','time')) |> 
    # filter(!is.na(setSize)) |> 
    filter(!grepl("^up_",gs_name)) |> 
    filter(!grepl("^down_",gs_name)) |> 
    # filter(time != 1) |>
    # filter(time == 1 & variable == 'total_RNA') |> 
    mutate(sig_pvalue = ifelse(p.adjust < 0.05,p.adjust,NA)) |> 
    ggplot(aes(x = value,
               y = gs_name,
               fill = sig_pvalue)) + 
    geom_density_ridges() + 
    ggpubr::theme_pubr() +
    geom_vline(xintercept = 0,linetype = 'dashed') + 
    facet_wrap(~time + variable,scales = 'free',nrow = 3) +
    scale_fill_viridis_c()



baseline_manipulations |> 
    filter(gs_name %in% c("common_up_injury",
                          "magenta_genes",
                          "the_genes_in306")) |> 
    mutate(present = TRUE) |> 
    unique() |> 
    pivot_wider(names_from = "gs_name",
                values_from = 'present',values_fill = FALSE) |> 
    tibble::column_to_rownames('ensgene') |> 
    eulerr::euler() |> 
    plot(,quantities = TRUE,fill = "transparent")

    

baseline_manipulations |> 
    filter(gs_name %in% c("common_up_injury",
                          "magenta_genes",
                          "the_genes_in306")) |> 
    mutate(present = TRUE) |> 
    unique() |> 
    pivot_wider(names_from = "gs_name",
                values_from = 'present',values_fill = FALSE) |> 
    left_join(new_ratio_bayesian_p_de,by = c("ensgene" = "gene")) |> 
    filter(time == 2) |> 
    filter(common_up_injury == TRUE & 
               magenta_genes == TRUE &
               the_genes_in306 == TRUE) |> 
    select(ensgene,gene_name,new_rna_sig,log2FoldChange,mean_bdnf_ntr,mean_control_ntr,bayesian_p) |> 
    View()


lil = full_enrichment |>   filter(ID %in% c("common_up_injury",
                                                  "magenta_genes",
                                                  "the_genes_in306")) |>  
    filter(type == 'total_RNA') |> 
    filter(time == 1) |> 
    select(ID,core_enrichment) |> 
    separate_rows(core_enrichment)

lil |> 
    mutate(present = TRUE) |> 
    unique() |> 
    pivot_wider(names_from = "ID",
                values_from = 'present',values_fill = FALSE) |> 
    tibble::column_to_rownames('core_enrichment') |> 
    eulerr::euler() |> 
    plot(,quantities = TRUE,fill = "transparent")
