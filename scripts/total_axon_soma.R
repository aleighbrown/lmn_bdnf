library(clusterProfiler)
library(enrichplot)
library(limma)
library(ggrepel)
library(rlang)
library(ComplexHeatmap)
source("~/Documents/GitHub/bdnf_4su/run_standard_limma.R")
library(dplyr)
library(gprofiler2)
library(rlang)
library(org.Hs.eg.db)
library(tidyverse)
# Van Nostrand, E.L., Freese, P., Pratt, G.A. et al. A large-scale binding and functional map of human RNA-binding proteins. Nature 583, 711–719 (2020). https://doi.org/10.1038/s41586-020-2077-3

rbps = readxl::read_excel('/Users/annaleigh/Documents/GitHub/bdnf_4su/data/41586_2020_2077_MOESM3_ESM.xlsx') |> 
    janitor::clean_names() |> 
    as.data.table()


# https://esbl.nhlbi.nih.gov/Databases/KSBP2/Targets/Lists/E3-ligases/
e3_ligases = readxl::read_excel('/Users/annaleigh/Documents/GitHub/bdnf_4su/data/Definite Ligase List.xlsx') |> janitor::clean_names() |> as.data.table()


# 10.7554/eLife.64943
kinases = readxl::read_excel('/Users/annaleigh/Documents/GitHub/bdnf_4su/data/elife-64943-supp1.xlsx') |> 
    janitor::clean_names() |> 
    as.data.table()

#https://depod.bioss.uni-freiburg.de/download.php
phophotases = readxl::read_excel('/Users/annaleigh/Documents/GitHub/bdnf_4su/data/PPases_in_DEPOD_201906.xls') |> 
    janitor::clean_names() |> 
    as.data.table()



amigoingmad = function(package = "dplyr", fix = TRUE, iteration = 0 ) {
    if (iteration > 1) {
        stop("Can't fix.")
    }
    conf = unique(conflicts())
    want_package = paste0("package:", package)
    conflicts_desired_package = conf[conf %in% ls(want_package)]
    conflict_envs = sapply(conflicts_desired_package, FUN = function(x) {
        environmentName(pryr::where(x, globalenv())) # relative to globalenv, not pkg env
    })
    is_good = conflict_envs == want_package
    potentially_bad_confs = conflicts_desired_package[!is_good]
    potentially_bad_envs = conflict_envs[!is_good]
    have_to_fix = rep(FALSE, length(potentially_bad_confs))
    for (i in seq_along(potentially_bad_confs)) {
        if (!identical(body(get(potentially_bad_confs[i], pos = want_package)),
                       body(get(potentially_bad_confs[i])))) {
            have_to_fix[i] = TRUE
        }
    }
    if (any(have_to_fix)) {
        message("The following functions don't have the environment you want.")
        print(data.frame(`function.` = potentially_bad_confs[have_to_fix],
                         environment = potentially_bad_envs[have_to_fix]),
              row.names = F)
        if (fix) {
            base::detach(name = want_package, character.only = TRUE)
            base::library(package, character.only = TRUE)
            message("Tried to fix this, calling myself again to make sure...")
            amigoingmad(package, fix, iteration + 1)
            message("Sanity restored!")
        }
    } else if (iteration == 0) {
        message("Everything looks normal. Maybe it's you.")
    }
}

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
amigoingmad()
e3_ligases = add_the_annotation(e3_ligases,'swiss_prot')
rbps = add_the_annotation(rbps,'x1')
kinases = add_the_annotation(kinases,'official_symbol')
phophotases = add_the_annotation(phophotases,'uni_prot_ac')



new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv") 

new_ratio_bayesian_p_de = new_ratio_bayesian_p_de %>% 
    dplyr::select(-gene_name) %>% 
    mutate(gene = gsub("\\..*", "", gene)) %>%  
    left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol), by = c("gene" = "ensgene")) %>% 
    unique()

bg_expressed_in_LMN = new_ratio_bayesian_p_de[baseMean > 15,gene]

metadata = readxl::read_excel("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/samplesheet.xlsx")
metadata = metadata |> 
    mutate(new_column = paste(condition,sample,sep = "_"))

total_proteins = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/report.pg_matrix.tsv") |> 
    janitor::clean_names() |> 
    as.data.table()

colnames(total_proteins) = gsub("c_beth_raw_drix067_apriil_tp_20240417_bg_drix067_tp_",'',colnames(total_proteins))
colnames(total_proteins) = gsub("(\\d)_.*", "\\1", colnames(total_proteins))


recode_vec <- setNames(metadata$sample, metadata$new_column)

raw_input <- total_proteins %>% 
    dplyr::rename(any_of(recode_vec))
raw_input = raw_input |> tibble::rowid_to_column()

# set up my crap database -------------------------------------------------
skin_enriched = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/skin_enriched.tsv")
skin_enhanced = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/skin_enhanced.tsv")

mapped_genes = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/gProfiler_hsapiens_05-04-2024_13-08-42.csv")

histone_uniprots = mapped_genes |>
    filter(grepl('histone',description)) |> pull(initial_alias)

beth_crap = fread('/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/beth_crap.csv')
contaims_more = c("K1C9", "K1C10", "K2C1",
                  "K2C2", "K2C5", "K2C6", "K2C7", "K2C8",
                  "K2C9", "K2CE","ALBU",
                  "K1C14")

histones = c("H2B","H2A1B","H37","NPM","H2A2","H3-3A", "H3-3B", 
             "H3-4", "H3-5", "H3C1", "H3C10", "H3C11", "H3C12", "H3C13",
             "H3C14", "H3C15", "H3C2", "H3C3", "H3C4", "H3C6", "H3C7", "H3C8",
             "H4","H2AX","H12","H1-5")

#get rid of proteins that are 'skin enriched' or 'skin enhanced'
skin_enriched_rows = raw_input |> separate_rows(protein_ids,sep = ";") |> filter(protein_ids %in% skin_enriched$Uniprot) |> pull(rowid) |> unique()
skin_enhanced_rows = raw_input |> separate_rows(protein_ids,sep = ";") |> filter(protein_ids %in% skin_enhanced$Uniprot) |> pull(rowid) |> unique()
#get rid of some more keratines
more_skin_rows = raw_input |> separate_rows(protein_names,sep = ";") |>   filter(grepl(paste0(contaims_more,collapse = '|'),protein_names)) |> pull(rowid) |> unique()
#get rid of proteins that are in Beth's crapome
beth_rows = raw_input |> separate_rows(protein_ids,sep = ";") |> filter(protein_ids %in% beth_crap$uniprot_id) |> pull(rowid) |> unique()
# #get rid of proteins that are histones
# histone_rows = raw_input |> separate_rows(protein_ids,sep = ";") |> filter(protein_ids %in% histone_uniprots) |> pull(rowid) |> unique()
# more_histone_rows = raw_input |> separate_rows(protein_names,sep = ";") |>   filter(grepl(paste0(histones,collapse = '|'),protein_names)) |> pull(rowid) |> unique()
#And now getting rid of all that garbage
rawest_input = raw_input
raw_input = raw_input |> 
    # filter(!(rowid %in% unique(c(skin_enhanced_rows,skin_enhanced_rows,more_skin_rows,beth_rows,histone_rows,more_histone_rows)))) |> 
    filter(!(rowid %in% unique(c(skin_enriched_rows,skin_enhanced_rows,more_skin_rows,beth_rows))))

raw_input %>%
    select(protein_ids, protein_names, 
           matches("control_soma"),matches('control_axon'),
           matches("^BDNF_soma"),matches('^BDNF_axon')) %>% 
    mutate(NA_counts = rowSums(select(., control_soma_a1:BDNF_axon_b5) %>% is.na(), na.rm = TRUE)) %>%
    group_by(protein_ids) |>
    slice(which.min(NA_counts))  |>
    ungroup() |>
    mutate(new_column = paste(protein_names,protein_ids,sep = "|")) |>
    tibble::column_to_rownames('new_column') |>
    select(control_soma_a1:BDNF_axon_b5) %>%
    # select(contains("axon")) %>%
    mutate(NA_counts = rowSums(is.na(.), na.rm = TRUE)) |>
    filter(NA_counts < 2) |> 
    select(-NA_counts) |>
    pheatmap::pheatmap(,show_rownames = FALSE,cluster_rows = FALSE,scale = 'row')



# Identifying axon samples with nuclear pore or polymerase ---------




for_clarification_tble = rawest_input |> 
    filter(genes %in% c("H1-0", "H1-10", "H1-2", "H1-2;H1-3;H1-4", "H1-2;H1-3;H1-4;H1-5", 
      "H1-2;H1-3;H1-4;H1-5;H1-6", "H1-2;H1-3;H1-4;H1-6", "H1-3;H1-4;H1-5", 
      "H1-4", "H1-5", "H1-6", "MAP2", "MAPT", "NFASC","MAP1B", "ATL1",
      "NUP107", "NUP133", "NUP153", "NUP155", "NUP160", "NUP205", "NUP214", 
      "NUP50", "NUP85", "NUP88", "NUP93", "NUP98", "POLA1", 
      "POLB", "POLD1", "POLD2", "POLD3", "POLG", "POLG2", "POLL", "TUBA1C"
    )) |> 
    select(genes,"control_soma_a1", "control_soma_a2", "control_soma_a3", "control_soma_a4", 
           "control_soma_a5", "BDNF_soma_c1", "BDNF_soma_c2", "BDNF_soma_c3", "iERKBDNF_soma_e1", 
           "iERKBDNF_soma_e2", "iERKBDNF_soma_e3", "iERKBDNF_soma_e4", "iERKBDNF_soma_e5", 
           "control_axon_d2", "control_axon_d3", "control_axon_d4", 
           "control_axon_d5", "BDNF_axon_b2", "BDNF_axon_b3", "BDNF_axon_b5", 
           "iERKBDNF_axon_f1", "iERKBDNF_axon_f2", "iERKBDNF_axon_f3", "iERKBDNF_axon_f4", 
           "iERKBDNF_axon_f5") |>
    rowwise() %>%
    mutate(observed_count = sum(!is.na(c_across(where(is.numeric))))) %>%
    ungroup() |> 
    filter(observed_count >=4) |> 
    select(-observed_count)  |> 
    filter(!is.na(genes)) |> 
    arrange(genes) |> 
    tibble::column_to_rownames('genes') |> 
    mutate_all(log2) 

my_sample_col = tibble(samples = colnames(for_clarification_tble)) |> 
    separate(samples,into = c("condition","compartment","replicate"),remove = FALSE) |> 
    tibble::column_to_rownames('samples') |> 
    select(-replicate)

my_colour = list(
    condition = c(control = "#F93FFF", BDNF = "#1B7DFD", iERKBDNF = "#EA8C55"),
    compartment = c(axon = "black", soma = "grey")
)
new_order = c("MAPT","MAP1B",'ATL1', "NFASC","TUBA1C","MAP2", "H1-0", "H1-10", "H1-2", "H1-2;H1-3;H1-4", "H1-2;H1-3;H1-4;H1-5", 
              "H1-2;H1-3;H1-4;H1-5;H1-6", "H1-2;H1-3;H1-4;H1-6", "H1-3;H1-4;H1-5", 
              "H1-4", "H1-5", "H1-6",  "NUP107", "NUP133", 
              "NUP153", "NUP155", "NUP160", "NUP205", "NUP214", "NUP50", "NUP85", 
              "NUP88", "NUP93", "NUP98", "POLA1", "POLB", "POLD1", "POLD2", 
              "POLD3", "POLG", "POLG2", "POLL")

for_clarification_tble <- for_clarification_tble[new_order, ]

pheatmap::pheatmap(for_clarification_tble,
                   annotation_colors = my_colour,
                   cluster_cols = FALSE, 
                   cluster_rows = FALSE,
                   annotation_col = my_sample_col,
                   scale = 'row',
                   show_colnames = TRUE)


# Proteome of the axon - Make PCA on AXON only to identify outliers -----------------------------------------
# axon_pca_input = raw_input %>%
#     select(protein_ids, protein_names, matches("_axon")) %>%
#     mutate(NA_counts = rowSums(select(., BDNF_axon_b1:control_axon_d5) %>% is.na(), na.rm = TRUE)) %>%
#     group_by(protein_ids) %>%
#     dplyr::slice(which.min(NA_counts))  |>
#     filter(NA_counts <= 4) |>
#     # select(-matches('b3')) |> #high amount of Nuclear Pore Proteins for axon
#     select(-NA_counts) |>
#     ungroup() |>
#     dplyr::relocate(protein_ids,protein_names) %>%
#     mutate(across(where(is.numeric), ~ ./median(.,na.rm = TRUE))) %>%
#     # mutate(across(where(is.numeric), ~ .* both_mean)) %>%
#     mutate(across(where(is.numeric), ~ log2(.))) %>%
#     mutate(new_column = paste(protein_names,protein_ids,sep = "|")) |>
#     tibble::column_to_rownames('new_column') |>
#     select(-protein_ids,-protein_names) %>%
#     na.omit()
# 
# axon_pca_input = raw_input %>%
#     dplyr::select(protein_ids, protein_names, matches("axon")) %>%
#     mutate(NA_counts = rowSums(dplyr::select(., control_axon_b1:iERKBDNF_axon_f5) %>% is.na(), na.rm = TRUE)) %>%
#     group_by(protein_ids) %>%
#     dplyr::slice(which.min(NA_counts))  |>
#     ungroup() |>
#     mutate(new_column = paste(protein_names,protein_ids,sep = "|")) |>
#     tibble::column_to_rownames('new_column') |>
#     filter(NA_counts == 0) |> 
#     # dplyr::select(-matches("d5|f2|b5|b2")) %>%
#     dplyr::select(control_axon_b1:iERKBDNF_axon_f5) %>% 
#     dplyr::select(-matches("iERK")) 
# 
axon_pca = axon_pca_input |>
    prcomp()
# 
axon_pca$rotation |>
    as.data.frame() |>
    tibble::rownames_to_column('sample') |>
    separate(sample,into = c("cond"),remove = FALSE) |>
    mutate(cond = ifelse(grepl("^ct",cond),"control",cond)) |>
    mutate(cond = ifelse(grepl("i",cond),"iERKBDNF",cond)) |>
    mutate(compart = ifelse(grepl("soma",sample),'soma','axon')) |>
    ggplot(aes(x = PC1, y = PC2,color = cond,shape = compart)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = sample)) +
    theme_bw()
# 
# axon_pca$rotation |>
#     as.data.frame() |>
#     tibble::rownames_to_column('sample') |>
#     separate(sample,into = c("cond"),remove = FALSE) |>
#     mutate(cond = ifelse(grepl("^ct",cond),"control",cond)) |> 
#     mutate(cond = ifelse(grepl("i",cond),"iERKBDNF",cond)) |> 
#     mutate(compart = ifelse(grepl("soma",sample),'soma','axon')) |> 
#     ggplot(aes(x = PC2, y = PC4,color = cond,shape = compart)) +
#     geom_point(size = 5) +
#     geom_text_repel(aes(label = sample)) +
#     theme_bw()
# 
# check_mew = axon_pca$x |> 
#     as.data.frame() |> tibble::rownames_to_column('protein') |> 
#     slice_min(PC2,n = 10) |> 
#     pull(protein)
# 
# 
# axon_pca_input |> 
#     tibble::rownames_to_column('protein') |>
#     filter(protein %in% check_mew) |> 
#     tibble::column_to_rownames('protein') |>
#     # dplyr::select(contains("axon")) |> 
#     pheatmap::pheatmap(scale = 'row')

# Using NAGuideR before PCA -----------------------------------------------
both_mean = raw_input |> 
    select(protein_ids,protein_names,"control_soma_a1", "control_soma_a2", "control_soma_a3", 
           "control_soma_a4", "control_soma_a5", "BDNF_soma_c1", "BDNF_soma_c2", 
           "BDNF_soma_c3","control_axon_b2", "control_axon_b3", "control_axon_b4", 
           "control_axon_b5", "BDNF_axon_d2", "BDNF_axon_d3", "BDNF_axon_d5") |> 
    na.omit() |> 
    melt() |> 
    summarize(m = mean(value)) |> pull(m)

axon_dt = raw_input %>%
    select(protein_ids,protein_names,"control_soma_a1", "control_soma_a2", "control_soma_a3", 
           "control_soma_a4", "control_soma_a5", "BDNF_soma_c1", "BDNF_soma_c2", 
           "BDNF_soma_c3","control_axon_b2", "control_axon_b3", "control_axon_b4", 
           "control_axon_b5", "BDNF_axon_d2", "BDNF_axon_d3", "BDNF_axon_d5") |> 
    select(protein_ids, protein_names, matches("_axon")) %>%
    mutate(NA_counts = rowSums(select(., control_axon_b2:BDNF_axon_d5) %>% is.na(), na.rm = TRUE)) %>%
    group_by(protein_ids) %>%
    dplyr::slice(which.min(NA_counts))  |> 
    filter(NA_counts <= 4) |> 
    # select(-matches('b3')) |> #high amount of Nuclear Pore Proteins for axon
    select(-NA_counts) |> 
    ungroup() |> 
    dplyr::relocate(protein_ids,protein_names) %>% 
    mutate(across(where(is.numeric), ~ ./median(.,na.rm = TRUE))) %>% 
    mutate(across(where(is.numeric), ~ .* both_mean)) %>% 
    mutate(across(where(is.numeric), ~ log2(.))) 


axon_dt |> 
    fwrite(here::here('data/axon_soma_proteomics/MEDIAN_NORMED_INPUT_NAGUIDER_AXON_TOTAL.csv'))

data.table(Samples = (colnames(axon_dt))) |> 
    filter(!(Samples %in% c("protein_ids","protein_names"))) |> 
    separate(Samples,remove = FALSE, into = c("Groups")) |> 
    fwrite(here::here('data/axon_soma_proteomics/samplesheet_correct_protein_axon_for_input_naguider.csv'))

# NAguideR::NAguideR_app()

soma_dt = raw_input %>%
    select(protein_ids,protein_names,"control_soma_a1", "control_soma_a2", "control_soma_a3", 
           "control_soma_a4", "control_soma_a5", "BDNF_soma_c1", "BDNF_soma_c2", 
           "BDNF_soma_c3","control_axon_b2", "control_axon_b3", "control_axon_b4", 
           "control_axon_b5", "BDNF_axon_d2", "BDNF_axon_d3", "BDNF_axon_d5") |> 
    select(protein_ids, protein_names, matches("_soma")) %>%
    mutate(NA_counts = rowSums(select(., control_soma_a1:BDNF_soma_c3) %>% is.na(), na.rm = TRUE)) %>%
    group_by(protein_ids) %>%
    dplyr::slice(which.min(NA_counts))  |> 
    filter(NA_counts <= 4) |> 
    # select(-matches('b3')) |> #high amount of Nuclear Pore Proteins for axon
    select(-NA_counts) |> 
    ungroup() |> 
    dplyr::relocate(protein_ids,protein_names) %>% 
    mutate(across(where(is.numeric), ~ ./median(.,na.rm = TRUE))) %>% 
    mutate(across(where(is.numeric), ~ .* both_mean)) %>% 
    mutate(across(where(is.numeric), ~ log2(.))) 
    
soma_dt |> 
    fwrite(here::here('data/axon_soma_proteomics/WITH_HISTONE_INPUT_SOMA_NAGUIDER.csv'))

data.table(Samples = colnames(soma_dt)) |> 
    filter(!(Samples %in% c("protein_ids","protein_names"))) |> 
    separate(Samples,remove = FALSE, into = c("Groups")) |> 
    fwrite(here::here('data/axon_soma_proteomics/samplesheet_correct_protein_soma_for_input_naguider.csv'))

# NAguideR::NAguideR_app()

# read in imputed values --------------------------------------------------
imputed_total_axon = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/All_Normed_first_imseqrob_AxonOriginal_Result_1721383562.17019.csv")
imputed_total_soma = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/IMPUTED_SOMA_NORMED_FIRST")


imputed_total_soma %>% 
    mutate(TotalSoma = rowSums(select_if(., is.numeric), na.rm = TRUE)) |> 
    left_join(imputed_total_axon %>% 
                  mutate(TotalAxon = rowSums(select_if(., is.numeric), na.rm = TRUE))) |> 
    ggplot(aes(x = TotalSoma,
               y = TotalAxon)) +
    geom_point() + 
    geom_abline() + 
    ggpubr::theme_pubr() + 
    geom_hline(linetype = 'dotted',yintercept = 140) 






# de - median- soma -----------------------------------------------
imputed_total_soma = imputed_total_soma %>% 
    column_to_rownames('V1')
soma_total_median_de = run_standard_limma(imputed_total_soma,baseline_grep = 'control',contrast_grep = 'BDNF',log2 = FALSE)

soma_total_median_de |> 
    mutate(new_column = sub('_', '|', un)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name))  |> 
    mutate(de_prot = case_when(p_value < 0.05 & log_fc > 0.5 ~ 1,
                               p_value < 0.05 & log_fc < -0.5 ~ 1,
                               TRUE ~ 0.2))    |>  
    mutate(plot_name = ifelse(de_prot == 1, gene_name, NA_character_)) |> 
    ggplot(aes(x = log_fc,
               y = -log10(p_value),
               alpha = de_prot)) + 
    geom_point() +
    ggpubr::theme_pubr() + 
    geom_hline(yintercept = -log10(0.05),linetype = 'dotted') + 
    geom_vline(xintercept = 0.5,linetype = "dotted") + 
    geom_vline(xintercept = -0.5,linetype = "dotted") + 
    ylab("-Log10 P-value") +
    xlab("Log2FC Protein Abundance BDNF v Control - Soma") + 
    geom_text_repel(aes(label = plot_name)) +
    theme(legend.position = 'none')


# de - median - all samples - axon - swapped ------------------------------

new_axon_imputed = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/All_Normed_first_imseqrob_AxonOriginal_Result_1721383562.17019.csv")

new_axon_imputed_f = new_axon_imputed %>%
    select(-matches("b1|d4|d1")) %>% 
    column_to_rownames('V1')
axon_pca = new_axon_imputed_f %>%
    prcomp()

axon_pca$rotation |>
    as.data.frame() |>
    tibble::rownames_to_column('sample') |>
    separate(sample,into = c("cond"),remove = FALSE) |>
    mutate(cond = ifelse(grepl("^ct",cond),"control",cond)) |>
    mutate(cond = ifelse(grepl("i",cond),"iERKBDNF",cond)) |>
    mutate(compart = ifelse(grepl("soma",sample),'soma','axon')) |>
    mutate(cond = fct_relevel(cond,'SwappedControl',"SwappedBDNF","iERKBDNF")) %>% 
    ggplot(aes(x = PC1, y = PC2,color = cond,shape = compart)) +
    geom_point(size = 5) +
    theme_bw() + 
    scale_color_manual(values = c("#F93FFF", "#1B7DFD","#EA8C55"))

erki_v_bdnf_swapped = run_standard_limma(new_axon_imputed_f,baseline_grep = 'SwappedBDNF',contrast_grep = 'iERKBDNF_',log2 = FALSE)
erki_v_control_swapped = run_standard_limma(new_axon_imputed_f,baseline_grep = 'SwappedControl',contrast_grep = 'iERKBDNF_',log2 = FALSE)

# de - median - all samples - soma ----------------------------------------

new_soma_imputed = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/SOMA_impseq_Original_Result_1721384471.75033.csv")
new_soma_imputed_f = new_soma_imputed %>%
    select(-matches("c5|c4")) |>
    column_to_rownames('V1')
soma_pca = new_soma_imputed %>%
    select(-matches("c5|c4")) |>
    column_to_rownames('V1') %>%
    prcomp()



soma_pca$rotation |>
    as.data.frame() |>
    tibble::rownames_to_column('sample') |>
    separate(sample,into = c("cond"),remove = FALSE) |>
    mutate(cond = ifelse(grepl("^ct",cond),"control",cond)) |>
    mutate(cond = ifelse(grepl("i",cond),"iERKBDNF",cond)) |>
    mutate(compart = ifelse(grepl("soma",sample),'soma','axon')) |>
    mutate(cond = fct_relevel(cond,'control',"BDNF","iERKBDNF")) %>% 
    ggplot(aes(x = PC1, y = PC2,color = cond,shape = compart)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = sample)) +
    theme_bw() +
    scale_color_manual(values = c("#F93FFF", "#1B7DFD","#EA8C55"))




erki_v_bdnf_soma = run_standard_limma(new_soma_imputed_f,baseline_grep = '^BDNF',contrast_grep = 'iERKBDNF_',log2 = FALSE)

erki_v_control_soma = run_standard_limma(new_soma_imputed_f,baseline_grep = 'control',contrast_grep = 'iERKBDNF_',log2 = FALSE)

bdnf_v_control_soma = run_standard_limma(new_soma_imputed_f,baseline_grep = 'control',contrast_grep = '^BDNF',log2 = FALSE)


    

# de - median- axon -----------------------------------------------
amigoingmad()
imputed_total_axon = imputed_total_axon %>% 
    column_to_rownames('V1')

axon_total_median_de = run_standard_limma(imputed_total_axon,baseline_grep = 'control',contrast_grep = '^BDNF',log2 = FALSE) 

axon_total_median_de |> 
    filter(un %in% de_bdnf_de_erki) %>% 
    mutate(new_column = sub('_', '|', un)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name))  |> 
    mutate(de_prot = case_when(adj_p_val < 0.1 & log_fc > 0.5 ~ 1,
                               adj_p_val < 0.1 & log_fc < -0.5 ~ 1,
                               TRUE ~ 0.2))    |>  
    mutate(plot_name = ifelse(de_prot == 1, gene_name, NA_character_)) |> 
    unique() %>% 
    ggplot(aes(x = log_fc,
               y = -log10(p_value),
               alpha = de_prot)) + 
    geom_point() +
    ggpubr::theme_pubr() + 
    geom_hline(yintercept = -log10(0.05),linetype = 'dotted') + 
    geom_vline(xintercept = 0.5,linetype = "dotted") + 
    geom_vline(xintercept = -0.5,linetype = "dotted") + 
    ylab("-Log10 P-value") +
    xlab("Log2FC Protein Abundance BDNF v Control - Axon") + 
    geom_text_repel(aes(label = plot_name)) + 
    theme(legend.position = 'none')



# GO soma BDNF effect ---------------------------------------------------
amigoingmad()
 


soma_bdnf = clusterProfiler::compareCluster(target~protein_effect_soma,
                                            data = soma_go_dt %>% 
                                                filter(protein_effect_soma != "NS"),
                                            universe = bg_expressed_in_LMN,
                                            keyType = 'ENSEMBL',
                                            OrgDb = org.Hs.eg.db,
                                            ont = 'ALL',
                                            readable = TRUE)

dotplot(soma_bdnf, x="protein_effect_soma") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(element_blank()) +
    theme(strip.text = element_text(size = 18)) +
    theme(axis.text.y  = element_text(size = 9))


# GO both -----------------------------------------------------------------

# explording the genes leading to downregulation of ribosome --------------

axon_bdnf@compareClusterResult %>% 
    as.data.table() %>% 
    filter(Cluster == 'downregulated')

# axon_go_dt2 = axon_total_median_de |> 
#     separate(un,into = c("uniprot"),remove = FALSE) %>% 
#     distinct(uniprot,log_fc,p_value) %>%  
#     mutate(protein_axon = case_when(log_fc > 0.5 & p_value < 0.05 ~ "upregulated",
#                                     log_fc < -0.5 & p_value < 0.05 ~ "downregulated",
#                                     TRUE ~ "NS")) 
# 
# axon_bdnf2 = clusterProfiler::compareCluster(target~protein_axon,data = axon_go_dt2 %>% filter(protein_axon != "NS"),
#                                             universe = mapped_uniprot,
#                                             keyType = 'UNIPROT',
#                                             OrgDb = org.Hs.eg.db,
#                                             ont = 'ALL',
#                                             readable = TRUE)

# GSEA axon BDNF effect ---------------------------------------------------
axon_gsea_dt = axon_go_dt |> 
    distinct(target,log_fc,p_value) %>% 
    mutate(fcsign = sign(log_fc),
           logP = -log10(p_value),
           metric= logP/fcsign) %>% 
    left_join(annotables::grch38 %>% select(ensgene,entrez),by = c("target" = 'ensgene')) %>% 
    filter(!is.na(entrez))


## feature 1: numeric vector
geneList = axon_gsea_dt$metric

## feature 2: named vector
names(geneList) = as.character(axon_gsea_dt$target)

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)

ego4 <- gseGO(geneList     = geneList,
              ont          = "MF",
              OrgDb='org.Hs.eg.db',
              keyType = "ENSEMBL",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.9,
              verbose      = FALSE)

ridgeplot(ego4,  showCategory = 20)

names(geneList) = as.character(axon_gsea_dt$entrez)
ak = gseKEGG(geneList, organism = "hsa", 
             keyType = "ncbi-geneid", 
             exponent = 1, nPerm = 1000, 
             minGSSize = 10, maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             pAdjustMethod = "BH", 
             verbose = TRUE, use_internal_data = FALSE, seed = FALSE)

ridgeplot(ak,  showCategory = 20)


# GSEA soma BDNF effect ---------------------------------------------------
soma_gsea_dt = soma_total_median_de |> 
    separate(un,into = c("uniprot"),remove = FALSE) %>% 
    add_the_annotation('uniprot') %>% 
    distinct(target,log_fc,p_value) %>% 
    mutate(fcsign = sign(log_fc),
           logP = -log10(p_value),
           metric= logP/fcsign) %>% 
    left_join(annotables::grch38 %>% select(ensgene,entrez),by = c("target" = 'ensgene')) %>% 
    filter(!is.na(entrez))

# compare soma and axon BDNF effect ---------------------------------------
soma_total_median_de |> 
    full_join(axon_total_median_de,by = c("un"),suffix = c("_soma",
                                                           "_axon")) |> 
    na.omit() |> 
    # filter(p_value_soma < 0.05 & p_value_axon < 0.05)  |> 
    # filter(abs(log_fc_soma) > 0.5) |> 
    # filter(abs(log_fc_axon) > 0.5) |> 
    mutate(new_column = sub('_', '|', un)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name))  |> 
    mutate(plot_name = ifelse(abs(log_fc_soma) > 0.5 & abs(log_fc_axon) > 0.5,gene_name,NA_character_)) %>% 
    # filter(genes %in% these_genes) |> 
    # left_join(just_ids) |>
    filter(p_value_axon < 0.05 & p_value_soma < 0.05) %>% 
    ggplot(aes(x = log_fc_soma,
               y = log_fc_axon)) + 
    geom_text_repel(aes(label = gene_name)) +
    geom_point() +
    geom_abline() + 
    ggpubr::theme_pubr() + 
    geom_hline(yintercept = 0,linetype = 'dotted') + 
    geom_vline(xintercept = 0,linetype = "dotted") + 
    ylab("Log2FC Protein Abundance BDNF v Control - Axon") +
    xlab("Log2FC Protein Abundance BDNF v Control - Soma")

soma_total_median_de |> 
    full_join(axon_total_median_de,by = c("un"),suffix = c("_soma",
                                                           "_axon")) |> 
    na.omit() |> 
    filter(p_value_soma < 0.05 & p_value_axon < 0.05)  |>
    filter(abs(log_fc_soma) > 0.5) |>
    filter(abs(log_fc_axon) > 0.5) |>
    mutate(new_column = sub('_', '|', un)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name))  |> 
    left_join(tmp) %>% 
    lef
    # left_join(just_ids) |>
    # filter(p_value_axon < 0.05 & p_value_soma < 0.05) |> 
    ggplot(aes(x = log_fc_soma,
               y = log_fc_axon)) + 
    geom_text_repel(aes(label = gene_name)) +
    geom_point() +
    geom_abline() + 
    ggpubr::theme_pubr() + 
    geom_hline(yintercept = 0,linetype = 'dotted') + 
    geom_vline(xintercept = 0,linetype = "dotted") + 
    ylab("Log2FC Protein Abundance BDNF v Control - Axon") +
    xlab("Log2FC Protein Abundance BDNF v Control - Soma")



# Top 10% Soma & Axon Control------------------------------------------------------------
minMax <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}


soma_dt = fread("data/axon_soma_proteomics/WITH_HISTONE_INPUT_SOMA_NAGUIDER.csv")
axon_dt = fread(here::here('data/axon_soma_proteomics/ALL_SAMPLES_MEDIAN_NORMED_INPUT_NAGUIDER_AXON_TOTAL.csv'))

reliably_detected = soma_dt %>% 
    full_join(axon_dt,by = c("protein_ids","protein_names"),
              suffix = c("_soma","_axon")) %>% 
    select(-matches('BDNF')) %>%
    filter(!grepl("Q5T750",protein_ids)) %>% 
    mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
    melt() %>% 
    separate(variable,into = c("cond",'compart','samp')) %>% 
    mutate(detected = value > 0) %>% 
    group_by(protein_ids,protein_names,compart) %>% 
    summarize(n_detect = sum(detected)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = 'compart',
                values_from = 'n_detect') %>% as.data.table() %>% 
    filter(axon >=3 | soma>=3) %>% 
    mutate(un = glue::glue('{protein_ids}_{protein_names}')) 
    


raw_axon_ranks_ctrl = axon_dt %>% 
    filter(protein_ids %in% reliably_detected$protein_ids) %>% 
    select(-matches('BDNF')) %>%
    # select(-matches('control')) %>% 
    mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% 
    arrange(-Total) %>% 
    mutate(rank_axon = rank((Total))) %>% 
    mutate(mmRank = minMax(Total))

raw_soma_ranks_ctrl = soma_dt %>% 
    filter(protein_ids %in% reliably_detected$protein_ids) %>% 
    select(-matches('BDNF')) %>%
    # select(-matches('control')) %>% 
    mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% 
    arrange(-Total) %>% 
    mutate(rank_axon = rank((Total))) %>% 
    mutate(mmRank = minMax(Total))


bth_rank_ctrl = soma_dt %>% 
    full_join(axon_dt,by = c("protein_ids","protein_names"),
              suffix = c("_soma","_axon")) %>% 
    select(-matches('BDNF')) %>%
    mutate(across(everything(), .fns = ~replace_na(.,0))) %>% 
    filter(protein_ids %in% reliably_detected$protein_ids) %>% 
    # select(-matches('BDNF')) %>%
    mutate(TotalSoma = rowSums(select(., contains("soma")), na.rm = TRUE)) %>% 
    mutate(TotalAxon = rowSums(select(., contains("axon")), na.rm = TRUE)) %>% 
    mutate(mmRank_soma = minMax(TotalSoma)) %>% 
    mutate(mmRank_axon = minMax(TotalAxon)) %>% 
    as.data.table()



n_soma_only = reliably_detected[axon == 0]$un %>% n_distinct()
n_axon_only = reliably_detected[soma == 0]$un %>% n_distinct()

bth_rank_ctrl %>% 
    ggplot(aes(x = mmRank_soma,
               y = mmRank_axon)) + 
    # geom_point(alpha = 0.3) + 
    geom_hex() +
    geom_abline() + 
    ylab("Axon normalised signal intensity rank") +
    xlab("Soma normalised signal intensity rank") +
    ggpubr::theme_pubr() +
    scale_fill_gradient(low="grey62",high="black") 
    


top_15_axon = raw_axon_ranks_ctrl %>%
    slice_max(Total,prop = 0.15) |>
    mutate(un = glue::glue('{protein_ids}_{protein_names}')) 


top_15_soma = raw_soma_ranks_ctrl %>%
    slice_max(Total,prop = 0.15) |>
    mutate(un = glue::glue('{protein_ids}_{protein_names}')) 



reliably_detected %>% 
    mutate(un = glue::glue('{protein_ids}_{protein_names}')) %>% 
    mutate(top15_axon = un %in% top_15_axon$un) %>% 
    mutate(top_15_soma = un %in% top_15_soma$un) %>% 
    mutate(detected_axon = axon >=3,
           detected_soma = soma >= 3) %>% 
    select(un,top_15_soma,top15_axon,detected_axon,detected_soma) %>% 
    unique() %>% 
    column_to_rownames('un') %>% 
    eulerr::euler()   |> 
    plot(quantities = TRUE)



tmp = reliably_detected %>% 
    separate_rows(protein_ids,sep = ";") %>% 
    add_the_annotation('protein_ids') %>% 
    mutate(rbp = target %in% rbps$x2) %>% 
    mutate(kinases = target %in% kinases$target) %>% 
    mutate(phophotases = target %in% phophotases$target) %>% 
    mutate(e3_ligases = target %in% e3_ligases$target) %>% 
    mutate(top15_axon = un %in% top_15_axon$un) %>% 
    mutate(top_15_soma = un %in% top_15_soma$un) %>% 
    distinct(un,target,rbp,kinases,phophotases,e3_ligases,top15_axon,top_15_soma,axon,soma) 


# GO on the top15 soma which are NOT top15 axon --------------------------------------
top_soma_not_top_axon = reliably_detected %>% 
    filter(un %in% top_15_soma$un) %>% 
    filter(!(un %in% top_15_axon$un)) %>% 
    separate(protein_ids,into = c("uniprot")) %>% 
    add_the_annotation('uniprot') %>% 
    unique()
go_top_soma_no_axon = clusterProfiler::enrichGO(top_soma_not_top_axon$target,
                                                universe = bg_expressed_in_LMN,
                                                keyType = 'ENSEMBL',
                                                OrgDb = org.Hs.eg.db,
                                                ont = 'MF',
                                                readable = TRUE)
dotplot(go_top_soma_no_axon)
cnetplot(go_top_soma_no_axon)
# GO on the top15 axon which are NOT top15 soma --------------------------------------
top_axon_not_top_soma = reliably_detected %>% 
    filter(un %in% top_15_axon$un) %>% 
    filter(!(un %in% top_15_soma$un)) %>% 
    separate(protein_ids,into = c("uniprot")) %>% 
    add_the_annotation('uniprot') %>% 
    unique()
go_top_axon_no_soma = clusterProfiler::enrichGO(top_axon_not_top_soma$target,
                                                universe = bg_expressed_in_LMN,
                                                keyType = 'ENSEMBL',
                                                OrgDb = org.Hs.eg.db,
                                                ont = 'ALL',
                                                readable = TRUE)
dotplot(go_top_axon_no_soma)
cnetplot(go_top_axon_no_soma)



# expression of top vs not ------------------------------------------------


raw_soma_ranks_ctrl %>%
    left_join(reliably_detected) %>% 
    mutate(detection = case_when(axon >=3 & soma >=3 ~  "Reliablly detected soma & axon",
                                 soma >=3 ~ "Reliably detected soma only")) %>% 
    filter(!is.na(detection)) %>% 
    ggplot(aes(x = detection, y = mmRank)) + 
    geom_boxplot() +
    ylab("Soma normalised signal intensity rank") + 
    theme_classic() + 
    ggpubr::stat_compare_means()


raw_axon_ranks_ctrl %>%
    left_join(reliably_detected) %>% 
    mutate(detection = case_when(axon >=3 & soma >=3 ~  "Reliablly detected soma & axon",
                                 axon >=3 ~ "Reliably detected axon only")) %>% 
    filter(!is.na(detection)) %>% 
    ggplot(aes(x = detection, y = mmRank)) + 
    geom_boxplot() + 
    ylab("Axonal normalised signal intensity rank") +
    theme_classic() + 
    ggpubr::stat_compare_means()

# PCA after imputation AXON----------------------------------------------------
imputed_total_axon = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/All_Normed_first_imseqrob_AxonOriginal_Result_1721383562.17019.csv")
imputed_total_axon = imputed_total_axon %>%
         select(-matches("b1|d4|d1")) %>% 
         column_to_rownames('V1')
axon_p = imputed_total_axon |> 
    prcomp()

v1 = summary(axon_p)$importance |> as.data.frame() |> slice(2) |> pull(PC1)
v2 = summary(axon_p)$importance |> as.data.frame() |> slice(2) |> pull(PC2)


axon_p$rotation |>
    as.data.frame() |>
    tibble::rownames_to_column('sample') |>
    separate(sample,into = c("cond","compart",'sample'),remove = FALSE) |> 
    mutate(cond = fct_relevel(cond,'control')) |> 
    ggplot(aes(x = PC1, y = PC2,fill = cond,shape = compart)) +
    geom_point(size = 5) +
    theme_bw() + 
    xlab(paste0("PC1", ', ', round(v1 * 100, digits = 2), '% variation')) + 
    ylab( paste0("PC2", ', ', round(v2 * 100, digits = 2), '% variation')) +
    scale_shape_manual(values = c('axon'=21, 'soma'=24)) +
    scale_fill_manual(values = c("#F93FFF", "#1B7DFD","#EA8C55"))



top_pc1 = axon_p$x |> 
    as.data.frame() |> tibble::rownames_to_column('protein') |> 
    slice_max(PC1,n = 20) |> 
    pull(protein)
bottom_pc1 = axon_p$x |> 
    as.data.frame() |> tibble::rownames_to_column('protein') |> 
    slice_min(PC1,n = 20) |> 
    pull(protein)
top_pc2 = axon_p$x |> 
    as.data.frame() |> tibble::rownames_to_column('protein') |> 
    slice_max(PC2,n = 20) |> 
    pull(protein)
bottom_pc2 = axon_p$x |> 
    as.data.frame() |> tibble::rownames_to_column('protein') |> 
    slice_min(PC2,n = 20) |> 
    pull(protein)

imputed_total_axon |> 
    rownames_to_column('V1') %>% 
    filter(V1 %in% top_pc1) |> 
    tibble::column_to_rownames('V1') %>% 
    pheatmap::pheatmap(scale = 'row')


# PCA after imputation SOMA -----------------------------------------------
imputed_total_soma = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/axon_soma_proteomics/SOMA_impseq_Original_Result_1721384471.75033.csv")
imputed_total_soma = new_soma_imputed %>%
    select(-matches("c5|c4")) |>
    column_to_rownames('V1')

soma_p = imputed_total_soma |> 
    prcomp()

v1 = summary(soma_p)$importance |> as.data.frame() |> slice(2) |> pull(PC1)
v2 = summary(soma_p)$importance |> as.data.frame() |> slice(2) |> pull(PC2)


soma_p$rotation |>
    as.data.frame() |>
    tibble::rownames_to_column('sample') |>
    separate(sample,into = c("cond","compart",'sample'),remove = FALSE) |> 
    mutate(cond = fct_relevel(cond,'control')) |> 
    ggplot(aes(x = PC1, y = PC2,fill = cond,shape = compart)) +
    geom_point(size = 5) +
    theme_bw() + 
    xlab(paste0("PC1", ', ', round(v1 * 100, digits = 2), '% variation')) + 
    ylab( paste0("PC2", ', ', round(v2 * 100, digits = 2), '% variation')) +
    scale_shape_manual(values = c('axon'=21, 'soma'=24)) +
    scale_fill_manual(values = c("#F93FFF", "#1B7DFD","#EA8C55"))


# PCA after imputation BOTH -----------------------------------------------
total_p = imputed_total_soma |> 
    left_join(imputed_total_axon,by = c("V1")) |> 
    na.omit()  |>
    select(-matches("b1|d4|d1")) %>%
    select(-matches("c5|c4")) |>
    select(V1,matches("control")) %>% 
    # mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
    tibble::column_to_rownames("V1") |> 
    prcomp()

v1 = summary(total_p)$importance |> as.data.frame() |> slice(2) |> pull(PC1)
v2 = summary(total_p)$importance |> as.data.frame() |> slice(2) |> pull(PC2)


total_p$rotation |>
    as.data.frame() |>
    tibble::rownames_to_column('sample') |>
    separate(sample,into = c("cond","compart",'sample'),remove = FALSE) |> 
    mutate(cond = fct_relevel(cond,'control')) |> 
    ggplot(aes(x = PC1, y = PC2,fill = cond,shape = compart)) +
    geom_point(size = 5) +
    theme_bw() + 
    scale_fill_manual(values = c(control = "#F93FFF", BDNF = "#1B7DFD",iERKBDNF = '#EA8C55')) + 
    xlab(paste0("PC1", ', ', round(v1 * 100, digits = 2), '% variation')) + 
    ylab( paste0("PC2", ', ', round(v2 * 100, digits = 2), '% variation')) +
    scale_shape_manual(values = c('axon'=21, 'soma'=24))


# Heatmap PC1 & PC2 BOTH-------------------------------------------------------


top_pc1 = total_p$x |> 
    as.data.frame() |> tibble::rownames_to_column('protein') |> 
    slice_max(PC1,n = 15) |> 
    pull(protein)

bottom_pc1 = total_p$x |> 
    as.data.frame() |> tibble::rownames_to_column('protein') |> 
    slice_min(PC1,n = 15) |> 
    pull(protein)



# Data processing pipeline
processed_data <- imputed_total_soma %>%
    left_join(imputed_total_axon, by = "V1") %>%
    select(-matches("b1|d4|d1")) %>%
    select(-matches("c5|c4")) |>
    select(V1,matches("control")) |>
    
    na.omit() %>%
    filter(V1 %in% c(bottom_pc1, top_pc1)) %>%
    mutate(new_column = sub('_', '|', V1)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name)) %>%
    column_to_rownames('gene_name') %>%
    select(-V1)

# Create annotation data
annotation_data <- imputed_total_soma %>%
    filter(V1 %in% c(bottom_pc1, top_pc1)) %>%
    mutate(annotation_label = ifelse(V1 %in% bottom_pc1, "Bottom PC1", "Top PC1")) %>%
    select(V1, annotation_label) %>%
    mutate(new_column = sub('_', '|', V1)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name)) %>%
    column_to_rownames('gene_name') %>%
    select(-V1)

# Ensure the row order of annotation matches the processed data
annotation_data <- annotation_data[rownames(processed_data), , drop = FALSE]

# Sort the rows based on the annotation_label
sorted_indices <- order(annotation_data$annotation_label, decreasing = TRUE)
processed_data <- processed_data[sorted_indices, ]
annotation_data <- annotation_data[sorted_indices, ]

# Scaling the data
scaled_data <- t(apply(processed_data, 1, scale))

# Create the row annotation
row_anno <- rowAnnotation(
    pc1_label = factor(annotation_data, levels = c("Top PC1", "Bottom PC1")),
    col = list(pc1_label = c("Bottom PC1" = "blue", "Top PC1" = "red")),
    show_annotation_name = FALSE  # This hides the annotation label
)

# Create column annotation
column_info <- tibble(colname = colnames(processed_data)) %>%
    separate(colname, into = c("condition", "compartment", "sample"), sep = "_", extra = "merge")

# Create the column annotation object
col_anno <- columnAnnotation(
    condition = column_info$condition,
    compartment = column_info$compartment,
    col = list(
        condition = c("control" = "#F93FFF", "BDNF" = "#1B7DFD",iERKBDNF = '#72f54aff'),
        compartment = c("soma" = "grey", "axon" = "black")
    )
)

# Creating the heatmap with both row and column annotations, and no row clustering
Heatmap(scaled_data, name = "expression", 
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = FALSE,  # No clustering on rows
        cluster_columns = TRUE, 
        left_annotation = row_anno,
        top_annotation = col_anno)


# HEATMAP - opposing and same ---------------------------------------------
sig2 = soma_imputed_results_total |> 
    left_join(axon_imputed_results_total,by = c("un"),suffix = c("_soma",
                                                                 "_axon")) |> 
    na.omit() |> 
    filter(p_value_soma < 0.05 & p_value_axon < 0.05)  |> 
    filter(abs(log_fc_soma) > 0.5) |> 
    filter(abs(log_fc_axon) > 0.5) |> 
    mutate(new_column = sub('_', '|', un)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name)) 

inter = sig2 %>%    # Rank by log_fc_soma and log_fc_axon
    mutate(rank_soma = rank(-abs(log_fc_soma)), 
           rank_axon = rank(-abs(log_fc_axon))) |> 
    mutate(keepReason = case_when(rank_soma <=10 ~ "TopSoma",
                                  rank_axon <= 10 ~ "TopAxon",
                                  (log_fc_soma * log_fc_axon) > 0 ~ "SameRegulation")) |> 
    filter(!is.na(keepReason)) |> 
    mutate(new_column = sub('_', '|', un)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name)) 

processed_data <- imputed_total_soma %>%
    left_join(imputed_total_axon, by = "V1") %>%
    mutate(new_column = sub('_', '|', V1)) %>%
    separate(new_column, into = c(NA, "gene_name"), sep = "\\|") %>%
    mutate(gene_name = gsub("_HUMAN", "", gene_name)) %>%
    na.omit() %>%
    filter(gene_name %in% inter$gene_name) %>%
    select(-V1) |> 
    column_to_rownames('gene_name')
# Create annotation data
annotation_data <- inter %>%
    select(gene_name, keepReason) %>%
    column_to_rownames('gene_name')
# Ensure the row order of annotation matches the processed data
annotation_data <- annotation_data[rownames(processed_data), , drop = FALSE]

# Set the factor levels and order by keepReason
annotation_data$keepReason <- factor(annotation_data$keepReason, levels = c("SameRegulation", "TopAxon", "TopSoma"))
sorted_indices <- order(annotation_data$keepReason)
processed_data <- processed_data[sorted_indices, ]
annotation_data <- annotation_data[sorted_indices, ]

# Scaling the data
scaled_data <- t(apply(processed_data, 1, scale))
scaled_data = processed_data
# Step 5: Create the row annotation object
row_anno <- rowAnnotation(
    keepReason = factor(annotation_data, levels = c("SameRegulation", "TopAxon", "TopSoma")),
    col = list(keepReason = c("SameRegulation" = "green", "TopAxon" = "red", "TopSoma" = "blue")),
    show_annotation_name = FALSE  # This hides the annotation label
)

# Create column annotation
column_info <- tibble(colname = colnames(processed_data)) %>%
    separate(colname, into = c("condition", "compartment", "sample"), sep = "_", extra = "merge")

# Create the column annotation object
col_anno <- columnAnnotation(
    condition = column_info$condition,
    compartment = column_info$compartment,
    col = list(
        condition = c("control" = "#F93FFF", "BDNF" = "#1B7DFD"),
        compartment = c("soma" = "grey", "axon" = "black")
    )
)
# Creating the heatmap with both row and column annotations, and no row clustering
Heatmap(scaled_data, name = "expression", 
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        cluster_rows = FALSE,  # No clustering on rows
        cluster_columns = FALSE, 
        left_annotation = row_anno,
        top_annotation = col_anno)





# plot 1433 expression --------------------------------------------------

imputed_total_axon[V1 == 'P31946;P61981_1433B_HUMAN;1433G_HUMAN'] %>% 
    select(-matches("b1|d4|d1")) %>%
    select(-matches("c5|c4")) |>
    melt() %>% 
    separate(variable,into = c("don")) %>% 
    mutate(don = fct_relevel(don,'control'))  %>% 
    ggplot(aes(x = don, y = value,)) + 
    geom_boxplot(outlier.color = 'NA',aes(fill = don)) +
    geom_jitter(height = 0,size = 2,aes(fill = don),pch = 21) + 
    theme_classic() +
    scale_fill_manual(values = c("#F93FFF", "#1B7DFD","#72f54aff")) + 
    ylab("1433 protein intensity") + 
    ggtitle('Axonal compartment') + 
    xlab(element_blank()) + 
    theme(legend.position = 'none') + facet_wrap(~V1)




imputed_total_soma[V1 == 'P31946;P61981_1433B_HUMAN;1433G_HUMAN'] %>% 
    select(-matches("b1|d4|d1")) %>%
    select(-matches("c5|c4")) |>
    melt() %>% 
    separate(variable,into = c("don")) %>% 
    mutate(don = fct_relevel(don,'control')) %>% 
    ggplot(aes(x = don, y = value)) + 
    geom_boxplot(outlier.color = 'NA',aes(fill = don)) +
    geom_jitter(height = 0,size = 2,aes(fill = don),pch = 21) + 
    theme_classic() +
    scale_fill_manual(values = c("#F93FFF", "#1B7DFD","#72f54aff")) + 
    ylab("1433 protein intensity") + 
    ggtitle('Soma compartment') +
    xlab(element_blank()) + 
    theme(legend.position = 'none')


