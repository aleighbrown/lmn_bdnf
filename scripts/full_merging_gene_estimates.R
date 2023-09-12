
source(here::here('new_ratio_functions.R'))

file_path = here::here("data", "full_run_grandslam")
suffix = "hr.tsv"
estimate_files = list.files(file_path,
                            pattern = suffix,
                            full.names = TRUE)

estimate_list_full = purrr::map(estimate_files,my_clean_reader)
samp_ids = base::tolower(purrr::simplify(purrr::map(estimate_files, basename)))
samp_ids = gsub(".tsv","",samp_ids)
##add on the name as an additional column
estimate_list_full = purrr::map2(estimate_list_full, samp_ids, ~cbind(.x, SampleID = .y))
estimate_list_full = purrr::map(estimate_list_full,fix_column_names)
###use the data.table rbindlist to make a long format
estimate_list_full = data.table::rbindlist(estimate_list_full)
estimate_list_full[,time := parse_number(time)]
estimate_list_full[,well := parse_number(well)]
names(estimate_list_full) = gsub("calmdsorted","",names(estimate_list_full))
estimate_list_full = estimate_list_full %>% 
    mutate(map_range = abs(`005quantile` - `095quantile`))

estimate_list_full = estimate_list_full[well != 4]



# calculating decay rate --------------------------------------------------

path = "/Users/annaleigh/Documents/GitHub/bdnf_4su/data/full_run_grandslam/"
grandslamfiles = list.files(path,
                            pattern = "hr.tsv",
                            full.names = TRUE)

grandslam_alpha = sapply(grandslamfiles, my_clean_reader, simplify = FALSE) %>% 
    purrr::reduce(., full_join, by = c("gene","symbol","length")) %>% 
    dplyr::select(contains("alpha"),contains("readcount"),contains("beta"),contains("map"),gene,symbol)

grandslam_alpha = grandslam_alpha[complete.cases(grandslam_alpha),] 

grandslam_alpha = grandslam_alpha %>% 
    data.table::melt(id.vars = c("gene","symbol")) %>% 
    mutate(variable = gsub("_calmd_sorted_","_",variable)) %>% 
    separate(variable, into = c("condition","sample","time","measure")) %>% 
    pivot_wider(names_from = "measure",values_from = "value") %>% 
    mutate(time = parse_number(time)) %>% as.data.table()


grandslam_alpha = grandslam_alpha %>% 
    group_by(condition,gene,time) %>%
    mutate(n_samples = n_distinct(sample)) %>% 
    filter(n_samples == 3) %>% 
    ungroup() %>% 
    as.data.table()


grandslam_alpha[,decay := mle_decay(alpha, beta,time), by = c("gene","condition")]


confidence_decay = grandslam_alpha %>% 
    group_by(gene,condition) %>% 
    nest() %>% 
    mutate(ci = map(data, getCI)) %>% 
    unnest_wider(ci) %>% 
    dplyr::select(-data)



half_lifes <- grandslam_alpha %>% 
    dplyr::select(gene,symbol,condition,decay) %>% 
    unique()

half_lifes = half_lifes %>% 
    left_join(confidence_decay, by = c("condition","gene")) %>% 
    mutate(half_life_map = log(2)/decay) %>% 
    mutate(half_life_upperCI = log(2)/lower) %>% 
    mutate(half_life_lowerCI = log(2)/upper) 


half_lifes %>% 
    filter(symbol == "CYP11B1") %>% 
    ggplot(aes(x = condition,
               y = half_life_map)) + 
    geom_pointrange(aes(ymin = half_life_lowerCI,
                        ymax = half_life_upperCI))

wide_lives = half_lifes %>%
    pivot_wider(id_cols = c("gene","symbol"),
                names_from = "condition",
                values_from = c("half_life_map",
                                "half_life_upperCI",
                                "half_life_lowerCI")) %>%
    as.data.table() 

wide_lives[,ci_range_control := abs(half_life_lowerCI_control - half_life_upperCI_control)]
wide_lives[,ci_range_bdnf := abs(half_life_upperCI_bdnf - half_life_lowerCI_bdnf)]

ranges = wide_lives[,.(half_life_lowerCI_bdnf,half_life_upperCI_bdnf,
                       half_life_lowerCI_control,half_life_upperCI_control)]

ci_dont_overlap = !((pmin(ranges[,1], ranges[,2]) <= pmax(ranges[,3], ranges[,4])) &
    (pmax(ranges[,1], ranges[,2]) >= pmin(ranges[,3], ranges[,4]))) %>% as.logical()


wide_lives %>% 
    mutate(log2halflife = log2(half_life_map_bdnf / half_life_map_control)) %>% 
    filter(ci_dont_overlap) %>% 
    arrange(-log2halflife) %>% View()


sig_half_life = wide_lives %>% 
    filter(ci_dont_overlap) %>% 
    mutate(bdnf_ci = abs(half_life_lowerCI_bdnf - half_life_upperCI_bdnf)) %>% 
    mutate(control_ci = abs(half_life_lowerCI_control - half_life_upperCI_control)) %>% 
    filter(bdnf_ci < 2 * half_life_map_bdnf) %>% 
    filter(control_ci < 2 * half_life_map_control) %>% 
    mutate(log2halflife = log2(half_life_map_bdnf / half_life_map_control)) %>% 
    arrange(-log2halflife)



# trying to compare new rna levesl at each time point ---------------------

tmp = grandslam_alpha %>% group_by(time,gene) %>% nest()

tmp2 = tmp %>%  mutate(values = map(data, bayesian_p_value)) 


bayesian_times = tmp2 %>% select(-data) %>% unnest()

n_samp_passing = estimate_list_full %>% 
    filter(map_range < 0.2) %>% 
    dplyr::select(time,gene,condition,SampleID) %>% 
    group_by(time,gene,condition) %>% mutate(n_samp = n_distinct(SampleID)) %>% 
    ungroup() %>% 
    dplyr::select(-SampleID) %>% 
    unique() %>% 
    pivot_wider(names_from = "condition",
                values_from = "n_samp",
                id_cols = c("time","gene"),
                names_prefix = "n_samp_passing_")


new_ratio_log2Fold = bayesian_times %>% 
    ungroup() %>% 
    left_join(n_samp_passing) %>% 
    dplyr::select(gene,time,bayesian_p,mean_diff,mean_bdnf_ntr,mean_control_ntr,
                  n_samp_passing_bdnf,n_samp_passing_control) %>% 
    left_join(full_de, by = c("gene" = "Geneid","time")) %>% 
    dplyr::select(-lfcSE,stat,pvalue)

# If 90% of  the sampled new RNA estimates was higher in the BDNF compared to Control, 
# then a gene was categorized as ‘bdnf_higher_new_rna’, if less than 10% of samples was higher in BDNF, 
# then a gene was categorized as ‘bdnf_lower_new_rna’, if the BDNF estimate was higher in 40% to 60% of the sampled new RNA estimates, 
# the gene was categorized as ‘bdnf_equals_control’. 
# For all categorizations the gene must have met sample passing requirements for 2 out of 3 replicates in both BDNF and Control conditions (n_samp_passing_bdnf >= 2 & n_samp_passing_control >= 2). Any gene that did not fall into this criteria were marked as ‘unclear’. 

new_ratio_log2Fold = new_ratio_log2Fold %>% 
    mutate(new_rna_sig = case_when(bayesian_p > 0.9 & 
                                       n_samp_passing_bdnf >= 2 & 
                                       n_samp_passing_control >=2 ~ "bdnf_lower_new_rna",
                                   bayesian_p < 0.1 & 
                                       n_samp_passing_bdnf >= 2 & 
                                       n_samp_passing_control >=2 ~ "bdnf_higher_new_rna",
                                   bayesian_p < 0.1 & 
                                       n_samp_passing_bdnf >= 2 & 
                                       n_samp_passing_control >=2 ~ "bdnfs_equal_controls",
                                   T ~ "unclear"
                                       )) %>% 
    mutate(total_rna_sig = case_when(padj < 0.1 & log2FoldChange > 0.75 ~ "upregulated",
                                     padj < 0.1 & log2FoldChange < -0.75 ~ "downregulated",
                                     T ~ "not_significant"))
