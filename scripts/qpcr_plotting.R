library(DESeq2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)



#FOS total and new, ACTD - 1 hour treatment

#DCTN2/1I2 total and qPCR translation
qpcr_data = fread("data/qpcr_results_bdnf.csv") %>% 
    filter(gene == 'DYNC1I2')

estimate_list_full = fread('data/estimate_list_full.csv')

qpcr_data = qpcr_data %>% 
    mutate(condition = fct_relevel(condition, 'No treatment')) %>% 
    mutate(value = log10(value))


qpcr_plot = ggboxplot(qpcr_data, x = "condition", y = "value")+
    geom_jitter(height = 0, width = 0.2, alpha = 0.6) +
    stat_compare_means(aes(label = after_stat(p.signif)),
                       method = 't.test',
                       hide.ns = TRUE,
                       comparisons = list(c("No treatment", "BDNF"),
                                          c("BDNF", "BDNF+CHX")),
                       p.adjust.method = "bonferroni",  
                       size = 10,
                       vjust = 1,
                       tip.length = 0) + 
    facet_wrap(~gene) +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x =  element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab(element_blank()) +
    theme(axis.title.y = element_text(size = 14)) +
    labs(
        y = str_wrap("qPCR levels normalised to GAPDH", width = 20)
    )


hour_two_featurecounts <- readRDS("~/Documents/GitHub/lmn_bdnf/data/hour_two_featurecounts.RDS")

rna_seq_plot = counts(hour_two_featurecounts$deseq_obj,normalized = TRUE) %>% 
    as.data.frame() %>% 
    rownames_to_column('gene') %>% 
    left_join(estimate_list_full %>% distinct(gene,symbol)) %>% 
    filter(symbol %in% qpcr_data$gene) %>%
    reshape2::melt(id.vars = c("gene",'symbol')) %>% 
    separate(variable, into = c("cond","rep","time")) %>% 
    mutate(cond = ifelse(cond == 'CONTROL','Control',"BDNF")) %>% 
    mutate(cond = fct_relevel(cond,"Control")) %>% 
    ggplot(aes(x = cond, y = value,color = cond)) + 
    geom_boxplot() + 
    geom_jitter(size = 3,height = 0) + 
    facet_wrap(~symbol) +
    scale_color_manual(values = c("BDNF" = "#0E7FFE", "Control" = "#FE02FF")) + 
    xlab(element_blank()) +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x =  element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab(element_blank()) +
    theme(axis.title.y = element_text(size = 14)) +
    labs(
        y = str_wrap("Normalised total RNA abundance", width = 20)
    )



other_bois = (rna_seq_plot + qpcr_plot)
other_bois

# FOS ---------------------------------------------------------------------

fos_qpcr = fread("data/fos_qpcr_results.csv")

fos_q_plot = fos_qpcr %>% 
    mutate(value = log10(value)) %>% 
    mutate(treatment = gsub("_"," ",treatment)) %>% 
    mutate(treatment = gsub("BDNF ActD","BDNF+ActD",treatment)) %>% 
    mutate(treatment = fct_relevel(treatment, 'No treatment','BDNF')) %>% 
    ggplot(aes(x = treatment, y = value)) + 
    geom_boxplot() + 
    geom_jitter(size = 3) + 
    facet_wrap(~gene) +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x =  element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab(element_blank()) +
    theme(axis.title.y = element_text(size = 14)) +
    labs(
        y = str_wrap("qPCR levels normalised to GAPDH", width = 20)
    ) +
    stat_compare_means(aes(label = after_stat(p.signif)),
                       method = 't.test',
                       hide.ns = TRUE,
                       comparisons = list(c("No treatment", "BDNF"),
                                          c("BDNF", "BDNF+ActD")),
                       p.adjust.method = "bonferroni",  
                       
                       size = 10,
                       vjust = 1,
                       tip.length = 0) 
    

fos_new_plot = estimate_list_full %>% 
    filter(symbol == 'FOS' & time == 1) %>% 
    mutate(condition = ifelse(condition == 'bdnf','BDNF',"Control")) %>% 
    mutate(condition = fct_relevel(condition,'Control')) %>% 
    ggplot(aes(x = condition, color = condition, y = map)) + 
    geom_boxplot() + 
    geom_jitter(size = 3,height = 0) + 
    scale_color_manual(values = c("BDNF" = "#0E7FFE", "Control" = "#FE02FF")) + 
    ylab("New-to-total RNA fraction") +
    xlab(element_blank()) +
    theme_bw() +
    facet_wrap(~symbol) +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x = element_text(size = 14)) +
    theme(axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::percent_format())

fossy = (fos_q_plot + fos_new_plot)

fossy
