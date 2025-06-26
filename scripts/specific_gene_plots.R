library(ggplot2)
library(tidyr)
library(magrittr)
library(dplyr)
library(data.table)
library(msigdbr)
library(reshape2)
library(gridExtra)  # for arranging plots

new_ratio_bayesian_p_de <- fread("data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(total_rna_sig = case_when(padj < 0.1 & log2FoldChange > 0.75 ~ " upregulated",
                                     padj < 0.1 & log2FoldChange < -0.75 ~ "downregulated",
                                     T ~ "not_significant"))

new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(log2Fold_newRNA = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    mutate(gene = gsub("\\..*", "", gene))



# Filter the base data
base_data <- new_ratio_bayesian_p_de %>% 
    filter(gene_name %in% c("ATF3","JUN","PTEN","STAT3")) %>% 
    select(gene_name, new_rna_sig, total_rna_sig, log2Fold_newRNA, log2FoldChange, time)

# Create New RNA plot
new_rna_plot <- base_data %>%
    mutate(
        significant = ifelse(new_rna_sig != "unclear", "*", ""),
        star_y = ifelse(significant == "*", log2Fold_newRNA + 0.1 * sign(log2FoldChange) + 0.05, NA)
    ) %>% 
    ggplot(aes(x = as.factor(time), y = log2Fold_newRNA)) + 
    geom_col(fill = "#ff7f0e") +  # Orange for new RNA
    geom_text(aes(y = star_y, label = significant), 
              size = 8, fontface = "bold", na.rm = TRUE) +
    facet_wrap(~gene_name, nrow = 1) +
    labs(
        title = "New RNA",
        x = "",
        y = "Log2 Fold Change"
    ) +
    theme_minimal() +
    theme(
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

# Create Total RNA plot  
total_rna_plot = base_data %>%
    mutate(
        significant = ifelse(total_rna_sig != "not_significant", "*", ""),
        star_y = ifelse(significant == "*", log2FoldChange + 0.1, NA)
    ) %>% 
    ggplot(aes(x = as.factor(time), y = log2FoldChange)) + 
    geom_col(fill = "#1f77b4") +  # Blue for total RNA
    geom_text(aes(y = star_y, label = significant), 
              size = 8, fontface = "bold", na.rm = TRUE) +
    facet_wrap(~gene_name,nrow = 1) +
    labs(
        title = "Total RNA",
        x = "Time (h)",
        y = "Log2 Fold Change"
    ) +
    theme_minimal() +
    theme(
        strip.text = element_blank(),  # Remove strip text since it's the same as above
    ) +
    scale_y_continuous(breaks=c(-0.5,0,0.5, 1, 1.5))

# Combine the plots
grid.arrange(new_rna_plot, total_rna_plot, ncol = 1, heights = c(1, 1))
