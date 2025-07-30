library(DESeq2)
library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)
# Function to perform 1-way ANOVA for each gene
perform_anova <- function(data) {
    genes <- unique(data$gene)
    results <- list()
    
    for(g in genes) {
        cat(paste("\n=== ANOVA for Gene:", g, "===\n"))
        gene_data <- data[data$gene == g, ]
        
        # Perform 1-way ANOVA
        anova_result <- aov(value ~ condition, data = gene_data)
        anova_summary <- summary(anova_result)
        
        # Extract key statistics
        f_statistic <- anova_summary[[1]][["F value"]][1]
        p_value <- anova_summary[[1]][["Pr(>F)"]][1]
        df_between <- anova_summary[[1]][["Df"]][1]
        df_within <- anova_summary[[1]][["Df"]][2]
        
        cat(paste("F-statistic:", round(f_statistic, 4), "\n"))
        cat(paste("p-value:", format(p_value, scientific = TRUE, digits = 4), "\n"))
        cat(paste("Degrees of freedom: between =", df_between, ", within =", df_within, "\n"))
        
        # Print full ANOVA table
        print(anova_summary)
        
        # Store results
        results[[g]] <- data.frame(
            gene = g,
            f_statistic = f_statistic,
            p_value = p_value,
            df_between = df_between,
            df_within = df_within,
            significant = p_value < 0.05
        )
        
        # Post-hoc analysis if significant
        if(p_value < 0.05) {
            cat("\nSignificant ANOVA! Performing post-hoc tests...\n")
            
            # Tukey's HSD for pairwise comparisons
            tukey_result <- TukeyHSD(anova_result)
            print(tukey_result)
            
            # Store Tukey results
            tukey_df <- as.data.frame(tukey_result$condition)
            tukey_df$comparison <- rownames(tukey_df)
            tukey_df$gene <- g
            results[[paste0(g, "_tukey")]] <- tukey_df
        } else {
            cat("\nNon-significant ANOVA. No post-hoc tests needed.\n")
        }
        
        # Check ANOVA assumptions
        cat("\n--- Checking ANOVA Assumptions ---\n")
        
        # 1. Normality of residuals (Shapiro-Wilk test)
        residuals <- residuals(anova_result)
        shapiro_test <- shapiro.test(residuals)
        cat(paste("Shapiro-Wilk test for normality: p =", 
                  format(shapiro_test$p.value, scientific = TRUE, digits = 4), "\n"))
        cat("(p > 0.05 suggests residuals are normally distributed)\n")
        
        # 2. Homogeneity of variances (Levene's test)
        library(car)
        levene_test <- leveneTest(value ~ condition, data = gene_data)
        cat(paste("Levene's test for equal variances: p =", 
                  format(levene_test$`Pr(>F)`[1], scientific = TRUE, digits = 4), "\n"))
        cat("(p > 0.05 suggests equal variances)\n")
        
        cat("\n" %+% paste(rep("-", 50), collapse = "") %+% "\n")
    }
    
    return(results)
}


#FOS total and new, ACTD - 1 hour treatment

#DCTN2/1I2 total and qPCR translation

qpcr_data = fread("data/qpcr_results_bdnf.csv")

estimate_list_full = fread('data/estimate_list_full.csv')

plot_data = qpcr_data %>% 
    mutate(condition = gsub("\\+","\n+\n",condition)) %>% 
    mutate(condition = fct_relevel(condition, 'No treatment'))

# Run ANOVA analysis
anova_results <- perform_anova(plot_data)

# Create summary table of main results
main_results <- do.call(rbind, anova_results[grep("_tukey", names(anova_results), invert = TRUE)])
rownames(main_results) <- NULL

cat("\n=== SUMMARY OF ANOVA RESULTS ===\n")
print(main_results)
significant_genes <- main_results$gene[main_results$significant]

if(length(significant_genes) > 0) {
    cat(paste("\nCreating post-hoc comparison plot for significant genes:", 
              paste(significant_genes, collapse = ", "), "\n"))
    
    # Filter data for significant genes only
    sig_data <- plot_data[plot_data$gene %in% significant_genes, ]
    
    if(nrow(sig_data) > 0) {
        posthoc_plot <- sig_data %>% 
            ggplot(aes(y = value, x = condition)) + 
            geom_boxplot() + 
            geom_jitter(width = 0.4) + 
            facet_wrap(~gene) + 
            theme_bw() +
            ylab("qPCR levels normalised to GAPDH") +
            xlab(element_blank()) +
            stat_compare_means(method = "anova", 
                               label.y = max(sig_data$value, na.rm = TRUE) * 1.2) +
            stat_compare_means(comparisons = list(c("No treatment", "BDNF"),
                                                  c("BDNF", "BDNF\n+\nCHX"),
                                                  c("BDNF", "BDNF\n+\nTorin")),
                               method = "t.test",  # or "tukey_hsd"
                               label = "p.signif") +
            theme(text = element_text(size = 18))
        
        

    }
}

posthoc_plot


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
    ylab("Normalised total RNA abundance") +
    xlab(element_blank()) +
    theme_bw() +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 18))



other_bois = (posthoc_plot | rna_seq_plot)


# FOS ---------------------------------------------------------------------

fos_qpcr = fread("/Users/annaleigh/Downloads/fos_qpcr_results.csv")

fos_q_plot = fos_qpcr %>% 
    filter(gene == 'FOS') %>% 
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
    )
    

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
