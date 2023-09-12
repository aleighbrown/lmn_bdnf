library(data.table)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(clusterProfiler)
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(tidyverse)
library(enrichplot)
my_get_wordcloud <- function(cluster, ggData, nWords){

    words <- ggData %>%
        filter(!is.na(label)) %>%
        pull(label) %>%
        gsub(" in ", " ", .) %>%
        gsub(" [0-9]+ ", " ", .) %>%
        gsub("^[0-9]+ ", "", .) %>%
        gsub(" [0-9]+$", "", .) %>%
        gsub(" [A-Za-z] ", " ", .) %>%
        gsub("^[A-Za-z] ", "", .) %>%
        gsub(" [A-Za-z]$", "", .) %>%
        gsub(" / ", " ", .) %>%
        gsub(" and ", " ", .) %>%
        gsub(" of ", " ", .) %>%
        gsub("to", " ", .) %>%
        gsub(",", " ", .) %>%
        gsub(" - ", " ", .)
    
    net_tot <- length(words)
    
    clusters <- unique(ggData$group)

    words_i = ggData |> filter(group == cluster) |> 
        filter(!is.na(label)) |> 
        pull(label) 
    words_i = words_i%>%
        gsub(" in ", " ", .) %>%
        gsub(" [0-9]+ ", " ", .) %>%
        gsub("^[0-9]+ ", "", .) %>%
        gsub(" [0-9]+$", "", .) %>%
        gsub(" [A-Za-z] ", " ", .) %>%
        gsub("^[A-Za-z] ", "", .) %>%
        gsub(" [A-Za-z]$", "", .) %>%
        gsub(" / ", " ", .) %>%
        gsub(" and ", " ", .) %>%
        gsub(" of ", " ", .) %>%
        gsub(",", " ", .) %>%
        gsub(" - ", " ", .)  %>%
        gsub("to", " ", .)
        
    
    sel_tot <- length(words_i)
    sel_w <- enrichplot:::get_word_freq(words_i)
    net_w_all <- enrichplot:::get_word_freq(words)
    net_w <- net_w_all[names(sel_w)]
    tag_size <- (sel_w/sel_tot)/(net_w/net_tot)
    tag_size <- tag_size[order(tag_size, decreasing = TRUE)]
    nWords <- min(nWords, length(tag_size))
    tag <- names(tag_size[seq_len(nWords)])
    
    # Order of words
    dada <- strsplit(words_i, " ")
    len <- vapply(dada, length, FUN.VALUE=1)
    rank <- NULL

    for(i in seq_len(sel_tot)) {

        rank <- c(rank, seq_len(len[i]))
    }

    word_data <- data.frame(word = unlist(dada), rank = rank)
    word_rank1 <- stats::aggregate(rank ~ word, data = word_data, sum)
    rownames(word_rank1) <- word_rank1[, 1]
    
    word_rank1 <- word_rank1[names(sel_w), ]
    # Get an average ranking order
    word_rank1[, 2] <- word_rank1[, 2]/as.numeric(sel_w)
    tag_order <- word_rank1[tag, ]
    tag_order <- tag_order[order(tag_order[, 2]), ]
    tag_clu_i <- paste(tag_order$word, collapse=" ")
    return(tag_clu_i)
}

# Writing out the table to be clustered  ----------------------------------


new_ratio_bayesian_p_de <- fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/data/new_ratio_bayesian_p_de.csv")
new_ratio_bayesian_p_de = new_ratio_bayesian_p_de |> 
    mutate(total_rna_sig = case_when(padj < 0.1 & log2FoldChange > 0.75 ~ "upregulated",
                                     padj < 0.1 & log2FoldChange < -0.75 ~ "downregulated",
                                     T ~ "not_significant"))

never_sig = new_ratio_bayesian_p_de |> 
    filter(total_rna_sig == "not_significant") |> 
    group_by(gene) |> 
    dplyr::count() |> 
    filter(n == 3) |> pull(gene)

log2clust = new_ratio_bayesian_p_de |> 
    filter(baseMean > 15) |> 
    filter(!(gene %in% never_sig)) |>
    # mutate(log2FoldNTR = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    # filter(n_samp_passing_bdnf >2 & n_samp_passing_control > 2) |> 
    dplyr::select(time,gene,gene_name,log2FoldChange) |> 
    pivot_wider(names_from = 'time',
                values_from = c('log2FoldChange'), 
                values_fill = NA) |> 
    na.omit() |> 
    # mutate(
    #        across(where(is.numeric), ~ as.numeric(scale(.)))
    # ) |>
    janitor::clean_names()


# just normal ole clustering plotting ---------------------------------------------------------------


set.seed(42)
df = log2clust |>  
    mutate(gene = glue::glue('{gene}_{gene_name}')) |> 
    unique() |> 
    select(-gene_name) |> 
    tibble::column_to_rownames('gene')


# Dissimilarity matrix
d <- dist(df, method = "euclidean")


# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
# ac <- function(x) {
#     agnes(df, method = x)$ac
# }
# 
# map_dbl(m, ac)
# map_dbl(m, ac)
# average    single  complete      ward 
# 0.9891487 0.9619502 0.9917806 0.9986569

hc5 <- hclust(d, method = "ward.D2" )

sub_grp <- cutree(hc5, k = 16)
avg_cluster = log2clust |> 
    mutate(cluster = sub_grp) %>%
    melt(id.vars = c("gene","gene_name","cluster")) |>
    mutate(time = readr::parse_number(as.character(variable))) |>
    janitor::clean_names() |>
    group_by(cluster,time) |>
    summarize(plot_value_median = median(value,na.rm = TRUE),
              plot_value_mean = mean(value,na.rm = TRUE),
              n = n(),
              sd = sd(value), 
              se = 2*sd) |> 
    mutate(plot_name = glue::glue("cluster {cluster}: {n} genes"))



long_cluster = log2clust |> 
    mutate(cluster = sub_grp) %>%
    # filter(gene_name == 'ARC')
    # filter(probability > 0.6) |> 
    melt(id.vars = c("gene","gene_name","cluster")) |> 
    mutate(time = parse_number(as.character(variable))) |> 
    add_count(cluster) |> 
    mutate(plot_name = glue::glue("cluster {cluster}: {n/3} genes")) |> 
    as.data.table()

# filter(cluster == cluster_n)



ggplot() + 
    # geom_point(aes(x = time,
    #                y = value,
    #                group = gene_name),
    #            alpha = 0.3,data = long_cluster) +
    # ggtitle(glue::glue("Cluster {cluster_n}")) + 
    ggpubr::theme_pubr() +
    geom_point(aes(x = time, y = plot_value_mean,
                   group = cluster),data = avg_cluster) + 
    geom_ribbon(aes(x = time, ymin = plot_value_mean - se,
                    ymax = plot_value_mean + se,
                  group = cluster),fill = '#C1EFFF',
                color = '#6689FF',data = avg_cluster,size = 1.2) + 

    geom_hline(yintercept = 0) + 
    theme(legend.position = 'none') + 
    geom_line(aes(x = time,
                  y = value,
                  group = gene_name),alpha = 0.3,data = long_cluster,color = 'red') +
    geom_line(aes(x = time, y = plot_value_mean,
                  group = cluster),color = 'blue',
              data = avg_cluster) + 
    geom_hline(yintercept = 0) +
    labs(y = bquote('Log'[2]~ 'fold change'),
         x = 'BDNF treatment (hr)') + 
    scale_x_continuous(breaks = c(0,1,2,6)) +
    expand_limits(x = 0, y = 0) +
    facet_wrap(~cluster,scales = 'free') 
     # geom_line(aes(x = time,
    #               y = value,
    #               group = gene_name),data = long_cluster[gene_name == 'ARC'], color = 'black') +
    # geom_line(aes(x = time,
    #               y = value,
    #               group = gene_name),data = long_cluster[gene_name == 'FOS'],color = 'black') +
    # geom_line(aes(x = time,
    #               y = value,
    #               group = gene_name),data = long_cluster[gene_name == 'JUN'],color = 'black') 
    # 



# plotting al the go of the clusters --------------------------------------


long_cluster = long_cluster |> 
    mutate(ens = gsub("\\..*", "",gene))
uni = new_ratio_bayesian_p_de |>  mutate(ens = gsub("\\..*", "",gene)) |> pull(ens) |> unique()
go_objects = data.table()
cluster_with_go = c()

pdf("~/Desktop/cluster_plots_all_ont_expressed_genes_universe.pdf")

for(i in unique(long_cluster$cluster)){
    rm(go_this_cluster)
    rm(go_this_cluster_flat)
    cluster_n = i
    
    gene_in_cluster = long_cluster |> 
        filter(cluster == cluster_n) |> 
        pull(ens) |> 
        unique()
    
    message("starting go")
    print(i)
    go_this_cluster <- enrichGO(gene          = gene_in_cluster,
                                universe = uni,
                                keyType = "ENSEMBL",
                                OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                ont           = "ALL",
                                pAdjustMethod = "BH",
                                readable      = TRUE)
    
    if(any(go_this_cluster@result$p.adjust < go_this_cluster@pvalueCutoff)){
        p = dotplot(go_this_cluster) + ggtitle(glue::glue("Cluster {cluster_n}"))
        print(p)
        s = pairwise_termsim(go_this_cluster)
        tr <- NULL
        # Wrap the code in a tryCatch block
        tryCatch({
            tr <- treeplot(s)
        }, error = function(e) {
            # Handle the error
            print(paste("An error occurred:", e))
            # Set tr to NULL or some other value to indicate an error
            tr <- 1
        })
        
        # Check if tr is NULL to determine if an error occurred
        if (!is.null(tr)) {
           print(tr)
        }
        c = cnetplot(go_this_cluster) + ggtitle(glue::glue("Cluster {cluster_n}"))
        
        print(c)
        message("Done go")
        go_this_cluster_flat = go_this_cluster@result |> 
            mutate(cluster = i)
        go_objects = rbind(go_objects,go_this_cluster_flat)
        cluster_with_go = c(cluster_with_go,i)
        
    }else{
        message("No GO found for ")
        print(i)
    }

}

go_objects |> fwrite('~/Desktop/go_flat_expressed_background.csv')
beepr::beep(4)
dev.off()


# semantic go specific clusters -------------------------------------------


gene_in_cluster = long_cluster |> filter(cluster == 16) |> pull(ens) |> unique()
clipr::write_clip(gene_in_cluster)
go_this_cluster <- enrichGO(gene          = gene_in_cluster,
                            keyType = "ENSEMBL",
                            OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "BH",
                            universe = uni,
                            readable      = TRUE)

dotplot(go_this_cluster)
s = pairwise_termsim(go_this_cluster) 
tr <- treeplot(s,nCluster = 4)
tr
my_get_wordcloud("cluster_2",tr$data,5)
acti# Plotting the log2NTR and NTR of genes in each cluster -------------------

new_ratio_bayesian_p_de |> 
    filter(new_rna_sig != 'unclear') |> 
    mutate(log2ntr = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    left_join(unique(long_cluster[,.(cluster,gene,gene_name)]),by = c("gene","gene_name")) |> 
    filter(!is.na(cluster)) |>  
    ggplot(aes(x = time,
               y = log2ntr,
               group = gene)) + 
    geom_line(alpha = 0.3) +
    facet_wrap(~cluster) +
    ggpubr::theme_pubr() + 
    geom_hline(yintercept = 0)




new_ratio_bayesian_p_de |> 
    filter(new_rna_sig != 'unclear') |> 
    mutate(log2ntr = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    left_join(unique(long_cluster[,.(cluster,gene,gene_name)]),by = c("gene","gene_name")) |> 
    filter(!is.na(cluster)) |> 
    group_by(time,new_rna_sig,cluster) |> 
    tally() |> 
    filter(cluster %in% c(4,9,7)) |> 
    ggplot(aes(x = as.factor(time),
               fill = n,
               y = new_rna_sig)) + 
    geom_tile() +
    ggpubr::theme_pubr() +
    facet_grid(cols = vars(cluster),
               rows = vars(new_rna_sig))



set.seed(42)
ntr_df = new_ratio_bayesian_p_de |> 
    filter(new_rna_sig != 'unclear') |> 
    mutate(log2ntr = log2(mean_bdnf_ntr / mean_control_ntr)) |> 
    select(gene,gene_name,time,log2ntr) |> 
    pivot_wider(names_from = 'time',
                values_from = 'log2ntr') |> 
    janitor::clean_names() |> 
    na.omit() |> 
    mutate(gene = glue::glue('{gene}_{gene_name}')) |> 
    unique() |> 
    select(-gene_name) |> 
    tibble::column_to_rownames('gene')
    

# Dissimilarity matrix
d2 <- dist(ntr_df, method = "euclidean")
hc_test <- hclust(d2, method = "ward.D2" )
sub_grp <- cutree(hc_test, k = 12)
avg_cluster = ntr_df |> 
    mutate(cluster = sub_grp) %>%
    melt(id.vars = c("gene","gene_name","cluster")) |>
    mutate(time = readr::parse_number(as.character(variable))) |>
    janitor::clean_names() |>
    group_by(cluster,time) |>
    summarize(plot_value_median = median(value,na.rm = TRUE),
              plot_value_mean = mean(value,na.rm = TRUE),
              n = n(),
              sd = sd(value), 
              se = 2*sd) |> 
    mutate(plot_name = glue::glue("cluster {cluster}: {n} genes"))



long_cluster_ntr = ntr_df |> 
    mutate(cluster = sub_grp) %>%
    # filter(gene_name == 'ARC')
    # filter(probability > 0.6) |> 
    tibble::rownames_to_column('id') |> 
    melt(id.vars = c("id","cluster")) |> 
    mutate(time = parse_number(as.character(variable))) |> 
    add_count(cluster) |> 
    mutate(plot_name = glue::glue("cluster {cluster}: {n/3} genes")) |> 
    as.data.table()


ggplot() + 
    # geom_point(aes(x = time,
    #                y = value,
    #                group = gene_name),
    #            alpha = 0.3,data = long_cluster) +
    # ggtitle(glue::glue("Cluster {cluster_n}")) + 
    ggpubr::theme_pubr() +
    # geom_point(aes(x = time, y = plot_value_mean,
    #                group = cluster),data = avg_cluster) + 
    # geom_ribbon(aes(x = time, ymin = plot_value_mean - se,
    #                 ymax = plot_value_mean + se,
    #                 group = cluster),fill = '#C1EFFF',
    #             color = '#6689FF',data = avg_cluster,size = 1.2) + 
    # 
    geom_hline(yintercept = 0) + 
    theme(legend.position = 'none') + 
    geom_line(aes(x = time,
                  y = value,
                  group = id),alpha = 0.3,data = long_cluster_ntr,color = 'red') +
    # geom_line(aes(x = time, y = plot_value_mean,
    #               group = cluster),color = 'blue',
    #           data = avg_cluster) + 
    geom_hline(yintercept = 0) +
    labs(y = bquote('Log'[2]~ 'new RNA fold change'),
         x = 'BDNF treatment (hr)') + 
    scale_x_continuous(breaks = c(0,1,2,6)) +
    expand_limits(x = 0, y = 0) +
    facet_wrap(~cluster) 

# DP_GP_cluster plotting ---------------------------------------------------------------
# This was the code to cluster ---------------------------------------------------------------

# DP_GP_cluster.py \
# -i /Users/annaleigh/Documents/GitHub/bdnf_4su/bdnf_test.txt \
# --clusterings /Users/annaleigh/Desktop/clustering/DP_GP_cluster/bdnf/no_mean_no_scale_clusterings.txt \
# --log_likelihoods /Users/annaleigh/Desktop/clustering/DP_GP_cluster/bdnf/no_mean_no_scale_log_likelihoods.txt \
# --criterion MAP \
# --post_process \
# --cluster_uncertainty_estimate \
# --plot --plot_types png \
# --sim_mat /Users/annaleigh/Desktop/clustering/DP_GP_cluster/bdnf/no_mean_no_scale_posterior_similarity_matrix.txt \
# --output /Users/annaleigh/Documents/GitHub/bdnf_4su/bdnf_opt_no_mean_no_scale




clust_t = fread("/Users/annaleigh/Documents/GitHub/bdnf_4su/bdnf_opt_no_mean_no_scale_optimal_clustering.txt")
# clust_t[probability > 0.6,.N,by = cluster] |> arrange(-N)
background = gsub("\\..*", "",new_ratio_bayesian_p_de$gene)
clust_t[,ens := gsub("\\..*", "",gene)]

certain_clustering = clust_t[probability > 0.6] 
#get an average behavior of the certain genes in each cluster
avg_cluster_movement = log2clust |> 
    left_join(clust_t) |> 
    select(-ens) |> 
    filter(probability > 0.6) |>
    melt(id.vars = c("gene","gene_name","cluster","probability")) |>
    mutate(time = readr::parse_number(as.character(variable))) |>
    janitor::clean_names() |>
    group_by(cluster,time) |>
    summarize(plot_value_median = median(value,na.rm = TRUE),
              plot_value_mean = mean(value,na.rm = TRUE),
              n = n(),
              sd = sd(value), 
              se = sd/sqrt(n))

the_clusters = unique(certain_clustering$cluster)
for(i in 1:length(the_clusters)){
    
}
cluster_n = the_clusters[i]
gene_in_cluster = certain_clustering |> filter(cluster == cluster_n) |> pull(ens)

go_this_cluster <- enrichGO(gene          = gene_in_cluster,
                            keyType = "ENSEMBL",
                            OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                            ont           = "ALL",
                            pAdjustMethod = "BH",
                            readable      = TRUE)

p = dotplot(go_this_cluster) + ggtitle(glue::glue("Cluster {cluster_n}"))
print(p)
cnetplot(go_this_cluster)
beepr::beep()




log2clust |> 
    left_join(certain_clustering) |> 
    select(-ens) |> 
    # filter(gene_name == 'ARC')
    # filter(probability > 0.6) |> 
    melt(id.vars = c("gene","gene_name","cluster","probability")) |>
    filter(cluster %in% c(1,7)) |> 
    # filter(cluster == cluster_n)
    ggplot(aes(x = variable,
               y = value,group = gene_name,color = as.character(cluster))) + 
    geom_point() +
    geom_line() + 
    geom_hline(yintercept = 0) +
    # ggtitle(glue::glue("Cluster {cluster_n}")) + 
    ggpubr::theme_pubr()


# GWENA --------------------------------------------------
hour_one_featurecounts <- readRDS("~/Library/CloudStorage/GoogleDrive-brown.annaleigh@gmail.com/My Drive/_BACKUP2/bdnf_4su/hour_one_featurecounts.RDS")
hour_two_featurecounts <- readRDS("~/Library/CloudStorage/GoogleDrive-brown.annaleigh@gmail.com/My Drive/_BACKUP2/bdnf_4su/hour_two_featurecounts.RDS")
hour_six_featurecounts <- readRDS("~/Library/CloudStorage/GoogleDrive-brown.annaleigh@gmail.com/My Drive/_BACKUP2/bdnf_4su/hour_six_featurecounts.RDS")
n1 = counts(hour_one_featurecounts$deseq_obj,normalized = TRUE) |> as.data.frame() |> tibble::rownames_to_column('geneId')
n2 = counts(hour_two_featurecounts$deseq_obj,normalized = TRUE) |> as.data.frame() |> tibble::rownames_to_column('geneId')
n6 = counts(hour_six_featurecounts$deseq_obj,normalized = TRUE) |> as.data.frame() |> tibble::rownames_to_column('geneId')
# Expression data  with genes as columns and samples as rows.
full_no = n1 |> 
    left_join(n2) |> 
    left_join(n6) |> 
    drop_na()

wide_no = full_no |> 
    mutate(geneId = gsub("\\.*","",geneId)) |> 
    tibble::column_to_rownames('geneId') |> 
    t()

wide_no_filtered <- filter_low_var(wide_no, pct = 0.7, type = "median")
ncol(wide_no_filtered)

net <- build_net(wide_no_filtered, cor_func = "spearman", 
                 n_threads = 2)

# Power selected :
net$metadata$power
#> [1] 8

# Fit of the power law to data ($R^2$) :
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]

modules <- detect_modules(wide_no_filtered, 
                          net$network, 
                          detailled_result = TRUE,
                          merge_threshold = 0.75,pam_respects_dendro = TRUE)

# Number of modules before merging :
length(unique(modules$modules_premerge))
#> [1] 6
# Number of modules after merging: 
length(unique(modules$modules))
#> [1] 3

layout_mod_merge <- plot_modules_merge(
    modules_premerge = modules$modules_premerge, 
    modules_merged = modules$modules)

ggplot2::ggplot(data.frame(modules$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
    ggplot2::ylab("Number of genes") +
    ggplot2::xlab("Module")

enrichment <- bio_enrich(modules$modules)
