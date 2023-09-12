library(data.table)
library(purrr)
library(tidyverse)
library(distr)

my_clean_reader <- function(file){
    df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
    return(df)
}
#this function computes the variances of a beta distribution,https://en.wikipedia.org/wiki/Beta_distribution#Variance
beta_var = function(alpha,beta){
    numerator = alpha * beta
    denom = ((alpha + beta)^2)*(alpha + beta + 1) 
    return(numerator / denom)
}

getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

sample_from_beta = function(alphas, betas, n = 1000){
    #find the variances of all the alphas
    beta_variances = mapply(beta_var, alphas, betas)
    beta_variances = beta_variances / sum(beta_variances)
    
    #weight each beta function by it's inverse variance, and then sum it to combine
    x_values=  seq(0,1.0,0.001)
    weighted_beta = function(weight,alpha,beta, x=  x_values){return(weight * dbeta(x, alpha, beta))}
    
    prob = rowSums(mapply(weighted_beta, beta_variances,alphas,betas))
    prob[is.infinite(prob)] <- 0
    #the sum of the prob should be approxiamately one for this to work
    if(sum(prob) != 1){
        prob = prob / (sum(prob))
    }
    #using the distr packagine
    D<-DiscreteDistribution(supp=x_values,prob= prob)
    rdist <- r(D)                 # function to create random variates from p
    sampled <- rdist(n)                 # sample from X ~ p
    return(sampled)
}

#this returns the probability that the beta distribution 1 is less than that of 2

bayesian_p_value = function(df_table, n = 1 * 10 ^ 5, return_central_values = TRUE){
    #read conditions
    conds = df_table %>% pull(condition) %>% unique()
    
    s1 = sample_from_beta(df_table %>% filter(condition == conds[1]) %>% pull(alpha), 
                          df_table %>% filter(condition == conds[1]) %>% pull(beta), n = n)
    
    s2 = sample_from_beta(df_table %>% filter(condition == conds[2]) %>% pull(alpha), 
                          df_table %>% filter(condition == conds[2]) %>% pull(beta), n = n)
    
    if(return_central_values){
        values = list((sum(sample(s1) < sample(s2))/ n),
                      (mean(sample(s2) - sample(s1))),
                      median(sample(s2) - sample(s1)), 
                      mean(log2(s2/s1)), 
                      median(log2(s2/s1)),
                      median(s1),
                      median(s2),
                      mean(s1),
                      mean(s2), 
                      getmode(s1),
                      getmode(s2),
                      conds[1],
                      conds[2])
        names(values) = c("bayesian_p","mean_diff"
                          ,"median_diff","mean_log2Fold",
                          "median_log2Fold", 
                          glue::glue("median_{conds[1]}_ntr"),
                          glue::glue("median_{conds[2]}_ntr"),
                          glue::glue("mean_{conds[1]}_ntr"),
                          glue::glue("mean_{conds[2]}_ntr"),
                          glue::glue("mode_{conds[1]}_ntr"),
                          glue::glue("mode_{conds[2]}_ntr"),
                          "s1_names",
                          "s2_names")
        values = as_tibble(values)
        return(values)
    }
    
    
    return(sum(sample(s1) < sample(s2))/ n)
}

make_comparison_plot = function(estimate_list_full, t = 2, gene_name = "EGR1"){
    
    b_alpha = estimate_list_full %>% 
        filter(time == t, 
               symbol == gene_name,
               condition == "bdnf") %>% 
        pull(alpha)
    
    b_beta = estimate_list_full %>% 
        filter(time == t, 
               symbol == gene_name,
               condition == "bdnf") %>% 
        pull(beta)
    
    ct_alpha = estimate_list_full %>% 
        filter(time == t, 
               symbol == gene_name,
               condition == "control") %>% 
        pull(alpha)
    
    ct_beta = estimate_list_full %>% 
        filter(time == t, 
               symbol == gene_name,
               condition == "control") %>% 
        pull(beta)
    
    bdnf = sample_from_beta(b_alpha,b_beta)
    
    ctrl = sample_from_beta(ct_alpha,ct_beta)
    
    my_plot = tibble(bdnf,ctrl) %>% 
        melt() %>% 
        ggplot() + 
        geom_density(aes(x = value, color = variable)) +
        ggtitle(glue::glue("{gene_name} - {t} hrs")) + 
        xlab("Estimate New RNA fraction")
    
    print(my_plot)
    return(my_plot)
    
}


all_time_comparison = function(estimate_list_full, gene_name = "EGR1"){
    
    one = make_comparison_plot(estimate_list_full,t = 1, gene_name = gene_name)
    two = make_comparison_plot(estimate_list_full,t = 2, gene_name = gene_name)
    six = make_comparison_plot(estimate_list_full,t = 6, gene_name = gene_name)
    this_plot = ggpubr::ggarrange(one,two,six,common.legend = T)
    return(this_plot)
}

total_rna_time_plot = function(full_de, symbol = "EGR1"){

    he_plot = full_de %>% filter(gene_name == symbol) %>% 
        mutate(sig_change = !(padj > 0.1 | is.na(padj))) %>% 
        ggplot(aes(x = time, y= log2FoldChange,fill = sig_change)) + 
        geom_hline(yintercept = 0,size = 2) + 
        geom_point(size = 3,pch = 21) + 
        ggtitle(glue::glue('{symbol} Total RNA BDNF v Control')) + 
        scale_fill_manual(values = c("black","#E69F00")) + 
        ggpubr::theme_pubr()
    
    return(he_plot)
}

loglik=function(d,a,b,t){
    sum((a-1)*log(1-exp(-t*d))-t*d*b)  
} 

mle_decay = function(alphas, betas, times, decays = c(0,5,10,50)) {
    
    decay_rate = optimize(loglik,decays,a=alphas,b=betas,t= times,maximum=T)$maximum
    return(decay_rate)
}

# approximate confidence intervals
getCI <- function(half_life_df,  confidence=0.95) {
    #adapted from https://academic.oup.com/bib/article/22/6/bbab219/6315814#312133059
    interval=c(1e-12,2)
    
    alphas = half_life_df$alpha
    betas = half_life_df$beta
    times = half_life_df$time
    
    optimized_decay_rate = unique(half_life_df$decay)
    
    threshold <- stats::qchisq(confidence, 1)/2
    
    objective <- function(x) loglik(optimized_decay_rate, a=alphas,b=betas,t= times) - loglik(x, a=alphas,b=betas,t= times) - threshold
    
    optimalObjective <- objective(optimized_decay_rate)
    
    ci <- c(NA, NA)
    
    if (optimalObjective * objective(interval[1]) < 0)
        ci[1] <- stats::uniroot(objective, c(interval[1], optimized_decay_rate))$root
    if (optimalObjective * objective(interval[2]) < 0)
        ci[2] <- stats::uniroot(objective, c(optimized_decay_rate, interval[2]))$root
    
    names(ci) = c("lower","upper")
    return(ci)
}
estimate_transcription_rate_off_new = function(new_rna,time,decay){
    transcript_rate = (-new_rna * decay) / ((exp(1)^-time * decay) - 1)
    return(transcript_rate)
}

estimate_transcription_rate_off_old = function(old_rna,time,decay){
    transcript_rate = (old_rna * decay)/((exp(1)^(-time * decay)))
    return(transcript_rate)
}
my_clean_reader <- function(file){
    df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
    return(df)
}

fix_column_names = function(df,suffix = "h.tsv"){
    #first we'll get the value in the SampleID column,this will be unique
    sample_id = unique(df$SampleID)
    names(df) = gsub("_","",gsub(sample_id,"",names(df)))
    vals = purrr::simplify(strsplit(sample_id,"_"))
    df$condition = vals[1]
    df$well = vals[2]
    
    df$time = vals[3]
    return(df)
}


draw_beta = function(alpha, beta, vector_names,xlim = c(0,1)){
    p = seq(0,1, length=10000)
    prob = dbeta(p, alpha ,  beta)
    prob = map2(alpha, beta, ~dbeta(x = p, shape1 = .x, shape2 = .y,))
    names(prob) = vector_names
    
    df = data.frame(p,prob) %>% melt(id.vars = 'p')
    
    plt = ggplot(df) + 
        geom_line(aes(x = p, y = value, color = variable),size = 1.5) +
        xlim(xlim)
    
    print(plt)
    return(plt)
}

draw_gene_beta = function(estimate_list_full = estimate_list_full, gene_name = "RN7SK", xlim = c(0,1)){
    df = estimate_list_full[symbol == gene_name]
    df = df %>% mutate(vector_name = glue::glue("{condition}_{time}"))
    
    alphas = df$alpha
    betas = df$beta
    draw_beta(alphas, betas, df$vector_name,xlim)
    
}

