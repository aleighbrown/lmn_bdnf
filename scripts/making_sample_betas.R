# Load required libraries
library(distr)
library(ggplot2)
#this function computes the variances of a beta distribution,https://en.wikipedia.org/wiki/Beta_distribution#Variance
beta_var = function(alpha,beta){
    numerator = alpha * beta
    denom = ((alpha + beta)^2)*(alpha + beta + 1) 
    return(numerator / denom)
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

make_df = function(a1,b1,a2,b2,group = 'a'){
    # Parameters for the first beta distribution
    alpha1 <- a1
    beta1 <- b1
    
    # Parameters for the second beta distribution
    alpha2 <- a2
    beta2 <- b2
    
    # Generate data points from the beta distributions
    x <- seq(0, 1, by = 0.01)
    y1 <- dbeta(x, shape1 = alpha1, shape2 = beta1)
    y2 <- dbeta(x, shape1 = alpha2, shape2 = beta2)
    
    # Create a data frame for plotting
    df <- data.frame(x = rep(x, 2), y = c(y1, y2), 
                     distribution = rep(c("Control", "BDNF"), each = length(x)))
    df = df |> mutate(group = group)
    return(df)
}



plot_beta_distributions <- function(data) {
    
    # Create a sequence of x values
    x_seq <- seq(0, 1, by = 0.01)
    
    
    combined_df_bdnf = sample_from_beta(data[condition == 'bdnf',alpha],data[condition == 'bdnf',beta])
    combined_df_control = sample_from_beta(data[condition == 'control',alpha],data[condition == 'control',beta])
    
    
    df_to_plot <- expand.grid(x = x_seq, SampleID = data$SampleID)
    df_to_plot <- df_to_plot %>%
        left_join(data, by = "SampleID") %>%
        rowwise() %>%
        mutate(y = dbeta(x, alpha, beta))
    
    p <- ggplot(df_to_plot, aes(x = x, y = y)) +
        scale_color_manual(values = c("bdnf" = "#0E7FFE", "control" = "#FE02FF")) +
        theme_minimal() +
        labs(y = "Density", x = "Value", title = "Beta Distributions")
    
    # Dynamically add geom_line calls for each unique SampleID
    for(sample_id in unique(df_to_plot$SampleID)) {
        p <- p + geom_line(data = subset(df_to_plot, SampleID == sample_id), aes(color = condition))
    }
    
    p = p + geom_density(inherit.aes = FALSE, data = tibble(combined_df_bdnf),aes(x = combined_df_bdnf),color = "#0E7FFE",size = 2.5)
    p = p + geom_density(inherit.aes = FALSE, data = tibble(combined_df_control),aes(x = combined_df_control),color = "#FE02FF",size = 2.5)
    
    
    return(p)
}


estimate_list_full = fread("estimate_list_full.csv")
draw_plot = function(t,s){
    df = estimate_list_full |> filter(time == t & symbol == s) |> dplyr::select(SampleID,time,condition,alpha,beta)
    plot_beta_distributions(df)
}

draw_plot(6,'COL4A2')
