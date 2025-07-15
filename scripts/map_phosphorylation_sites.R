# Load required libraries
library(dplyr)
library(stringr)


# Function to find positions in string1 corresponding to characters that precede '(' in string2
find_positions_before_parentheses <- function(str1, str2) {
    # Find all positions of '(' in string2
    paren_positions <- gregexpr("\\(", str2)[[1]]
    
    if (paren_positions[1] == -1) {
        return(numeric(0))  # No parentheses found
    }
    
    
    positions_in_str1 <- numeric()
    
    # For each parenthesis, we need to map back to string1
    for (paren_pos in paren_positions) {
        # Count how many characters from string1 appear before this parenthesis in string2
        
        # Extract the part of string2 before this parenthesis
        before_paren <- substr(str2, 1, paren_pos - 1)
        
        # Remove all previous parenthetical insertions from this substring
        clean_before <- before_paren
        while (grepl("\\([^)]*\\)", clean_before)) {
            clean_before <- gsub("\\([^)]*\\)", "", clean_before)
        }
        
        
        # The length of this cleaned string is the position in string1
        pos_in_str1 <- nchar(clean_before)
        positions_in_str1 <- c(positions_in_str1, pos_in_str1)
        
        if (pos_in_str1 > 0 && pos_in_str1 <= nchar(str1)) {
        }
    }
    
    return(positions_in_str1)
}

map_phospho_site <- function(phospho_prob, full_seq) {
    
    # Step 1: Extract probabilities and find sites with highest probability
    # Extract all patterns like "LETTER(NUMBER)"
    matches <- str_extract_all(phospho_prob, "[STY]\\([0-9.]+\\)")[[1]]
    
    if (length(matches) == 0) {
        return(list(list(position = NA, site = NA, probability = NA, peptide = NA, error = "No matches found")))
    }
    
    # Extract probabilities and letters
    probs <- as.numeric(str_extract(matches, "(?<=\\()[0-9.]+(?=\\))"))
    clean_peptide <- str_remove_all(phospho_prob, "\\([0-9.]+\\)")
    
    # extract the positions within the string that could have the site
    all_possible_sites = find_positions_before_parentheses(clean_peptide,phospho_prob)
    
    # Find ALL sites with the highest probability (handles ties)
    max_prob <- max(probs)
    max_prob_indices <- which(probs == max_prob)
    these_sites_in_peptide = all_possible_sites[max_prob_indices]
    # find which of the possible have that max prob
    
    # Step 2: Clean the peptide sequence (remove numbers and parentheses)
    
    # Step 3: Find where this peptide maps in the full sequence
    peptide_start <- str_locate(full_seq, clean_peptide)[1]
    
    if (is.na(peptide_start)) {
        return(list(list(position = NA, site = NA, probability = max_prob, 
                         peptide = clean_peptide, error = "Peptide not found in sequence")))
    }
    
    # Step 5: Create results for all sites with maximum probability
    results <- list()
    
    for (i in 1:length(these_sites_in_peptide)) {

        site_pos_in_peptide = these_sites_in_peptide[i]
        # Calculate absolute position in full sequence
        absolute_position <- peptide_start + site_pos_in_peptide - 1
        
        results[[i]] <- list(
            position = absolute_position,
            probability = max_prob,
            peptide = clean_peptide,
            peptide_start = peptide_start,
            site_pos_in_peptide = site_pos_in_peptide
        )
    } 
    return(results)
}


# Main function to map phosphorylation sites
map_phosphorylation_sites <- function(df) {
    
    # Apply the function to the dataframe and expand rows for multiple sites
    expanded_data <- list()
    
    for (i in 1:nrow(df)) {
        results <- map_phospho_site(df$phospho_sty_probabilities[i], df$seq[i])
        
        # For each result, create a row with all original data plus phospho info
        for (j in 1:length(results)) {
            row_data <- df[i, ]  # Copy original row
            
            # Add phosphorylation mapping results
            row_data$phospho_position <- ifelse(is.null(results[[j]]$position), NA, results[[j]]$position)
            row_data$phospho_site <- ifelse(is.null(results[[j]]$site), NA, results[[j]]$site)
            row_data$phospho_probability <- ifelse(is.null(results[[j]]$probability), NA, results[[j]]$probability)
            row_data$clean_peptide <- ifelse(is.null(results[[j]]$peptide), NA, results[[j]]$peptide)
            row_data$site_number <- j  # Add indicator for which site this is
            
            expanded_data[[length(expanded_data) + 1]] <- row_data
        }
    }
    
    # Combine all rows into final dataframe
    df_expanded <- do.call(rbind, expanded_data)
    
    return(df_expanded)
}


