# Title: Analysis of comparmentalised immune responses 
# Version: 1.0
# Date: 2025-07-24
# Author: Dr Alexander J Keeley
# Inputs: 
#   - data/OF_blood_correlation_df.RDS
#   - data/compartmentalisation_ratio_df.RDS
#   - scripts/setup_environment.R
#   - data/df_for_pca.RDS
# Outputs:
#   - Correlation matrix visualisation
#   - Dunn’s test comparisons across age groups
#   - Boxplots of compartmentalisation scores by antigen and age
#   - PCA analysis and plots

# Description:
# This script performs an exploratory baseline analysis of immune compartmentalisation in participants with no evidence of Strep A disease events. It includes:
# 1. Calculating pairwise correlation coefficients across blood IgG and oral fluid IgA levels at baseline.
# 2. Reordering and visualising the correlation matrix by antigen and sample type.
# 3. Extracting correlation coefficient ranges within each sample type group.
# 4. Calculating a compartmentalisation ratio (oral IgA : blood IgG) and evaluating its association with age.
# 5. Performing Dunn’s post hoc tests to detect age-group differences and visualising the findings with boxplots and significance brackets.




# Load librarian
library(librarian)

# Load required packages using librarian
shelf(tidyverse, mfp, CorrMixed, gridExtra,corrplot) 

# Setup analysis environment
source("scripts/setup_environment.R")


#################################################
#### Baseline correlation coefficients ####
################################################


# Load baseline datafram: IgG and IgA levels at study baseline provided no Strep A disease event occurred at the time
correlation_df <- readRDS("data/OF_blood_correlation_df.RDS")

# Full correlation matrix across all variables
cor_matrix <- cor(correlation_df %>% ungroup() %>% select(-pid), use = "pairwise.complete.obs")

# View the correlation matrix
print(cor_matrix)

cor_matrix

# Reorder columns logically by group and antigen
col_order <- c(
    "DNAseB Blood IgG", "GAC Blood IgG", "SLO Blood IgG", "SpyAD Blood IgG", "SpyCEP Blood IgG",
    "GAC OF IgA", "SLO OF IgA", "SpyAD OF IgA", "SpyCEP OF IgA"
)

cor_matrix <- cor_matrix[col_order, col_order]

# Plot the reordered correlation matrix
corrplot::corrplot.mixed(cor_matrix, upper = "color", lower = "number",
                         tl.col = "black", tl.srt = 45, lower.col = "black", number.cex = .8,
                         tl.pos = "lt")  # "lt" for left and top labels

# Convert matrix to long-form dataframe
cor_data_long <- as.data.frame(as.table(cor_matrix))

# Rename columns for clarity
colnames(cor_data_long) <- c("row", "col", "coef")

# Extract group and antigen information
cor_data_long <- cor_data_long %>%
    separate(row, into = c("Antigen1", "Group1"), sep = " ", extra = "merge") %>%
    separate(col, into = c("Antigen2", "Group2"), sep = " ", extra = "merge")

# View the resulting dataframe
print(cor_data_long)

# Remove self-correlations for comparison
cor_data_filtered <- cor_data_long %>%
    filter(Antigen1 != Antigen2)

# Calculate the range of correlation coefficients for each group
range_df <- cor_data_filtered %>%
    filter(Group1 == Group2) %>%
    group_by(Group1) %>%
    summarise(min_coef = min(coef, na.rm = TRUE), max_coef = max(coef, na.rm = TRUE))

# Extract the ranges for each group
blood_igg_range <- range_df %>% filter(Group1 == "Blood IgG")
of_iga_range <- range_df %>% filter(Group1 == "OF IgA")


# Construct sentence for the manuscript 
sentence <- sprintf(
    "When comparing antibody titres within individuals from the same group, a strong correlation was observed between antibodies to all antigens. Specifically, the correlation coefficients ranged from %.2f to %.2f for Blood IgG titres, and %.2f to %.2f for OF IgA titres.",
    blood_igg_range$min_coef, blood_igg_range$max_coef,
    of_iga_range$min_coef, of_iga_range$max_coef
)

# Print the sentence
print(sentence)


###############################################
#### Compartmentalisation scores by age #######
###############################################


# Load data containing Z-score transformed IgA and IgG levels at baseline
# Compartmentalisation score = ratio of OF IgA to Blood IgG Z scores

compartmentalisation_ratio_df <- readRDS("data/compartmentalisation_ratio_df.RDS")


# Perform Dunn’s test to assess differences in log-transformed compartmentalisation scores by age group

dunn_results <- compartmentalisation_ratio_df %>%
    mutate(log_score = log(Compartmentalization_Score)) %>%
    group_by(Antigen) %>%
    rstatix::dunn_test(log_score ~ age_grp, p.adjust.method = "bonferroni") %>%
    filter(p.adj < 0.05) 

# Stack y.positions so brackets don't overlap
dunn_results <- dunn_results %>%
    group_by(Antigen) %>%
    mutate(
        y.position = (5)
        + 0.7 * row_number()
    )

# Plot compartmentalisation score by age group. 

plot_IgA_age_group <- ggplot(compartmentalisation_ratio_df, aes(x = age_grp, y = log(Compartmentalization_Score))) +
    geom_boxplot(outliers = F)+
    geom_jitter(aes(color = Antigen), alpha = 0.5, width = 0.2) +
    facet_wrap(~Antigen, nrow = 1) +
    theme_minimal() +
    ggpubr::stat_pvalue_manual(
        dunn_results,
        label = "p.adj.signif",
        hide.ns = TRUE
    ) +
    labs(
        #  title = "Compartmentalization Score (Ratio of Z score transformed antibody levels) vs. Age",
        x = "Age",
        y = "log OF IgA:Blood IgG Z score ratio"
    ) +
    guides(col = "none") +
    scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
    theme_universal() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 1))

# Display the plot
plot_IgA_age_group


###############################################
#### Compartmentalisation PCA analysis ########
###############################################


# Inputs:
#   
# Outputs:
#   - PCA plots of individuals and variable loadings
#   - Combined figure panels (scores and loadings by group)

# Description:
# This script performs principal component analysis (PCA) on all paired oral fluid IgA and blood IgG titres across antigens. It:
# 1. Performs PCA on the compartmentalised immune response matrix.
# 2. Visualises PCA scores stratified by age group (individual projections).
# 3. Visualises PCA loadings coloured by sample type (oral fluid vs blood).
# 4. Combines the score and loading plots for figure inclusion in the manuscript.

shelf(mixOmics)# PCA with enhanced functionality



# Load PCA input dataframe

df_for_pca <- readRDS("data/df_for_pca.RDS")


# Perform multivariate PCA using mixOmics

pca_result <- pca(X = df_for_pca %>%
                        select(-age_grp),
                      scale = T, 
                        ncomp = 5)

# Extract age group metadata
pca_age_group <- df_for_pca %>% select(age_grp)



#### PCA Scores (Individuals) ########


# Extract individual PC coordinates (scores)

scores_df <- as_tibble(pca_result$variates$X)
scores_df$age_grp <- pca_age_group$age_grp  # Add age group

# Plot PC1 vs PC2 coloured by age group
figb <- ggplot(scores_df, aes(x = PC1, y = PC2, color = age_grp)) +
    geom_point(size = 3, alpha = 0.2) +
    labs(x = "PC1", y = "PC2") +
    theme_minimal() +
    stat_ellipse(type = "norm", linetype = 2)


#### PCA Loadings (Variables) ########


# Extract variable loadings
loadings_df <- as_tibble(pca_result$loadings$X)
loadings_df$variable <- rownames(pca_result$loadings$X)

# Extract antigen and group from variable names
loadings_df <- loadings_df %>%       mutate(
    group = word(variable, -1),                     # last word
    Antigen = word(variable, 1, -2, sep = fixed(" ")) # everything but last word
)

shelf(ggrepel)

figc <- 
    loadings_df %>%
    mutate(group = recode(group,
                          "OF_IgA" = "Oral fluid IgA",
                          "Blood_IgG" = "Blood IgG")) %>%
    
    ggplot( aes(x = PC1, y = PC2, col = group, label = Antigen)) +
    geom_segment(
        aes(xend = 0, yend = 0),
        arrow = arrow(length = unit(0.2, "cm"), ends = "first"),
        alpha = 0.7
    ) +
    geom_text_repel(size = 3, max.overlaps = Inf, col = "black", nudge_x = 0.1) +
    coord_equal() +
    labs(x = "PC1", y = "PC2") +
    theme_minimal()

figc

shelf(cowplot)


#### Combine Figures for Output #######

updated_fig <- plot_grid(
    figb +
        labs(col = "Age group")
    
    
    , figc +
        labs(col = "Compartment"), 
    labels = c("A", "B"))

updated_fig





