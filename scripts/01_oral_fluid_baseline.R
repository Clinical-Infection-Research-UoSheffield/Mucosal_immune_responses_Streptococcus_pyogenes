# Title: Baseline Analysis of IgA Levels in Participants Without Strep A Disease Events from SpyCATS household cohort
# Version: 1.0
# Date: 2025-07-21
# Author: Dr Alexander J Keeley
# Inputs: data/
# Outputs: Plots and summary statistics describing baseline antibody titres across age groups and antigens

# Description:

# This script performs a baseline analysis of IgA titres in participants without Strep A disease events. The analysis focuses on:
# 1. Extracting baseline IgA levels from participants with no Strep A positive disease events.
# 2. Fitting fractional polynomial models to predict IgG titres based on age and generating residual and fitted values.
# 3. Calculating centile values and adding prediction bands based on these centiles.
# 4. Generating plots for each antigen showing the distribution of titres and prediction bands.
# 5. Performing correlation matrix and analysis of baseline oral fluid IgA and blood IgG levels 

# Requirements:

# 
library(librarian)
# Load required packages using librarian
shelf(tidyverse, mfp, CorrMixed, gridExtra) 

# Setup environment
source("scripts/setup_environment.R")


########### Read IgA data ##################

# To protect the confidentiality of individual participants, 
# this dataframe has been modified before public release.
# Specifically, we have replaced the actual ages with a "pseudo-age".
# For each participant, a random age is generated within their corresponding age group.
# This approach maintains the overall age group distribution and supports reproducibility of analyses,
# while significantly reducing the risk of re-identifying individual participants.

df <-readRDS("data/baseline_IgA_blood_no_disease_titres.RDS")

# create an empty dataframe to append the output of your for loop
centile_df <- data.frame()


##### update 10/07/2024

# Create an empty list to store plots

plots <- list()

plots2 <- list()
# loop through each antigen

for (a in unique(df$Antigen)) {        
    
    
    # create a datframe for that antigen
    ctrl <- df %>%
        filter(Antigen == a) %>%
        arrange(age)                                
    
    common_y_limits <- c(min(df$titre -1), max(df$titre))
    
    
    # fit a fractional polynomial model 
    mfp_model <- mfp(titre ~ fp(age), data = ctrl, family =  gaussian())  
    # add the residuals to filtered dataframe 
    ctrl$.resid <- residuals(mfp_model)             
    # add the fitted values to the dataframe
    ctrl$.fitted <- fitted(mfp_model)              
    
    # Calculate RMSE
    rmse <- sqrt(mean(ctrl$.resid^2))       
    
    # Calculate centiles
    centile_values <- quantile(ctrl$.resid / rmse, c(0.025, 0.50, 0.80, 0.975))
    
    # Add the upper prediction bands based on centiles and update centile dataframe
    ctrl <- ctrl %>%
        mutate(
            upp25 = .fitted + centile_values[1] * rmse,
            upp50 = .fitted + centile_values[2] * rmse,
            upp80 = .fitted + centile_values[3] * rmse,
            upp975 = .fitted + centile_values[4] * rmse
        )
    
    # Append the centile data for this antigen to the centile dataframe
    centile_df <- rbind(centile_df, ctrl)
    
    # Plot
    p <- ggplot(ctrl, aes(x = age, y = titre, col = Antigen)) +
        scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
        guides(color = "none") +
        geom_point(alpha = 0.5) +
        geom_line(aes(y = .fitted), linetype = "dashed", color = "black") +
        geom_line(aes(y = upp25), linetype = "dashed", color = "red", alpha = 0.5) +
        geom_line(aes(y = upp80), linetype = "dashed", color = "blue") +
        geom_line(aes(y = upp975), linetype = "dashed", color = "red", alpha = 0.5) +
        labs(x = "Age", title = paste0(a), y="Titre (RLU/mL)") +
        theme_minimal() +
        scale_y_continuous(limits =common_y_limits, 
                           #breaks = c(1:6), 
                           labels = log10_to_exp) +
        theme_universal(base_size = plot_basesize)
    
    # Store the plot in the list with the antigen name as the key
    
    plots[[a]] <- p 
    
    
    
}

# Now, arrange the plots into a grid layout

# If you have many plots, adjust the number of rows and columns as needed

plot_grid <- do.call(grid.arrange, c(plots, ncol = 4))

plots2 <- list()
# loop through each antigen

for (a in unique(df$Antigen)) {        
    
    
    # create a datframe for that antigen
    ctrl <- df %>%
        filter(Antigen == a) %>%
        arrange(age)                                
    
    
    common_y_limits <- c(min(df$titre -1), max(df$titre))
    
    
    # fit a fractional polynomial model 
    mfp_model <- mfp(titre ~ fp(age), data = ctrl, family =  gaussian())  
    # add the residuals to filtered dataframe 
    ctrl$.resid <- residuals(mfp_model)             
    # add the fitted values to the dataframe
    ctrl$.fitted <- fitted(mfp_model)              
    
    # Calculate RMSE
    rmse <- sqrt(mean(ctrl$.resid^2))       
    
    # Calculate centiles
    centile_values <- quantile(ctrl$.resid / rmse, c(0.025, 0.50, 0.80, 0.975))
    
    # Add the upper prediction bands based on centiles and update centile dataframe
    ctrl <- ctrl %>%
        mutate(
            upp25 = .fitted + centile_values[1] * rmse,
            upp50 = .fitted + centile_values[2] * rmse,
            upp80 = .fitted + centile_values[3] * rmse,
            upp975 = .fitted + centile_values[4] * rmse
        )
    
    # Append the centile data for this antigen to the centile dataframe
    centile_df <- rbind(centile_df, ctrl)
    
    
    # Plot
    p2 <- 
        
        ctrl %>%
        filter(age <=15) %>%
        ggplot(aes(x = age, y = titre, col = Antigen)) +
        scale_colour_manual(values = c("SpyCEP" = "#FDC086", "SpyAD" = "#d19c2f", "SLO" = "#386CB0", "GAC" = "#7FC97F", "DNAseB" = "#BEAED4")) +
        guides(color = "none") +
        geom_point(alpha = 0.5) +
        geom_line(aes(y = .fitted), linetype = "dashed", color = "black") +
        geom_line(aes(y = upp25), linetype = "dashed", color = "red", alpha = 0.5) +
        geom_line(aes(y = upp80), linetype = "dashed", color = "blue") +
        geom_line(aes(y = upp975), linetype = "dashed", color = "red", alpha = 0.5) +
        labs(x = "Age", y = "Titre (RLU/mL)") +
        theme_minimal() +
        scale_y_continuous(limits =common_y_limits,
                           #breaks = c(1:6), 
                           labels = log10_to_exp) + 
        theme_universal(base_size = plot_basesize)
    
    # Store the plot in the list with the antigen name as the key
    
    plots2[[a]] <- p2 
    
    
    
}

plot_grid2 <- do.call(grid.arrange, c(plots2, ncol = 4))

joined_plot <- grid.arrange(plot_grid, plot_grid2, ncol =1)

joined_plot_IgA <- joined_plot



