# Title: Event risk modelling using oral fluid IgA levels with mixed-effects logistic regression analysis
# Version: 1.0
# Date: 2025-07-24
# Author: Dr Alexander J Keeley
# Inputs:
#   - R_objects/OF_IgA_titres.RDS
#   - R_objects/SpyCATS_incidence_df.RDS
# Outputs:
#   - Event risk probability plots by titre threshold
#   - Mixed-effects logistic regression tables (crude and adjusted)
#   - Forest plots of odds ratios across antigens and covariates
#   - Combined figures suitable for publication

# Description:
# This script evaluates how oral fluid IgA  titres relate to the 45 day risk of a Strep A event using mixed-effects logistic regression. It:
# 1. Estimates the probability of a future event within a 45-day window for iterative titre thresholds.
# 2. Fits mixed effects logistic to model the probability of an event as a function of IgA levlel, accounting for participant and household clustering.
# 3. Visualizes model-predicted risk and empirical incidence curves across antigen targets.
# 4. Performs adjusted regression analyses including age group, sex, and household size.
# 5. Produces publication-ready figures including combined forest plots and titre-risk prediction panels.



# Setup analysis environment
source("scripts/setup_environment.R")
source("scripts/load_functions.R")


# Set output directory

output_dir <- "R_output/"


# load packages:

shelf(officer,flextable)


##########################################################
########## Event Probability at next visit Estimation ####
##########################################################

# Function: plot_probabilites

# Description:
# This function estimates and visualizes the probability of experiencing a Strep A disease event at the next scheduled visit 
# across a range of IgA titre thresholds. The probability curve is overlaid with the proportion of the population contributing 
# to each threshold bin.

# Inputs:
# - path_to_titre.df: Path to RDS file containing titre data (log10 RLU/mL format expected).
# - sample: Character label for sample type (e.g., "OF").
# - class: Character label for antibody class (e.g., "IgA").
# - next_event_window: Time window (in days) to define "next event" (default = 45 - given sampling framework in study).
# - var_name: Name of titre variable (default = "titre").
# - antigen: Specific antigen to filter for (e.g., "SLO").
# - slice_window: Titre threshold increment (default = 0.1).

# Output:
# - A named list containing:
#   - plot: A ggplot object visualizing:
#       - The probability of an event within the time window vs titre threshold.
#       - Alpha shading scaled to the proportion of individuals in each threshold bin.




plot_probabilites <- function(path_to_titre.df, sample, class, next_event_window = 45, var_name = "titre", antigen, df = 1, slice_window = 0.1, protective_threshold = 0.67) {
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- read_titres(path_to_titre.df, var_name) %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
        group_by(Antigen) %>%
        ungroup()
    
    # Filter for specific antigen
    antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid)) 
    
    pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")
    
    # Prepare event data
    fun_df <- pos_incidence_zero %>%
        select(pid, date, gas_event) %>%
        arrange(pid, date) %>%
        group_by(pid) %>%
        mutate(
            next_date = lead(date),
            check = next_date - date,
            event_next_n = case_when(
                next_date - date <= next_event_window & lead(gas_event) == 1 ~ 1,
                next_date  - date <= next_event_window & lead(gas_event) == 0 ~ 0,
                next_date -date > next_event_window ~ NA
            )) %>%
        select(-next_date) %>%
        ungroup()
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age)) %>%
        group_by(pid) %>%
        fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    # Remove rows with NA in titre or event_next_n
    final_df_4 <- na.omit(final_df_3, cols = c("titre", "event_next_day"))
    
    # Calculate and adjust thresholds for analysis
    min_titre <- min(final_df_4$titre, na.rm = TRUE)
    max_titre <- max(final_df_4$titre, na.rm = TRUE)
    
    # Round down/up to the nearest multiple of slice_window
    rounded_min <- floor(min_titre / slice_window) * slice_window
    rounded_max <- ceiling(max_titre / slice_window) * slice_window
    
    thresholds <- seq(from = rounded_min, to = rounded_max, by = slice_window)
    
    # Initialize a data frame to store results
    probabilities <- data.frame(Threshold = numeric(), Probability = numeric(), percent = numeric())
    
    for (threshold in thresholds) {
        # Subset the data for rows where titre >= threshold
        subset_df <- final_df_4 %>% filter(titre >= threshold)
        
        if (nrow(subset_df) > 0) {
            prob <- sum(subset_df$event_next_n == 1) / nrow(subset_df)
        } else {
            prob <- NA
        }
        
        n_obs <- nrow(subset_df) / nrow(final_df_4) * 100
        probabilities <- rbind(probabilities, data.frame(Threshold = threshold, Probability = prob, percent= n_obs))
    }
    
    max_threshold <- max(probabilities$Threshold[probabilities$percent >= 97.5])
    
    
    ########## Plot it #############
    plot <- ggplot() +
        geom_col(data = probabilities, aes(x = Threshold, y = Probability, alpha = percent), fill = "#138aa8") +
        #    ylim(0, 0.15) +
        xlim(max_threshold, max(thresholds)) +
        labs(title = paste(antigen),
             x = "IgA level (log10 RLU/mL)",
             y = paste0("Proportion experiencing event within ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        theme_universal(base_size = plot_basesize)
    
    results_list <- list("plot" = plot)
    
    return(results_list)
}



# OF

Antigen_list <- c("GAC", "SLO" , "SpyAD",  "SpyCEP")

plot_list <- list()

# Loop through each Anitgen and apply the function: plot_probabilites


for (ag in Antigen_list) {
    
    results <- plot_probabilites(path_to_titre.df = "R_objects/OF_IgA_titres.RDS",
                                 sample = "OF",
                                 class = "IgA",
                                 var_name = "log_RLU_titres",
                                 antigen = ag,
                                 slice_window = 0.1)
    
    
    
    # Extract and store the table and plot
    plot_list[[ag]] <- results$plot
}

# Combine the plots 
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 5) # adjust ncol as needed

# Print the combined plot

print(combined_plot)


################################################
########## mixed effect regression #############
################################################


# This function performs a regression analysis  
#

# Input:
# - Antibody titre data including participant IDs, dates, and titres.
# - Event incidence data ("R_objects/SpyCATS_incidence_df.RDS") with dates of events and participant IDs.
# - Demographic data

# Description:
# - Data cleaning and preparation: Merges antibody titre data with demographic and event data, filters by antigen, and calculates variable (next at next visit within a given time period).
# - Analysis: Calculates titre thresholds, fits a generalized linear mixed-effects model to evaluate the effect of titres on event_next_n. 

# Output:
# - Regression table summarizing the effects of titre levels on event probability.
# - Visualization plot showing the relationship between titre levels and event probability, with mixed effects regression + confidence intervals plotted.

#### how many people: 


mixed_effects_protection_glmer <- function(path_to_titre.df, sample, class, next_event_window = 45, var_name = "titre", antigen, df = 1, slice_window = 0.1, protective_threshold = 0.67) {
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- read_titres(path_to_titre.df, var_name) %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
        group_by(Antigen) %>%
        ungroup()
    
    # Filter for specific antigen
    antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid)) 
    
    pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")
    
    # Prepare event data
    fun_df <- pos_incidence_zero %>%
        select(pid, date, gas_event) %>%
        arrange(pid, date) %>%
        group_by(pid) %>%
        mutate(
            next_date = lead(date),
            check = next_date - date,
            event_next_n = case_when(
                next_date - date <= next_event_window & lead(gas_event) == 1 ~ 1,
                next_date  - date <= next_event_window & lead(gas_event) == 0 ~ 0,
                next_date -date > next_event_window ~ NA
            )) %>%
        select(-next_date) %>%
        ungroup()
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age)) %>%
        group_by(pid) %>%
        fill(titre) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    # Remove rows with NA in titre or event_next_day
    final_df_4 <- na.omit(final_df_3, cols = c("titre", "event_next_day"))
    
    # Calculate and adjust thresholds for analysis
    min_titre <- min(final_df_4$titre, na.rm = TRUE)
    max_titre <- max(final_df_4$titre, na.rm = TRUE)
    
    # Round down/up to the nearest multiple of slice_window
    rounded_min <- floor(min_titre / slice_window) * slice_window
    rounded_max <- ceiling(max_titre / slice_window) * slice_window
    
    thresholds <- seq(from = rounded_min, to = rounded_max, by = slice_window)
    
    # Initialize a data frame to store results
    probabilities <- data.frame(Threshold = numeric(), Probability = numeric(), percent = numeric())
    
    for (threshold in thresholds) {
        # Subset the data for rows where titre >= threshold
        subset_df <- final_df_4 %>% filter(titre >= threshold)
        
        if (nrow(subset_df) > 0) {
            prob <- sum(subset_df$event_next_n == 1) / nrow(subset_df)
        } else {
            prob <- NA
        }
        
        n_obs <- nrow(subset_df) / nrow(final_df_4) * 100
        probabilities <- rbind(probabilities, data.frame(Threshold = threshold, Probability = prob, percent= n_obs))
    }
    
    max_threshold <- max(probabilities$Threshold[probabilities$percent >= 97.5])
    
    # Fit the logistic regression model
    model <- lme4::glmer(event_next_n ~ titre + (1 | pid) + (1 | hid), data = final_df_3, family = binomial)
    
    # Extract AIC
    model_aic <- AIC(model)
    print(paste(antigen,":", AIC(model)))
    
    tb1 <- model %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) 
    
    tb1$table_body$label <- paste0(antigen)
    
    # Extract Odds Ratios, CIs, and p-value
    model_summary <- summary(model)
    ORs <- exp(model_summary$coefficients[, "Estimate"])
    CIs <- exp(confint(model, parm = "titre", method = "Wald"))
    p_value <- coef(summary(model))["titre", "Pr(>|z|)"]  # Extract p-value for titre
    
    # Check if p-value is significant (less than 0.05) and format accordingly
    p_label <- if (p_value < 0.05) {
        paste0("p = ", sprintf("%.3f", p_value), "*")
    } else {
        paste0("p = ", sprintf("%.3f", p_value))
    }
    
    # Create a sequence of titre values for predictions
    new_data <- data.frame(titre = seq(min(final_df_3$titre, na.rm = TRUE), 
                                       max(final_df_3$titre, na.rm = TRUE), length.out = 100))
    new_data$prob <- predict(model, newdata = new_data, type = "response", re.form = NA)
    
    conf_int <- predict(model, newdata = new_data, type = "link", re.form = NA, se.fit = TRUE)
    new_data$lower <- plogis(conf_int$fit - 1.96 * conf_int$se.fit)
    new_data$upper <- plogis(conf_int$fit + 1.96 * conf_int$se.fit)
    
    
    ### 
    x_min <- max_threshold  # You already use this
    x_max <- max(thresholds)
    
    
    ########## Plot it #############
    plot1 <- ggplot(data = probabilities, aes(x = Threshold, y = Probability)) +
        geom_col(fill = "#138aa8") +
        xlim(max_threshold, max(thresholds)) +
        labs(title = paste(antigen),
             x = paste0("IgA threshold (log10 RLU/mL)"),
             y = paste0("Proportion with event within ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        xlim(x_min, x_max) +   # unified x-scale
        theme_universal(base_size = plot_basesize)
    
    
    
    
    plot2 <- ggplot(new_data, aes(x = titre, y = prob)) +
        geom_line(color = "blue") +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
        ylim(0, 0.15) +
        xlim(x_min, x_max) +   # unified x-scale
        labs(title = paste(antigen),
             x = paste0("OF IgA level (log10 RLU/mL)"),
             y = paste0("Probability of event in ", next_event_window, " days")) +
        guides(fill = "none", alpha = "none") +
        theme_minimal() +
        annotate("text", x = max_threshold + 0.5, y = 0.10, 
                 label = paste0("OR: ", round(ORs[2], 2), 
                                "\n95% CI: [", round(CIs[1, 1], 2), ", ", round(CIs[1, 2], 2), "]", 
                                "\n", p_label), 
                 hjust = 0,
                 size = 5) +
        theme_universal(base_size = plot_basesize)
    
    
    plot3 <- final_df_3 %>%
        ggplot(
            aes(
                x = titre
            )
        ) +
        geom_density(fill = StrepA_colscheme[antigen], alpha = 0.8) +
        labs(
            x = "OF IgA level (log10 RLU/mL)") +
        
        theme_minimal() +
        xlim(x_min, x_max) +   # unified x-scale
        theme_universal(base_size = plot_basesize)
    
    plot3
    
    plot <- cowplot::plot_grid(plot1,plot3,plot2,
                               ncol = 1,
                               rel_heights = c(1,0.4,1))
    
    results_list <- list("table" = tb1, "plot" = plot, "AIC_glmer" = model_aic)
    
    return(results_list)
}



# OF IgA

Antigen_list <- c("GAC", "SLO" , "SpyAD",  "SpyCEP")
#Antigen_list <- c("SpyAD",  "SpyCEP")


table_list <- list()
plot_list <- list()


# Initialize a dataframe to store AIC values
AIC_glmer <- data.frame(Antigen = character(), AIC_glmer = numeric(), stringsAsFactors = FALSE)



# Loop through each Anitgen and apply the function: mixed_effects_protection_glmer
for (ag in Antigen_list) {
    
    results <- mixed_effects_protection_glmer(path_to_titre.df = "R_objects/OF_IgA_titres.RDS",
                                              sample = "OF",
                                              class = "IgA",
                                              var_name = "log_RLU_titres",
                                              antigen = ag,
                                              slice_window = 0.1)
    
    
    
    # Extract and store the table and plot
    table_list[[ag]] <- results$table
    plot_list[[ag]] <- results$plot
    # Append the AIC to the dataframe
    AIC_glmer <- rbind(AIC_glmer, data.frame(Antigen = ag, AIC_glmer = results$AIC_glmer))
    
}

# Combine the gtsummary tables
combined_table <- tbl_stack(table_list)

# Combine the plots 
combined_plot <- patchwork::wrap_plots(plot_list, ncol = 4) # adjust ncol as needed



# Print the combined table and plot
print(combined_table)
print(combined_plot)
print(AIC_glmer)





###############################################################
############ Adjusted Mixed-Effects Logistic Model ############
###############################################################

# Function: glmer_adjusted

# Description:
# This function fits an adjusted mixed-effects logistic regression model evaluating the relationship between antibody titres 
# and the short-term probability of a Strep A event, while controlling for key demographic covariates (age group, sex, and household size).

# Inputs:
# - path_to_titre.df: Path to RDS file containing titre data (e.g., OF IgA titres).
# - sample: String indicating the sample type (e.g., "OF").
# - class: String for antibody class (e.g., "IgA").
# - var_name: Column name for titre variable (given both AUC and RLU titres were generated) (default = "titre").
# - antigen: Antigen target to model (e.g., "SpyCEP").
# - next_event_window: Time in days to define a "next event" (default = 45).

# Model details:
# - Random intercepts for participant (`pid`) and household (`hid`) account for repeated measures and shared exposure risk.
# - Age group is relevelled to use "Over 18 years" as the reference category.

# Output:
# - A named list containing:
#   - tbl: A gtsummary regression table (exponentiated coefficients with bolded p-values).
#   - or_df: A tidy dataframe with:
#       - Odds ratios (ORs), 95% confidence intervals, and p-values
#       - Harmonised term names and groupings for use in forest plots



glmer_adjusted <- function(path_to_titre.df, sample, class, var_name = "titre", antigen, next_event_window = 45) {
    
    # Read titre data and merge with start dates and age, adjust visit dates, and categorize titres
    fun_titres <- read.titres(path_to_titre.df, var_name) %>%
        left_join(start_dates) %>%
        left_join(age) %>%
        mutate(visit_date = as.numeric(ymd(as.character(visit_date))) - entry_1) %>%
        group_by(Antigen) %>%
        ungroup()
    
    # Filter for specific antigen
    antigen_df <- fun_titres %>% filter(Antigen == antigen & !is.na(pid)) 
    
    pos_incidence_zero <- readRDS("R_objects/SpyCATS_incidence_df.RDS")
    
    # Prepare event data
    fun_df <- pos_incidence_zero %>%
        select(pid, date, gas_event) %>%
        arrange(pid, date) %>%
        group_by(pid) %>%
        mutate(
            next_date = lead(date),
            check = next_date - date,
            event_next_n = case_when(
                next_date - date <= next_event_window & lead(gas_event) == 1 ~ 1,
                next_date  - date <= next_event_window & lead(gas_event) == 0 ~ 0,
                next_date -date > next_event_window ~ NA
            )) %>%
        select(-next_date) %>%
        ungroup()
    
    
    final_df_3  <-fun_df %>%
        left_join(antigen_df %>% 
                      select(pid, date = visit_date, titre, age_grp, sex, hhsize)) %>%
        group_by(pid) %>%
        fill(titre,age_grp,sex,hhsize) %>%
        ungroup() %>%
        mutate(hid = substring(pid, 1, 3))
    
    final_df_5 <-
        final_df_3
    
    final_df_5$age_grp <- relevel(final_df_5$age_grp, ref = "Over 18 years")
    
    # Fit the logistic regression model
    model_above_only <- lme4::glmer(event_next_n ~ titre
                                    +age_grp 
                                    +sex
                                    +hhsize
                                    + (1 | pid) + (1 | hid), data = final_df_5, family = binomial)
    
    tb1 <- model_above_only %>%
        tbl_regression(exponentiate = TRUE) %>%
        add_nevent() %>%
        bold_p(t = 0.05) %>%
        modify_header(list(label ~ paste(sample, class, "titre"))) 
    
    
    # Extract the summary of the model
    model_summary <- summary(model_above_only)
    
    # Calculate confidence intervals
    conf_intervals <- confint(model_above_only, parm = "beta_", method = "Wald")  # Wald CIs are common for GLMMs
    
    # Create a dataframe with hazard ratios (exponentiated coefficients), confidence intervals, and p-values
    or_df <- data.frame(
        term = rownames(model_summary$coefficients),  # Variable names
        estimate = exp(model_summary$coefficients[, "Estimate"]),  # OR (exp(coef))
        conf.low = exp(conf_intervals[, 1]),  # Lower bound of 95% CI (exponentiated)
        conf.high = exp(conf_intervals[, 2]),  # Upper bound of 95% CI (exponentiated)
        p.value = model_summary$coefficients[, "Pr(>|z|)"]  # P-values
    ) %>%
        mutate(
            # Add descriptive labels for the variables
            variable_label = case_when(
                grepl("titre", term) ~ "Titre",
                grepl("age_grp", term) ~ "Age group",
                grepl("hhsize", term) ~ "Household size",
                grepl("sexFemale", term) ~ "Sex",
                TRUE ~ NA_character_  # Default to NA if no match
            ),
            # Simplify term names for clarity
            term = case_when(
                grepl("titre", term) ~ "Titre",
                grepl("age_grp< 2 years", term) ~ "< 2 years",
                grepl("age_grp12-18 years", term) ~ "12-18 years",
                grepl("age_grp2-4 years", term) ~ "2-4 years",
                grepl("age_grp5-11 years", term) ~ "5-11 years",
                grepl("sexFemale", term) ~ "Female",
                grepl("hhsize", term) ~ "Household size",
                TRUE ~ term  # Retain term if it doesn't match any pattern
            )
        )
    
    # Add an antigen column to the dataframe
    or_df$antigen <- antigen
    
    
    results_list <- list("tbl" = tb1, "or_df" = or_df)
    
    return(results_list)
}


tbl_list <- list()
df_list <- list()
adjusted_model_output <- tibble()

Antigen_list <- c("SpyAD", "SpyCEP")
# Loop through each Anitgen and apply the function: mixed_effects_protection_glmer
for (ag in Antigen_list) {
    
    results <- glmer_adjusted(path_to_titre.df = "R_objects/OF_IgA_titres.RDS",
                                        sample = "OF",
                                        class = "IgA",
                                        var_name = "log_RLU_titres",
                                        antigen = ag)
    
    
    
    # Extract and store the table and plot
    tbl_list[[ag]] <- results$tbl
    df_list[[ag]] <- results$or_df
    
    
    adjusted_model_output <- bind_rows(adjusted_model_output,  results$or_df)
}

gtsummary::tbl_merge(tbl_list)

adjusted_model_output


######################################
### arrange data for a forest plot ###
######################################

new_rows <- unique(adjusted_model_output$antigen) %>%
    expand.grid(antigen = ., term = c("Over 18 years (ref)", "Male (ref)")) %>%
    mutate(
        estimate = NA, 
        conf.low = NA, 
        conf.high = NA,
        p.value = NA,
        variable_label = case_when(
            term == "Below threshold (ref)" ~ "Relation to titre threshold",
            term == "Over 18 years (ref)" ~ "Age group",
            term == "Male (ref)" ~ "Sex"
        )
    )

final_fp_df<- bind_rows(new_rows,adjusted_model_output)

as.vector(levels(as.factor(final_fp_df$term)))

levels(final_fp_df$term) <- gsub("hhsize", "Household size", levels(final_fp_df$term))


final_fp_df <-final_fp_df %>%
    mutate(
        term = factor(term, levels = c(
            "Titre", 
            "Over 18 years (ref)",
            "< 2 years",
            "2-4 years", 
            "5-11 years", 
            "12-18 years", 
            "Male (ref)",
            "Female",
            "Household size" 
        ))
    )

final_fp_df  %>%
    filter(term != "(Intercept)") %>% # Exclude the intercept
    ggplot(aes(x = estimate, y = factor(term, levels = c(
        "Titre", 
        "Over 18 years (ref)",
        "< 2 years",
        "2-4 years", 
        "5-11 years", 
        "12-18 years", 
        "Male (ref)",
        "Female",
        "Household size"
    )))) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    #  facet_wrap(~ antigen, scales = "free") + # Facet by antigen
    labs(
        title = "",
        x = "Odds ratios with 95% CIs",
        y = ""
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    # Add a point for the reference level at x = 1, y = position of reference variable
    
    theme_minimal() +
    geom_point(data = new_rows %>% filter(!antigen %in% c("DNAseB", "GAC")),
               aes(x = 1, y = term), color = "blue", size = 3) +
    # Separate age_grp and event_type under subheadings on y-axis, including reference level
    facet_grid(cols = vars(antigen), rows = vars(variable_label), scales = "free", space = "free_y") +
    scale_x_log10() + 
    theme(
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18, angle = 0), 
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        axis.text.x =element_text(size=12, angle = 90),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))


# Ensure 'term' and 'variable_label' are ordered as intended
final_fp_df <- final_fp_df %>%
    mutate(
        term = factor(term, levels = c(
            "Titre", "Over 18 years (ref)", 
            "12-18 years", "5-11 years", "2-4 years", "< 2 years", 
            "Male (ref)", "Female", "Household size"
        )),
        variable_label = forcats::fct_relevel(variable_label, 
                                              "Titre", "Age group", "Sex", "Household size")
    )

# Plot the forrest plot 

figure <- final_fp_df %>%
    filter(term != "(Intercept)") %>% # Exclude the intercept
    ggplot(aes(x = estimate, y = term)) +  # term already defined with levels
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    labs(
        title = "",
        x = "Odds ratios with 95% CIs",
        y = ""
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    theme_minimal() +
    geom_point(
        data = new_rows %>% filter(!antigen %in% c("DNAseB", "GAC")),
        aes(x = 1, y = term), color = "blue", size = 3
    ) +
    facet_grid(
        cols = vars(antigen),
        rows = vars(forcats::fct_relevel(variable_label, 
                                         "Titre", "Age group", "Sex", "Household size")), 
        scales = "free", 
        space = "free_y"
    ) +
    scale_x_log10() +  # Log scale for odds ratios
    theme(
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18, angle = 0),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)
    ) #+
# theme_universal(base_size = plot_basesize)

compiled_plot <-
    
    plot_grid(
        combined_plot,
        figure, 
        rel_heights = c(0.7,0.2),
        ncol = 1)



compiled_plot

