# Install and load Reticulate + Miniconda  and import Python module 
pipeline <- reticulate::import("pipeline") 


## Control pipeline

### Run the group level pipeline
res_con <- pipeline$group_pipeline(
  
  # input/output paths
  vhdr_files = "/Users/dariagref/Desktop/Data/raw",
  log_files = "/Users/dariagref/Desktop/Data/log",
  output_dir = "/Users/dariagref/Desktop/Data/output_con",
  
  # downsampling
  downsample_sfreq = 250.0,
  
  # interpolate bad channels
  bad_channels = "auto",
  
  # re-referencing default
  
  #filtering
  highpass_freq = 0.1,
  lowpass_freq = 30.0,
  
  # epoching
  triggers = c(201:208, 211:218),
  components = list(
    "name" = list("N2", "P3b"),
    "tmin" = list(0.25, 0.4),
    "tmax" = list(0.35, 0.55),
    "roi" = list(
      c("FC1", "FC2", "C1", "C2", "Cz"),
      c("CP3", "CP1", "CPz", "CP2", "CP4", "P3", "Pz", "P4", "PO3", "POz", "PO4")
    )
  ),
  
  # baseline correction default
  
  #deactivate peak-to-peak rejection
  reject_peak_to_peak = NULL,
  
  # combine trials/evokeds
  average_by = c("n_b")
)

### Extract results
trials_con <- res_con[[1]]
evokeds_con <- res_con[[2]]
config_con <- res_con[[3]]

# Single trial data frame
print(trials_con)

# Single trial N2 mean amplitudes
N2_con <- ggplot(trials_con, aes(x = N2, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_con$N2), max(trials_con$N2)))

# Single trial P3b mean amplitudes
P3b_con <- ggplot(trials_con, aes(x = P3b, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_con$P3b), max(trials_con$P3b)))


## Rejection Pipeline

### Run the group level pipeline
res_rej <- pipeline$group_pipeline(
  
  # input/output paths
  vhdr_files = "/Users/dariagref/Desktop/Data/raw",
  log_files = "/Users/dariagref/Desktop/Data/log",
  output_dir = "/Users/dariagref/Desktop/Data/output_rej",
  
  # downsampling
  downsample_sfreq = 250.0,
  
  # interpolate bad channels
  bad_channels = "auto",
  
  # re-referencing default
  
  #filtering
  highpass_freq = 0.1,
  lowpass_freq = 30.0,
  
  # epoching
  triggers = c(201:208, 211:218),
  components = list(
    "name" = list("N2", "P3b"),
    "tmin" = list(0.25, 0.4),
    "tmax" = list(0.35, 0.55),
    "roi" = list(
      c("FC1", "FC2", "C1", "C2", "Cz"),
      c("CP3", "CP1", "CPz", "CP2", "CP4", "P3", "Pz", "P4", "PO3", "POz", "PO4")
    )
  ),
  
  # baseline correction default
  
  # peak-to-peak rejection
  reject_peak_to_peak = 100.0,
  
  # combine trials/evokeds
  average_by = c("n_b")
)

### Extract results
trials_rej <- res_rej[[1]]
evokeds_rej <- res_rej[[2]]
config_rej <- res_rej[[3]]

# Single trial data frame
print(trials_rej)

# Single trial N2 mean amplitudes
N2_rej <- ggplot(trials_rej, aes(x = N2, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_rej$N2), max(trials_rej$N2)))

# Single trial P3b mean amplitudes
P3b_rej <- ggplot(trials_rej, aes(x = P3b, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_rej$P3b), max(trials_rej$P3b)))


## ICA Pipeline

### Run the group level pipeline
res_ica <- pipeline$group_pipeline(
  
  # input/output paths
  vhdr_files = "/Users/dariagref/Desktop/Data/raw",
  log_files = "/Users/dariagref/Desktop/Data/log",
  output_dir = "/Users/dariagref/Desktop/Data/output_ica",
  
  # downsampling
  downsample_sfreq = 250.0,
  
  # VEOG and HEOG channels
  veog_channels = "auto",
  heog_channels = "auto",
  
  # interpolate bad channels
  bad_channels = "auto",
  
  # re-referencing default
  
  # preprocessing option
  ica_method = "fastica",
  ica_n_components = 0.99,
  
  #filtering
  highpass_freq = 0.1,
  lowpass_freq = 30.0,
  
  # epoching
  triggers = c(201:208, 211:218),
  components = list(
    "name" = list("N2", "P3b"),
    "tmin" = list(0.25, 0.4),
    "tmax" = list(0.35, 0.55),
    "roi" = list(
      c("FC1", "FC2", "C1", "C2", "Cz"),
      c("CP3", "CP1", "CPz", "CP2", "CP4", "P3", "Pz", "P4", "PO3", "POz", "PO4")
    )
  ),
  
  # baseline correction default
  
  # peak-to-peak rejection
  reject_peak_to_peak = 200.0,
  
  # combine trials/evokeds
  average_by = c("n_b")
)

### Extract results
trials_ica <- res_ica[[1]]
evokeds_ica <- res_ica[[2]]
config_ica <- res_ica[[3]]

# Single trial data frame
print(trials_ica)

# Single trial N2 mean amplitudes
N2_ica <- ggplot(trials_ica, aes(x = N2, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_ica$N2), max(trials_ica$N2)))

# Single trial P3b mean amplitudes
P3b_ica <- ggplot(trials_ica, aes(x = P3b, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_ica$P3b), max(trials_ica$P3b)))


## BESA Pipeline

### Run the group level pipeline
res_besa <- pipeline$group_pipeline(
  
  # input/output paths
  vhdr_files = "/Users/dariagref/Desktop/Data/raw",
  log_files = "/Users/dariagref/Desktop/Data/log",
  output_dir = "/Users/dariagref/Desktop/Data/output_besa",
  
  # downsampling
  downsample_sfreq = 250.0,
  
  # VEOG and HEOG channels
  veog_channels = "auto",
  heog_channels = "auto",
  
  # interpolate bad channels
  bad_channels = "auto",
  
  #re-referencing default
  
  # preprocessing option
  besa_files = "/Users/dariagref/Desktop/Data/cali",
  
  #filtering
  highpass_freq = 0.1,
  lowpass_freq = 30.0,
  
  # epoching
  triggers = c(201:208, 211:218),
  components = list(
    "name" = list("N2", "P3b"),
    "tmin" = list(0.25, 0.4),
    "tmax" = list(0.35, 0.55),
    "roi" = list(
      c("FC1", "FC2", "C1", "C2", "Cz"),
      c("CP3", "CP1", "CPz", "CP2", "CP4", "P3", "Pz", "P4", "PO3", "POz", "PO4")
    )
  ),
  
  # baseline correction default
  
  # peak-to-peak rejection
  reject_peak_to_peak = 200.0,
  
  # combine trials/evokeds
  average_by = c("n_b")
)

### Extract results
trials_besa <- res_besa[[1]]
evokeds_besa <- res_besa[[2]]
config_besa <- res_besa[[3]]

# Single trial data frame
print(trials_besa)

# Single trial N2 mean amplitudes
N2_besa <-ggplot(trials_besa, aes(x = N2, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_besa$N2), max(trials_besa$N2)))

# Single trial P3b mean amplitudes
P3b_besa <- ggplot(trials_besa, aes(x = P3b, fill = n_b)) +
  geom_density(alpha = 0.5) +
  theme_classic(base_size = 30) +
  coord_cartesian(xlim = c(min(trials_besa$P3b), max(trials_besa$P3b)))

# Combined single trial mean amplitudes plot

library(cowplot)

n2_plot <- plot_grid(N2_con, N2_rej, N2_ica, N2_besa, labels = c('Control', 'Rejection', 'ICA', 'MSEC'), 
                     label_size = 20, label_colour = "#95979c", scale = 0.9)
p3b_plot <- plot_grid(P3b_con, P3b_rej, P3b_ica, P3b_besa, labels = c('Control', 'Rejection', 'ICA', 'MSEC'),
                      label_size = 20, label_colour = "#95979c", scale = 0.9)
n2_plot
p3b_plot




# Calculate the data quality metrics

## Combine dataframes

library(dplyr)

combined_trials <- bind_rows(
  trials_con %>%
    mutate(Pipeline = "Control"),
  trials_rej %>%
    mutate(Pipeline = "Rejection"),
  trials_ica %>%
    mutate(Pipeline = "ICA"),
  trials_besa %>%
    mutate(Pipeline = "MSEC")
)

combined_trials <- bind_rows(
  trials_con %>%
    mutate(Pipeline = "Control", participant_id = as.integer(participant_id)),
  trials_rej %>%
    mutate(Pipeline = "Rejection", participant_id = as.integer(participant_id)),
  trials_ica %>%
    mutate(Pipeline = "ICA", participant_id = as.integer(participant_id)),
  trials_besa %>%
    mutate(Pipeline = "MSEC", participant_id = as.integer(participant_id))
)

## Group by participant, condition and pipeline and calculate mean, sd, and SME

grouped_data <- combined_trials %>%
  group_by(participant_id, n_b, Pipeline) 

summary_stats <- grouped_data %>%
  na.omit() %>%
  summarize(
    mean_N2 = mean(N2, na.rm = TRUE),
    sd_N2 = sd(N2, na.rm = TRUE),
    mean_P3b = mean(P3b, na.rm = TRUE),
    sd_P3b = sd(P3b, na.rm = TRUE),
    n = n(),
    SME_N2 = sd_N2/sqrt(n),
    SME_P3b = sd_P3b/sqrt(n),
    .groups = 'keep'
  )

View(summary_stats)


## SME Plot N2

library(dplyr)

# Add the pipeline_n_b levels to summary_stats
summary_stats <- summary_stats %>%
  mutate(pipeline_n_b = paste(Pipeline, n_b, sep = "\n"),
         pipeline_n_b = factor(pipeline_n_b, levels = c("Control\nblurr", "Control\nnormal", "Rejection\nblurr", "Rejection\nnormal", "ICA\nblurr", "ICA\nnormal", "MSEC\nblurr", "MSEC\nnormal")))

# Define the desired colors for each level
color_mapping <- c(
  "Control\nblurr" = "#E69F00",
  "Control\nnormal" = "#FFCF00",
  "Rejection\nblurr" = "#009E73",
  "Rejection\nnormal" = "#81CFA2",
  "ICA\nblurr" = "#0072B2",
  "ICA\nnormal" = "#4DA1E5",
  "MSEC\nblurr" = "#CC79A7",
  "MSEC\nnormal" = "#FCBBD8"
)

# Plot
summary_stats %>%
  mutate(pipeline_n_b = paste(Pipeline, n_b, sep = "\n"),
         pipeline_n_b = factor(pipeline_n_b, levels = names(color_mapping))) %>%
  arrange(participant_id, pipeline_n_b) %>%
  ggplot(aes(x = pipeline_n_b, y = SME_N2, color = pipeline_n_b, group = participant_id)) +
  geom_line(position = position_jitter(width = 0.4, seed = 42), alpha = 0.1) +
  geom_point(size = 1.8, shape = 20, position = position_jitter(width = 0.4, seed = 42)) +
  labs(title = "N2", x = "Pipeline", y = "SME") +
  theme_classic(base_size = 14) +
  scale_color_manual(values = color_mapping) +  # Use the color_mapping here
  guides(fill = "none") +
  theme(legend.position = "none")

## SME Plot P3b

# Define the desired colors for each level
color_mapping <- c(
  "Control\nblurr" = "#E69F00",
  "Control\nnormal" = "#FFCF00",
  "Rejection\nblurr" = "#009E73",
  "Rejection\nnormal" = "#81CFA2",
  "ICA\nblurr" = "#0072B2",
  "ICA\nnormal" = "#4DA1E5",
  "MSEC\nblurr" = "#CC79A7",
  "MSEC\nnormal" = "#FCBBD8"
)

# Plot
summary_stats %>%
  mutate(pipeline_n_b = paste(Pipeline, n_b, sep = "\n"),
         pipeline_n_b = factor(pipeline_n_b, levels = names(color_mapping))) %>%
  arrange(participant_id, pipeline_n_b) %>%
  ggplot(aes(x = pipeline_n_b, y = SME_P3b, color = pipeline_n_b, group = participant_id)) +
  geom_line(position = position_jitter(width = 0.4, seed = 42), alpha = 0.1) +
  geom_point(size = 1.8, shape = 20, position = position_jitter(width = 0.4, seed = 42)) +
  labs(title = "P3b", x = "Pipeline", y = "SME") +
  theme_classic(base_size = 14) +
  scale_color_manual(values = color_mapping) +  # Use the color_mapping here
  guides(fill = "none") +
  theme(legend.position = "none")


## SME Linear Mixed Effects Models

library(lmerTest)

# Convert "Pipeline" to a factor variable with the desired levels
summary_stats$Pipeline <- factor(summary_stats$Pipeline, levels = c("Control", "Rejection", "ICA", "MSEC"))

# Convert "n_b" to a factor variable
summary_stats$n_b <- factor(summary_stats$n_b, levels = c("normal", "blurr"))

# Assign the forward difference coding to the "Pipeline" factor
contrasts(summary_stats$Pipeline) <- MASS::contr.sdif(4)

# Contrasts n_b factor
contrasts(summary_stats$n_b) <- contr.sum(2)/2

# Fit the linear mixed-effects models with forward difference coding

### model_N2 <- lmer(SME_N2 ~ n_b + Pipeline + (1 + n_b + Pipeline | participant_id), data = summary_stats)
### Model failed to converge, exclude random slope for visibility
### model_N2 <- lmer(SME_N2 ~ n_b + Pipeline + (1 + Pipeline | participant_id), data = summary_stats)
### Model failed to converge again, exclude random slope for pipeline
model_N2 <- lmer(SME_N2 ~ n_b + Pipeline + (1 | participant_id), data = summary_stats)

### model_P3b <- lmer(SME_P3b ~ n_b + Pipeline + (1 + n_b + Pipeline | participant_id), data = summary_stats)
### Model failed to converge, exclude random slope for visibility
### model_P3b <- lmer(SME_P3b ~ n_b + Pipeline + (1 + Pipeline | participant_id), data = summary_stats)
### Model failed to converge again, exclude random slope for pipeline
model_P3b <- lmer(SME_P3b ~ n_b + Pipeline + (1 | participant_id), data = summary_stats)

# Print the summary of the models
summary(model_N2)
summary(model_P3b)

# Calculate mean values separately by conditions
library(emmeans)

emmeans(model_N2, specs = ~Pipeline + n_b) 
emmeans(model_P3b, specs = ~Pipeline + n_b)




## Number of trials needed to reach significance

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

# Define function
calculate_percentage_significant <- function(data, dependent_variable, num_trials, num_repetitions, significance_threshold) {
  percentage_significant_list <- vector("numeric", length(num_trials))
  
  for (num_trial in num_trials) {
    model_significance <- vector("logical", num_repetitions)
    
    for (i in 1:num_repetitions) {
      grouped_trials <- data %>%
        group_by(participant_id, n_b)
      
      sampled_trials <- grouped_trials %>%
        sample_n(size = num_trial, replace = TRUE, seed = i)
      
      tryCatch({
        model_sampled <- lmer(paste0(dependent_variable, " ~ n_b + (1 + n_b | participant_id)"), data = sampled_trials)
        model_summary <- summary(model_sampled)
        model_p_values <- model_summary$coefficients[, "Pr(>|t|)"]
        
        if (isSingular(model_sampled)) {
          model_significance[i] <- FALSE
        } else {
          if (any(grepl("Model failed to converge", warnings()))) {
            model_significance[i] <- FALSE
          } else {
            model_significance[i] <- any(model_p_values < significance_threshold)
          }
        }
      }, error = function(e) {
        model_significance[i] <- FALSE
      })
    }
    
    percentage_significant <- mean(model_significance) * 100
    percentage_significant_list[num_trial/100] <- percentage_significant
  }
  
  return(percentage_significant_list)
}

# Define the parameters
num_trials <- seq(from = 0, to = 100, by = 10)
num_repetitions <- 1000
significance_threshold <- 0.05

# Call the function for each dataset and dependent variable
percentage_significant_list_con_n2 <- calculate_percentage_significant(trials_con, "N2", num_trials, num_repetitions, significance_threshold)
percentage_significant_list_con_p3b <- calculate_percentage_significant(trials_con, "P3b", num_trials, num_repetitions, significance_threshold)

percentage_significant_list_rej_n2 <- calculate_percentage_significant(trials_rej, "N2", num_trials, num_repetitions, significance_threshold)
percentage_significant_list_rej_p3b <- calculate_percentage_significant(trials_rej, "P3b", num_trials, num_repetitions, significance_threshold)

percentage_significant_list_ica_n2 <- calculate_percentage_significant(trials_ica, "N2", num_trials, num_repetitions, significance_threshold)
percentage_significant_list_ica_p3b <- calculate_percentage_significant(trials_ica, "P3b", num_trials, num_repetitions, significance_threshold)

percentage_significant_list_besa_n2 <- calculate_percentage_significant(trials_besa, "N2", num_trials, num_repetitions, significance_threshold)
percentage_significant_list_besa_p3b <- calculate_percentage_significant(trials_besa, "P3b", num_trials, num_repetitions, significance_threshold)

## N2
# Combine dataframes

library(dplyr)

# Create data frames for each pipeline's results 

df_con_n2 <- data.frame(pipeline = "con", 
                        trials = num_trials, 
                        power = percentage_significant_list_con_n2)

df_rej_n2 <- data.frame(pipeline = "rej", 
                        trials = num_trials, 
                        power = percentage_significant_list_rej_n2)

df_ica_n2 <- data.frame(pipeline = "ica", 
                        trials = num_trials, 
                        power = percentage_significant_list_ica_n2)

df_besa_n2 <- data.frame(pipeline = "besa", 
                         trials = num_trials, 
                         power = percentage_significant_list_besa_n2)

# Combine the individual data frames into one big data frame
combined_df_n2 <- bind_rows(df_con_n2, df_rej_n2, df_ica_n2, df_besa_n2)

# Print the combined data frame
print(combined_df_n2)

## P3b
# Create data frames for each pipeline's results

df_con_p3b <- data.frame(pipeline = "con", 
                         trials = num_trials, 
                         power = percentage_significant_list_con_p3b)

df_rej_p3b <- data.frame(pipeline = "rej", 
                         trials = num_trials, 
                         power = percentage_significant_list_rej_p3b)

df_ica_p3b <- data.frame(pipeline = "ica", 
                         trials = num_trials, 
                         power = percentage_significant_list_ica_p3b)

df_besa_p3b <- data.frame(pipeline = "besa", 
                          trials = num_trials, 
                          power = percentage_significant_list_besa_p3b)

# Combine the individual data frames into one big data frame
combined_df_p3b <- bind_rows(df_con_p3b, df_rej_p3b, df_ica_p3b, df_besa_p3b)

# Print the combined data frame
print(combined_df_p3b)

# Plot N2

ggplot(combined_df_n2, aes(x = trials, y = power, color = pipeline)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Trials",
       y = "Power (%)",
       title = "Power of Different Pipelines - N2",
       color = "Pipeline") +
  scale_color_manual(values = c("#E69F00", "#009E73", "#0072B2", "#CC79A7")) +  # Custom colors
  theme_minimal()

# Plot P3b

combined_df_p3b <- read.csv("/Users/dariagref/Desktop/Data/combined_df_p3b.csv", header = TRUE, sep = "," )

# Create a line plot using ggplot2
ggplot(combined_df_p3b, aes(x = trials, y = power, color = pipeline)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Trials",
       y = "Power (%)",
       title = "Power of Different Pipelines - P3b",
       color = "Pipeline") +
  scale_color_manual(values = c("#E69F00", "#009E73", "#0072B2", "#CC79A7")) +  # Custom colors
  theme_minimal()