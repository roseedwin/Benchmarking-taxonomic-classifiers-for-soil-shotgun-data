library(yardstick)
library(tidyverse)
library(reshape2)
library(dplyr)
library(kableExtra)

cat(paste(rep("_", nchar(" STARTING YOUR SCRIPT FOR species LEVEL analysis ")), collapse = ""), "\n")

df_long <-  read.table("Output/Species_results.tsv", header =T,  sep="\t")
df_long <-  df_long %>% separate(sample_name, into = c("sample_name", "run", "tool"), sep = "\\.")
df_long2 <-   df_long %>%  group_by(species, sample_name, run, tool) %>%
  summarise(rel.abundance = sum(abundance))
df_long2$tool <- gsub("Kraken_gtdb", "Kraken.gtdb", df_long2$tool)

cat("\n Phew, that took a while to get your stupid relative abundance sorted :P \n")

expected <-  df_long2 %>% filter(tool == "insilico") %>%  #filter everything insilico
  rename(expected = rel.abundance)  #rename rel.abundance as expected

observed <-  df_long2 %>% filter(tool != "insilico")%>%  #filter everything not insilico
  rename(observed = rel.abundance) #rename rel.abundance as observed

#pull into different data frames based on cut-off and add column threshold with cut-off value
data_0.0001 <- observed %>% filter_at(vars(5), any_vars(. > 0.0001)) %>% mutate(threshold = "0.0001")
data_0.0005 <- observed %>% filter_at(vars(5), any_vars(. > 0.0005)) %>% mutate(threshold = "0.0005")
data_0.001 <- observed %>% filter_at(vars(5), any_vars(.> 0.001)) %>% mutate(threshold = "0.001")
data_0.005 <- observed %>% filter_at(vars(5), any_vars(.> 0.005)) %>% mutate(threshold = "0.005")
data_0.01 <- observed %>% filter_at(vars(5), any_vars(. > 0.01)) %>% mutate(threshold = "0.01")
data_0.05 <- observed %>% filter_at(vars(5), any_vars(.> 0.05)) %>% mutate(threshold = "0.05")
data_0.1 <- observed %>% filter_at(vars(5), any_vars(. > 0.1)) %>% mutate(threshold = "0.1")
data_0.5 <- observed %>% filter_at(vars(5), any_vars(. > 0.5)) %>% mutate(threshold = "0.5")
data_0 <- observed %>% filter_at(vars(5), any_vars(. > 0)) %>% mutate(threshold = "0")

#Checking the lowest and highest value in each data frame

cat("\n Okay, lets just check if your data frames with different thresholds are secure. \n")
# Get a list of objects in the environment
objects <- ls()
# Filter the objects that start with "data_"
data_frames <- objects[grep("^data_", objects)]
# Iterate over each data frame and check the range of rel.abundance
for (df_name in data_frames) {
  df <- get(df_name)
  min_value <- min(df$observed)
  max_value <- max(df$observed)

  cat("Data frame:", df_name, "\n")
  cat("Minimum value:", min_value, "\n")
  cat("Maximum value:", max_value, "\n\n")
}

#combining all data.frame back
observed <- rbind(data_0,data_0.0001, data_0.0005, data_0.001, data_0.005, data_0.01, data_0.05, data_0.1, data_0.5)

cat("\n Thanks for your patience, We will now process each tool \n")
cat("\n Starting with Kaiju \n")

#Processing each tool differently
KJ_observed <-  subset(observed, tool == "Kaiju")
KJ_observed <-  KJ_observed %>%
  select(species, tool,threshold, observed, sample_name, run) %>%
  pivot_wider(names_from = threshold, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

result <- full_join(KJ_observed, expected, by = c("run", "sample_name", "species"))

# Fill NA values in tool.x with "Kaiju"
result$tool.x <- ifelse(is.na(result$tool.x), "Kaiju", result$tool.x)
# Fill NA values in tool.y with "insilico"
result$tool.y <- ifelse(is.na(result$tool.y), "insilico", result$tool.y)
# Replace values in the specified columns
result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] <-
  ifelse(result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] > 0, "present", "absent")

df_longer.KJ <- result %>%
  pivot_longer(cols = 5:13, names_to = "threshold", values_to = "observed") %>%
  select(species, tool.x, sample_name, run, expected, threshold, observed)  %>%
  mutate(expected = ifelse(is.na(expected) == T , "absent", expected)) %>%
  mutate(observed = ifelse(is.na(observed) == T , "absent", observed))

cat("..Moving to Metaphlan3 \n")
MP3_observed <- subset(observed, tool == "Metaphlan3")
MP3_observed <-  MP3_observed %>%
  select(species, tool,threshold, observed, sample_name, run) %>%
  pivot_wider(names_from = threshold, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

result <- full_join(MP3_observed, expected, by = c("run", "sample_name", "species"))
# Fill NA values in tool.x with "Metaphlan3"
result$tool.x <- ifelse(is.na(result$tool.x), "Metaphlan3", result$tool.x)
# Fill NA values in tool.y with "insilico"
result$tool.y <- ifelse(is.na(result$tool.y), "insilico", result$tool.y)
# Replace values in the specified columns
result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] <-
  ifelse(result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] > 0, "present", "absent")
df_longer.MP3 <- result %>%
  pivot_longer(cols = 5:13, names_to = "threshold", values_to = "observed") %>%
  select(species, tool.x, sample_name, run, expected, threshold, observed)  %>%
  mutate(expected = ifelse(is.na(expected) == T , "absent", expected)) %>%
  mutate(observed = ifelse(is.na(observed) == T , "absent", observed))

cat("...Moving to Metaphlan4 \n")
MP4_observed <- subset(observed, tool == "Metaphlan4")
MP4_observed <-  MP4_observed %>%
  select(species, tool,threshold, observed, sample_name, run) %>%
  pivot_wider(names_from = threshold, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

result <- full_join(MP4_observed, expected, by = c("run", "sample_name", "species"))
# Fill NA values in tool.x with "Metaphlan4"
result$tool.x <- ifelse(is.na(result$tool.x), "Metaphlan4", result$tool.x)
# Fill NA values in tool.y with "insilico"
result$tool.y <- ifelse(is.na(result$tool.y), "insilico", result$tool.y)
# Replace values in the specified columns
result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] <-
  ifelse(result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] > 0, "present", "absent")
df_longer.MP4 <- result %>%
  pivot_longer(cols = 5:13, names_to = "threshold", values_to = "observed") %>%
  select(species, tool.x, sample_name, run, expected, threshold, observed)  %>%
  mutate(expected = ifelse(is.na(expected) == T , "absent", expected)) %>%
  mutate(observed = ifelse(is.na(observed) == T , "absent", observed))

cat("....Moving to Kraken \n")
KR_observed <-  subset(observed, tool == "Kraken")
KR_observed <-  KR_observed %>%
  select(species, tool,threshold, observed, sample_name, run) %>%
  pivot_wider(names_from = threshold, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

result <- full_join(KR_observed, expected, by = c("run", "sample_name", "species"))
# Fill NA values in tool.x with "Kraken"
result$tool.x <- ifelse(is.na(result$tool.x), "Kraken", result$tool.x)
# Fill NA values in tool.y with "insilico"
result$tool.y <- ifelse(is.na(result$tool.y), "insilico", result$tool.y)
# Replace values in the specified columns
result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] <-
  ifelse(result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] > 0, "present", "absent")
df_longer.KR <- result %>%
  pivot_longer(cols = 5:13, names_to = "threshold", values_to = "observed") %>%
  select(species, tool.x, sample_name, run, expected, threshold, observed)  %>%
  mutate(expected = ifelse(is.na(expected) == T , "absent", expected)) %>%
  mutate(observed = ifelse(is.na(observed) == T , "absent", observed))

cat(".....And finally Kraken.gtdb \n")
KG_observed <-  subset(observed, tool == "Kraken.gtdb")
KG_observed <-  KG_observed %>%
  select(species, tool,threshold, observed, sample_name, run) %>%
  pivot_wider(names_from = threshold, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

result <- full_join(KG_observed, expected, by = c("run", "sample_name", "species"))
# Fill NA values in tool.x with "Kraken.gtdb"
result$tool.x <- ifelse(is.na(result$tool.x), "Kraken.gtdb", result$tool.x)
# Fill NA values in tool.y with "insilico"
result$tool.y <- ifelse(is.na(result$tool.y), "insilico", result$tool.y)
# Replace values in the specified columns
result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] <-
  ifelse(result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] > 0, "present", "absent")
df_longer.KG <- result %>%
  pivot_longer(cols = 5:13, names_to = "threshold", values_to = "observed") %>%
  select(species, tool.x, sample_name, run, expected, threshold, observed)  %>%
  mutate(expected = ifelse(is.na(expected) == T , "absent", expected)) %>%
  mutate(observed = ifelse(is.na(observed) == T , "absent", observed))

cat("....Moving to Kraken22 \n")
KR22_observed <-  subset(observed, tool == "Kraken22")
KR22_observed <-  KR22_observed %>%
  select(species, tool,threshold, observed, sample_name, run) %>%
  pivot_wider(names_from = threshold, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

result <- full_join(KR22_observed, expected, by = c("run", "sample_name", "species"))
# Fill NA values in tool.x with "Kraken22"
result$tool.x <- ifelse(is.na(result$tool.x), "Kraken22", result$tool.x)
# Fill NA values in tool.y with "insilico"
result$tool.y <- ifelse(is.na(result$tool.y), "insilico", result$tool.y)
# Replace values in the specified columns
result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] <-
  ifelse(result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] > 0, "present", "absent")
df_longer.KR22 <- result %>%
  pivot_longer(cols = 5:13, names_to = "threshold", values_to = "observed") %>%
  select(species, tool.x, sample_name, run, expected, threshold, observed)  %>%
  mutate(expected = ifelse(is.na(expected) == T , "absent", expected)) %>%
  mutate(observed = ifelse(is.na(observed) == T , "absent", observed))

cat("....Moving to Kraken22 \n")
KR22_observed <-  subset(observed, tool == "Kraken22")
KR22_observed <-  KR22_observed %>%
  select(species, tool,threshold, observed, sample_name, run) %>%
  pivot_wider(names_from = threshold, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

result <- full_join(KR22_observed, expected, by = c("run", "sample_name", "species"))
# Fill NA values in tool.x with "Kraken22"
result$tool.x <- ifelse(is.na(result$tool.x), "Kraken22", result$tool.x)
# Fill NA values in tool.y with "insilico"
result$tool.y <- ifelse(is.na(result$tool.y), "insilico", result$tool.y)
# Replace values in the specified columns
result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] <-
  ifelse(result[, c("0", "0.0001", "0.001", "0.0005", "0.005", "0.01", "0.05", "0.1", "0.5", "expected")] > 0, "present", "absent")
df_longer.KR22 <- result %>%
  pivot_longer(cols = 5:13, names_to = "threshold", values_to = "observed") %>%
  select(species, tool.x, sample_name, run, expected, threshold, observed)  %>%
  mutate(expected = ifelse(is.na(expected) == T , "absent", expected)) %>%
  mutate(observed = ifelse(is.na(observed) == T , "absent", observed))

#Pool all the tools together
observed_expected <-  rbind(df_longer.KG,df_longer.KJ, df_longer.MP3, df_longer.MP4, df_longer.KR, df_longer.KR22) %>%
  rename(tool = tool.x)

#observed_expected <-  read.table("OUTPUT_FILES/observed_expected.tsv", header =T,  sep="\t")

#Make threshold a character for easy downstream analysis
observed_expected$threshold <- as.character(observed_expected$threshold)
observed_expected$threshold <- gsub("1e-04", "0.0001", observed_expected$threshold)
observed_expected$threshold <- gsub("5e-04", "0.0005", observed_expected$threshold)

#Using Amy's script from here on
cat("\n Almost done: Creating the observed_expected2 \n")
observed_expected2 <-  observed_expected %>%
  mutate(tool = factor(tool, levels =c("Kraken", "Kraken.gtdb", "Kaiju", "Metaphlan3", "Metaphlan4", "Kraken22"))) %>%
  mutate(threshold =factor(threshold, levels = c("0", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5"))) %>%
  mutate(run =factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) %>%
  mutate(expected=factor(expected, levels=c("present", "absent"))) %>%
  mutate(observed=factor(observed, levels=c("present", "absent"))) %>%
  dplyr::rename(SampleID = sample_name) %>%
  unite(lib_id, c("run", "tool", "threshold"), sep ="_") %>%
  mutate(lib_id=as.factor(lib_id)) %>%
  group_by(lib_id)

observed_expected2 <-  observed_expected %>%
  mutate(tool = factor(tool, levels =c("Kraken", "Kraken.gtdb", "Kaiju", "Metaphlan3", "Metaphlan4", "Kraken22"))) %>%
  mutate(threshold =factor(threshold, levels = c("0","0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5"))) %>%
  mutate(run =factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) %>%
  mutate(expected=factor(expected, levels=c("present", "absent"))) %>%
  mutate(observed=factor(observed, levels=c("present", "absent"))) %>%
  dplyr::rename(SampleID = sample_name) %>%
  unite(lib_id, c("run", "tool", "threshold"), sep ="_") %>%
  mutate(lib_id=as.factor(lib_id)) %>%
  group_by(lib_id) %>%
  group_split() %>%
  setNames(sort(unique(observed_expected2$lib_id)))

cat("..Moving to confusion table \n")
`%nin%` = Negate(`%in%`)
confusion_table2 <- lapply(observed_expected2 , function(df) {
  df  %>%
    select(-SampleID) %>%
    mutate(expected=factor(expected,levels=c("present", "absent"))) %>%
    mutate(observed=factor(observed,levels=c("present", "absent")))
}
)

cm <-lapply(confusion_table2, yardstick::conf_mat, truth = expected,estimate = observed)

metrics_cm <- lapply(cm, function(df) {df %>%summary() %>%
    select(-.estimator) %>%
    filter(.metric %in%
             c("spec", "sens", "accuracy", "precision", "recall", "f_meas"))})

cat("...The output cm. \n")

output_cm <- metrics_cm %>%
  dplyr::bind_rows(., .id = "lib_id")%>%
  mutate_if(is.numeric, round, digits = 2) %>%
  tidyr::pivot_wider(names_from= .metric, values_from= .estimate)


write.table(output_cm, "Output/output_Species.tsv", sep="\t", quote = FALSE, row.names = F)

cat("\nAnd DONE, I have written the output file for you. \n")
cat("You are WELCOME idiot!! \n")
cat(paste(rep("=", nchar(" And DONE, I have written the output file for you. ")), collapse = ""), "\n")

#Repeat the same for genus level data

#Plots generated 
library(tidyverse)
library(reshape2)
library(dplyr)

categories <- c("Genus", "Species")

# Create a list to store the data frames
df_list <- list()

# Read and process each TSV file
for (category in categories) {
  # Generate the file name
  filename <- paste0("Output/output_", category, ".tsv")

  # Read the TSV file
  df <- read.delim(filename, header = TRUE, sep = "\t")

  # Add the category name as a new column
  df$category <- category

  # Add the processed data frame to the list
  df_list[[category]] <- df
}

# Combine the data frames into a single data frame
combined_df <- do.call(rbind, df_list)
combined_df <-  separate(combined_df, lib_id, into = c("run", "tool", "threshold"), sep = "_")
combined_df$tool <- gsub( "Kraken.gtdb", "custom.Kraken",combined_df$tool)
combined_df$tool <- gsub( "Metaphlan", "MetaPhlAn",combined_df$tool)

df1 <- combined_df  %>%
  mutate(FPR = (1 - spec))

df1$category <- factor(df1$category, levels = c("Species","Genus"))
df1$threshold <- factor(df1$threshold, levels = c("0", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1"))
levels(df1$threshold)[levels(df1$threshold) == "0"] <- "No cut-off"

classifier_colors <- c(
  "custom.Kraken" = "#70A494",
  "Kaiju" = "#BD925A",
  "MetaPhlAn3" =  "#DE8A5A",
  "MetaPhlAn4" =  "#CA562C",
  "Kraken" = "#798234"
)


# Calculate F1 average and standard deviation
df_summary <- df1 %>%
  filter(tool != "Kraken") %>%
  mutate(tool = ifelse(tool == "Kraken22", "Kraken", tool)) %>%
  group_by(threshold, tool, category) %>%
  summarize(
    avg_f_meas = mean(f_meas),
    sd_f_meas = sd(f_meas)
  )

plot <- ggplot(df_summary, aes(x = threshold, y = avg_f_meas, color = tool, group = tool)) +
  geom_line() +
  geom_ribbon(aes(ymin = avg_f_meas - sd_f_meas, ymax = avg_f_meas + sd_f_meas, fill = tool),
              alpha = 0.5,
              colour = NA) +
  facet_wrap(~ category, scales = "free_x") +
  xlab("Relative abundance threshold") +
  ylab("Average F1 score") +
  scale_color_manual(values = classifier_colors, name = "Classifier") +
  scale_fill_manual(values = classifier_colors, name = "Classifier") +  # Add scale_fill_manual() for fill color
  theme_bw() +
  theme(
    strip.text = element_text(margin = margin(0, 0, 10, 0), face = "bold", size = 17),
    axis.text = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 17),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold")
  )
plot
plot <- plot + theme(axis.title.x = element_blank())

ggsave("Plots/F1score.png", plot, width = 14, height = 8, dpi = 300)

##########################################################################
library(cowplot)
classifier_colors <- c(
  "Kraken" = "#798234",
  "custom.Kraken" = "#70A494",
  "Kaiju" = "#BD925A",
  "MetaPhlAn3" =  "#DE8A5A",
  "MetaPhlAn4" =  "#CA562C"
)

df_summary_prec <- df1 %>%
  filter(tool != "Kraken") %>%
  filter(category != "Species") %>%
  mutate(tool = ifelse(tool == "Kraken22", "Kraken", tool)) %>%
  group_by(threshold, tool) %>%
  summarize(
    avg_prec = mean(precision),
    sd_prec = sd(precision)
  )

df_summary_spec <- df1 %>%
  filter(tool != "Kraken") %>%
  filter(category != "Species") %>%
  mutate(tool = ifelse(tool == "Kraken22", "Kraken", tool)) %>%
  group_by(threshold, tool) %>%
  summarize(
    avg_spec = mean(sens),
    sd_spec = sd(sens),
    .groups = 'drop'  # This line avoids the warning message about regrouping output
  )

# Precision plot
prec_plot <- ggplot(df_summary_prec, aes(x = threshold, y = avg_prec, color = tool, group = tool)) +
  geom_line() +
  geom_ribbon(aes(ymin = avg_prec - sd_prec, ymax = avg_prec + sd_prec, fill = tool),
              alpha = 0.5,
              colour = NA) +
  ggtitle("Species") + # Correctly setting the title here
  xlab("Relative abundance threshold") +
  ylab("Average Precision") +
  scale_color_manual(values = classifier_colors, name = "Classifier") +
  scale_fill_manual(values = classifier_colors) +
  theme_bw() +
  theme(
    strip.text = element_text(margin = margin(0, 0, 10, 0), face = "bold", size = 17),
    axis.text = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 17),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face ="bold", hjust = 0.5, size = 17),
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.text = element_text(size = 14, face = "bold"),  # Increase legend text size and make it bold
    legend.title = element_text(size = 16, face = "bold"),  # Increase legend title size and make it bold
    legend.position = "none")
prec_plot

# Sensitivity plot
spec_plot <- ggplot(df_summary_spec, aes(x = threshold, y = avg_spec, color = tool, group = tool)) +
  geom_line() +
  geom_ribbon(aes(ymin = avg_spec - sd_spec, ymax = avg_spec + sd_spec, fill = tool),
              alpha = 0.5,
              colour = NA)+
  ggtitle("Species") +
  xlab("Relative abundance threshold") +
  ylab("Average Sensitivity") +
  scale_color_manual(values = classifier_colors, name = "Classifier") +
  scale_fill_manual(values = classifier_colors) +
  theme_bw() +
  theme(
    strip.text = element_text(margin = margin(0, 0, 10, 0), face = "bold", size = 17),
    axis.text = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 17),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face ="bold", hjust = 0.5, size = 17),
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.text = element_text(size = 14, face = "bold"),  # Increase legend text size and make it bold
    legend.title = element_text(size = 16, face = "bold"),  # Increase legend title size and make it bold
    legend.position = "none")

# Combine the two plots into a single row
plot_combined <- plot_grid(prec_plot, spec_plot, align = "h")


# Get the legend
legend <- get_legend(
  ggplot(df_summary_prec, aes(x = threshold, y = avg_prec, color = tool)) +
    geom_line() +
    geom_ribbon(aes(ymin = avg_prec - sd_prec, ymax = avg_prec + sd_prec, fill = tool),
                alpha = 0.5,
                colour = NA) +
    scale_color_manual(values = classifier_colors, name = "Classifier") +
    scale_fill_manual(values = classifier_colors, name = "Classifier") +  # Add scale_fill_manual() for fill color
    theme(
      legend.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 16, face = "bold")
    )
)

# Combine the plots and the legend
final_plot <- plot_grid(plot_combined, legend, rel_widths = c(5, 1))

print(final_plot)
ggsave("Plots/prec.sens.png", final_plot, width = 14, height = 8, dpi = 600)

