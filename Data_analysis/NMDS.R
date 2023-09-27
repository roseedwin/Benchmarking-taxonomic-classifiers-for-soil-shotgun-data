#Prepping data for NMDS
library(dplyr)
library(tidyr)

df_long <-  read.table("Output/Species_results.tsv", header =T,  sep="\t")
df_long <-  df_long %>% separate(sample_name, into = c("sample_name", "run", "tool"), sep = "\\.")
df_long2 <-   df_long %>%  group_by(species, sample_name, run, tool) %>%
  summarise(rel.abundance = sum(abundance))

df_long2$tool <- gsub("Kraken_gtdb", "Kraken.gtdb", df_long2$tool)

observed <-  df_long2 %>% filter(tool != "insilico")
observed1 <- observed %>%
  filter(tool %in% c("Kaiju", "Kraken.gtdb", "Metaphlan3", "Metaphlan4")) %>%
  filter(rel.abundance > 0.001) %>%
  rename(observed = rel.abundance)

observed2 <- observed %>%
  filter(tool %in% c("Kraken", "Kraken22")) %>%
  filter_at(vars(5), any_vars(.> 0.005)) %>%
  rename(observed = rel.abundance)

observed <- rbind(observed1, observed2)

observed <-  observed %>%
  select(species, tool, observed, sample_name, run) %>%
  pivot_wider(names_from = tool, values_from = observed) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0)))

expected <-  df_long2 %>% filter(tool == "insilico") %>%
  rename(expected = rel.abundance) %>%
  filter(expected != 0)

observed_expected <- full_join( expected, observed)

result <- observed_expected %>%
  mutate(expected = ifelse(is.na(expected), 0, expected)) %>%
  mutate(tool = ifelse(is.na(tool), "insilico", tool))  %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = 5:11, names_to = "classifier", values_to = "abundance")

result2 <- result %>%
  select(-4) %>%
  mutate(classifier = gsub("sample", "", classifier)) %>%
  unite(lib_id,sample_name, run,  classifier, sep = "_")

pivot_df <- result2 %>%
  pivot_wider(names_from = lib_id, values_from = abundance)

# Convert to matrix
mat <- as.matrix(pivot_df[, -1])
# Set row names
rownames(mat) <- pivot_df$species
mat[is.na(mat)] <- 0
write.csv(mat,"Output/mat_species.csv")

#######       NMDS

#load necessary package
library(vegan)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
data <- read.csv("Output/mat_species.csv", row.names = 1)
data <- t(data)
map <- read.csv("map2.csv", row.names =1)

meta.mds.sl <- metaMDS(data, distance = "bray")
plot(meta.mds.sl$points)
#write.csv(meta.mds.sl$points, "Output/mds_spec.csv")

df <- read.csv("Output/mds_spec.csv", header = TRUE)
df2 <- df[, -1]
rownames(df2) <- df[, 1]
df3 <- df[, 1]
df <- df2
df <- merge(df, map, by = 'row.names')

run <- df[, 4]
sample <- df[, 5]
tool <- df[, 6]

df <- df %>%
  filter(tool != "Kraken") %>%
  mutate(tool = ifelse(tool == "Kraken22", "Kraken", tool)) %>%
  mutate(tool = ifelse(tool == "Kraken.gtdb", "custom.Kraken", tool))
df$tool <- gsub("Metaphlan", "MetaPhlAn", df$tool)

classifier_colors <- c(
  "Kraken" = "#798234",
  "custom.Kraken" = "#70A494",
  "Kaiju" = "#BD925A",
  "MetaPhlAn3" =  "#DE8A5A",
  "MetaPhlAn4" =  "#CA562C",
  "expected" = "#008080"
)

### Plot the distance matrix with ggplot2 ###
Sp <- ggplot(df, aes(MDS1, MDS2, label = tool)) +
  theme_bw() +
  theme(
    legend.title = element_text(face = "bold", size = 12.5),
    axis.title.x = element_text(face = "bold", size = 12.5),
    axis.title.y = element_text(face = "bold", size = 12.5),
    plot.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 12.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12.5),
    axis.title = element_text(size = 12.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(panel.background = element_rect(colour = "black")) +
  ggtitle("Bray-Curtis MDS plot (coloured by classifier)") +
  geom_point(aes(color = tool), size = 1.5, shape = 21, fill = "white", stroke = 0.8) +
  stat_ellipse(aes(color = tool, fill = tool, alpha = 0.5), show.legend = FALSE) +
  scale_color_manual(values = classifier_colors, name = "Classifier") +
  scale_fill_manual(values = classifier_colors, name = "Classifier") +
  guides(color = guide_legend(
    override.aes = list(shape = 19, fill = "white", size = 3)),
    fill = guide_legend(
      override.aes = list(shape = 19, fill = "white", size = 3)))

library(tibble)
# Move the rownames to a column called "rowname"
data <- data %>%
  as.data.frame() %>%
  rownames_to_column(var = "rowname")
# Filter out rows where "rowname" ends with "_Kraken"
data <- data %>%
  filter(!grepl("_Kraken$", rowname))
# Convert data back to a matrix
mat1 <- as.matrix(data[, -1])  # Exclude the first column which contains the old rownames
# Set the rownames of the new matrix
rownames(mat1) <- data[, "rowname"]  # Use the old rownames stored in the first column of data
# Calculate dissimilarity matrix using Bray-Curtis distance
dissimilarity <- vegdist(mat1, method = "bray")
# Perform PERMANOVA
permanova_result <- adonis2(dissimilarity ~ tool, df, permutations = 999)

# Print PERMANOVA results
print(permanova_result)
# Extract R-squared and p-value for "tool" term
r2_tool <- permanova_result$R2[1]
p_value_tool <- permanova_result$Pr[1]

# Format R-squared as a string with two decimal places
r2_label <- paste0("PERMANOVA: R", "\u00B2", " = ", format(round(r2_tool, 2), nsmall = 2))

# Format p-value with significance stars
p_value_label <- ifelse(p_value_tool < 0.001, "p < 0.001 ***",
                        ifelse(p_value_tool < 0.01, paste0("p = ", format(p_value_tool, nsmall = 3), " **"),
                               ifelse(p_value_tool < 0.05, paste0("p = ", format(p_value_tool, nsmall = 2), " *"),
                                      paste0("p = ", format(p_value_tool, nsmall = 2)))))


# Print the labels
cat(r2_label, "\n")
cat(p_value_label, "\n")

# Update the plot
Sp <- Sp +
  geom_text(
    aes(label = r2_label),
    x = max(df$MDS1), y = min(df$MDS2),
    hjust = 0.9, vjust = 0,
    size = 4,  color = "black"
  ) +
  geom_text(
    aes(label = p_value_label),
    x = max(df$MDS1), y = min(df$MDS2),
    hjust = 1, vjust = 1.5,
    size = 4,  color = "black"
  )
print(Sp)

ggsave("Plots/NMDS.png", Sp, width = 10, height = 6, dpi = 600)






