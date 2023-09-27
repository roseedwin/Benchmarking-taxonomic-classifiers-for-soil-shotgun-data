#CLEAN ALL THE FILES AND CREATE TSV FILES FOR SPECIES LEVEL
#SAVE ALL THE FILES IN FOLDER OUTPUT_FILES
################################################################################################
library(dplyr)
library(tidyr)
library(reshape2)
################################################################################################
# DATA CLEANING: ABUNDANCE PROFILES
################################################################################################
# create an empty list to store the data frames
df_list <- list()

# use a for loop to read each CSV file and add a "Run" column with the corresponding run number
for (i in 1:10) {
  filename <- paste0("data/Run", i, "_abundance.txt")
  df <- read.csv(filename, header = FALSE, sep = "")
  df$Run <- paste0("Run", i)
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
AB <- do.call(rbind, df_list)

#Rename the column_names
colnames(AB)[1] <- "sample"
colnames(AB)[2] <- "accessionversion"
colnames(AB)[3] <- "abundance"

df_taxonomy <-  read.csv("df_taxonomy2.csv", header =T)

#Seperate taxonomy file with all the data relating to levels of taxonomy classification and matching accession IDs was made separately
AB_new <- dplyr::left_join(AB,df_taxonomy, by="accessionversion")
AB_new['Classifier'] = 'insilico'

#df_na <- AB_new [apply(is.na(AB_new ), 1, any), ]
#unique_accession_df_na <- unique(df_na$accessionversion)
AB_new[is.na(AB_new)] <- "NULL"

AB_new$sample_name = paste(AB_new$sample, AB_new$Classifier, sep=".")

expected_results <- AB_new
expected_results$abundance <- expected_results$abundance * 100

write.csv(expected_results,"expected_results.csv", row.names = FALSE)
expected_results <-  read.csv("expected_results.csv", header =T)

################################################################################################
###########################.  TAXONOMY CLASSIFIERS .############################################

################################################################################################
#DATA CLEANING: Kaiju Results (species level)
################################################################################################
df_list <- list()

# use a for loop to read each file, add a "Run" column with the corresponding run number, and store the data frame in the list
for (i in 1:10) {
  filename <- paste0("data/species_Run", i, "_Kaiju.txt")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = FALSE, sep = "\t", dec = ".")
  df$Run <- run_num
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
KJ_results <- do.call(rbind, df_list)

KJ_results <- KJ_results[,-c(2,4,5)]

colnames(KJ_results)[1] <- "sample"
colnames(KJ_results)[2] <- "abundance"
colnames(KJ_results)[3] <- "species"

KJ_results['Classifier'] = 'Kaiju'

KJ_results$Run <- gsub("^0", "", KJ_results$Run)
#unique(KJ_results$Run)
KJ_results$sample_name = paste(KJ_results$sample, KJ_results$Run,  KJ_results$Classifier, sep=".")
KJ_results <- KJ_results[,-c(1,4,5)]

###unclassified calculations ###
unclassified <- subset(KJ_results, species == 'unclassified')
avg_unclassified <-  mean(unclassified$abundance)

#Removed unclassified
#KJ_results <-  subset(KJ_results, species != "unclassified")
#KJ_results <- KJ_results %>%
#  group_by(sample_name, species) %>%
#  summarize(abundance= sum(abundance)) %>%
#  mutate(abundance = (abundance/ sum(abundance)) *100)

#Removed "cannot be assigned to a (non-viral) species"
KJ_results$species <- sub("cannot be assigned to a \\(non-viral\\) species", "cannot_be_assigned_to_nonviral_species", KJ_results$species)
#KJ_results <-  subset(KJ_results, species != "cannot be assigned to a (non-viral) species")
#KJ_results <- KJ_results %>%
#  group_by(sample_name, species) %>%
#  summarize(abundance= sum(abundance)) %>%
#  mutate(abundance = (abundance/ sum(abundance)) *100)

ncols <- max(stringr::str_count(KJ_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

KJ_results <- tidyr::separate(
  data = KJ_results,
  col = species,
  sep = " ",
  into = colmn,
  remove = FALSE
)

KJ_results$species = paste(KJ_results$`col 1`, KJ_results$`col 2`,   sep=" ")
KJ_results <- KJ_results[c(1,2,14)]

KJ_results$species <- gsub("[^A-Za-z0-9 ]", "", KJ_results$species)

#Check sum
sum_abundance <- aggregate(KJ_results$abundance, by =list(KJ_results$sample_name), FUN = sum)
colnames(sum_abundance)  <- c("sample_name", "total_abundance")
sum_abundance

KJ_results$species <- sub("unclassified NA", "unclassified", KJ_results$species)
KJ_results$species <- sub("cannotbeassignedtononviralspecies NA", "cannotbeassignedtononviralspecies", KJ_results$species)
KJ_results$species <- gsub("[^[:alnum:][:blank:]]", "", KJ_results$species)

################################################################################################
#DATA CLEANING: Metaphlan4 Results (species level)
################################################################################################

# Create an empty list to store the data frames
df_list <- list()

# Use a for loop to read each TSV file and perform data manipulations
for (i in 1:10) {
  filename <- paste0("data/Run", i, "_Metaphlan4.tsv")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = TRUE, sep = "\t", dec = ".")
  df_classified <- df[grep("s__", df$clade_name), ]
  df_unclassified <- df[grep("UNCLASSIFIED", df$clade_name), ]

  df_classified <- filter(df_classified, !grepl("\\|t__", df_classified$clade_name))
  df_classified$clade_name <- gsub(".*\\|s__", "", df_classified$clade_name)
  df_classified$clade_name <- gsub("_", " ", df_classified$clade_name)
  colnames(df_classified) <- gsub("_profile", "", colnames(df_classified))
  df_classified <- df_classified[, -2]
  df_classified <- pivot_longer(df_classified, cols = -c(clade_name), names_to = "sample_names", values_to = "abundance")
  df_classified$Run <- run_num

  df_unclassified$clade_name <- gsub("UNCLASSIFIED", "unclassified", df_unclassified$clade_name)
  colnames(df_unclassified) <- gsub("_profile", "", colnames(df_unclassified))
  df_unclassified <- df_unclassified[, -2]
  df_unclassified <- pivot_longer(df_unclassified, cols = -c(clade_name), names_to = "sample_names", values_to = "abundance")
  df_unclassified$Run <- run_num

  df_list[[i]] <- rbind(df_classified, df_unclassified)
}

# combine all the data frames into a single data frame using rbind
MP4_results <- do.call(rbind, df_list)

#avg_unclassified <-  mean(MP4_unclassified$abundance)

colnames(MP4_results)[1] <-  "species"

MP4_results['Classifier'] = 'Metaphlan4'

MP4_results$Run <- gsub("^0", "", MP4_results$Run)
#unique(MP4_results$Run)
MP4_results$sample_name = paste(MP4_results$sample_names, MP4_results$Run,  MP4_results$Classifier, sep=".")
MP4_results <- MP4_results[,c(1,3,6)]

MP4_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", MP4_results$species)

ncols <- max(stringr::str_count(MP4_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

#Check sum
sum_abundance <- aggregate(MP4_results$abundance, by =list(MP4_results$sample_name), FUN = sum)
colnames(sum_abundance)  <- c("sample_name", "total_abundance")
sum_abundance
MP4_results$species <- gsub("[^[:alnum:][:blank:]]", "", MP4_results$species)
################################################################################################
#DATA CLEANING: Metaphlan3 Results (species level)
################################################################################################

# Create an empty list to store the data frames
df_list <- list()

# Use a for loop to read each TSV file and perform data manipulations
for (i in 1:10) {
  filename <- paste0("data/Run", i, "_Metaphlan3.tsv")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = TRUE, sep = "\t", dec = ".", skip =1)
  df_classified <- df[grep("s__", df$clade_name), ]
  df_unclassified <- df[grep("UNKNOWN", df$clade_name), ]

  df_classified <- filter(df_classified, !grepl("\\|t__", df_classified$clade_name))
  df_classified$clade_name <- gsub(".*\\|s__", "", df_classified$clade_name)
  df_classified$clade_name <- gsub("_", " ", df_classified$clade_name)
  colnames(df_classified) <- gsub("_profile", "", colnames(df_classified))
  df_classified <- df_classified[, -2]
  df_classified <- pivot_longer(df_classified, cols = -c(clade_name), names_to = "sample_names", values_to = "abundance")
  df_classified$Run <- run_num

  df_unclassified$clade_name <- gsub("UNKNOWN", "unclassified", df_unclassified$clade_name)
  colnames(df_unclassified) <- gsub("_profile", "", colnames(df_unclassified))
  df_unclassified <- df_unclassified[, -2]
  df_unclassified <- pivot_longer(df_unclassified, cols = -c(clade_name), names_to = "sample_names", values_to = "abundance")
  df_unclassified$Run <- run_num

  df_list[[i]] <- rbind(df_classified, df_unclassified)
}

# combine all the data frames into a single data frame using rbind
MP3_results <- do.call(rbind, df_list)

#Unclassified
for (i in 1:10)  {
  filename <- paste0("data/Run", i, "_Metaphlan3.tsv")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = T, sep = "\t", dec = "." ,skip = 1)
  df <- df[grep("UNKNOWN", df$clade_name),]
  colnames(df) <- gsub("_profile","", colnames(df))
  df <- df[, -2]
  df <- df %>% pivot_longer(-c(clade_name), names_to="sample_names", values_to = "abundance")
  df$Run <- run_num
  df_list[[i]] <- df
}
#MP3_unclassified <- do.call(rbind, df_list)
#avg_unclassified <-  mean(MP3_unclassified$abundance)

colnames(MP3_results)[1] <-  "species"

MP3_results['Classifier'] = 'Metaphlan3'

MP3_results$Run <- gsub("^0", "", MP3_results$Run)
unique(MP3_results$Run)
MP3_results$sample_name = paste(MP3_results$sample_names, MP3_results$Run,  MP3_results$Classifier, sep=".")
MP3_results <- MP3_results[,c(1,3,6)]

MP3_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", MP3_results$species)

ncols <- max(stringr::str_count(MP3_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

#Check sum
sum_abundance <- aggregate(MP3_results$abundance, by =list(MP3_results$sample_name), FUN = sum)
colnames(sum_abundance)  <- c("sample_name", "total_abundance")
sum_abundance
MP3_results$species <- gsub("[^[:alnum:][:blank:]]", "", MP3_results$species)

################################################################################################
#DATA CLEANING: Kraken_gtdb Results (species level)
################################################################################################

# create an empty list to store the data frames
df_list <- list()

# use a for loop to read each file, add a "Run" column with the corresponding run number, and store the data frame in the list
for (i in 1:10){
  filename <- paste0("data/species_Run", i, "_Kraken_gtdb.txt")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = F, sep = "\t", dec = ".")
  df$Run <- run_num
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
KG_results <- do.call(rbind, df_list)

# print the combined data frame
KG_results

KG_results <- KG_results[,-c(3:7)]

colnames(KG_results)[1] <- "sample"
colnames(KG_results)[2] <- "species"
colnames(KG_results)[3] <- "abundance"

KG_results['Classifier'] = 'Kraken_gtdb'

KG_results$Run <- gsub("^0", "", KG_results$Run)
unique(KG_results$Run)
KG_results$sample_name = paste(KG_results$sample, KG_results$Run,  KG_results$Classifier, sep=".")
KG_results <- KG_results[,-c(1,4,5)]

KG_results$abundance <- KG_results$abundance * 100

ncols <- max(stringr::str_count(KG_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

KG_results <-
  tidyr::separate(
    data = KG_results,
    col = species,
    sep = " ",
    into = colmn,
    remove = FALSE
  )

KG_results$species = paste(KG_results$`col 1`, KG_results$`col 2`,   sep=" ")
KG_results <- KG_results[c(1,13,14)]

KG_results$species <- gsub("[^[:alnum:][:blank:]]", "", KG_results$species)

################################################################################################
#DATA CLEANING: Kraken Results (species level) 
################################################################################################

# create an empty list to store the data frames
df_list <- list()

# use a for loop to read each file, add a "Run" column with the corresponding run number, and store the data frame in the list
for (i in 1:10) {
  filename <- paste0("data/species_Run", i, "_Kraken22.txt")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = F, sep = "\t", dec = ".")
  df$Run <- run_num
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
KR22_results <- do.call(rbind, df_list)

# print the combined data frame
KR22_results

KR22_results <- KR22_results[,-c(3:7)]

colnames(KR22_results)[1] <- "sample"
colnames(KR22_results)[2] <- "species"
colnames(KR22_results)[3] <- "abundance"

KR22_results['Classifier'] = 'Kraken22'

KR22_results$Run <- gsub("^0", "", KR22_results$Run)
unique(KR22_results$Run)
KR22_results$sample_name = paste(KR22_results$sample, KR22_results$Run,  KR22_results$Classifier, sep=".")
KR22_results <- KR22_results[,-c(1,4,5)]

KR22_results$abundance <- KR22_results$abundance * 100

ncols <- max(stringr::str_count(KR22_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

KR22_results <-
  tidyr::separate(
    data = KR22_results,
    col = species,
    sep = " ",
    into = colmn,
    remove = FALSE
  )

KR22_results$species = paste(KR22_results$`col 1`, KR22_results$`col 2`,   sep=" ")
KR22_results <- KR22_results[c(1,12,13)]

KR22_results$species <- gsub("[^[:alnum:][:blank:]]", "", KR22_results$species)

#Check sum
#sum_abundance <- aggregate(KR_results$abundance, by =list(KR_results$sample_name), FUN = sum)
#colnames(sum_abundance)  <- c("sample_name", "total_abundance")
#sum_abundance

observed_results <- rbind(KG_results, KJ_results, KR_results, MP4_results, MP3_results, KR22_results )

################################################################################################
#SPECIES LEVEL FINAL : MERGED RESULTS
################################################################################################

Species_results <- observed_results

ncols <- max(stringr::str_count(Species_results$species, " ")) + 1
colmn <- paste("col",1:ncols)
Species_results <- tidyr::separate(
  data = Species_results,
  col = species,
  sep = " ",
  into = colmn,
  remove = FALSE
)
Species_results$species = paste(Species_results$`col 1`, Species_results$`col 2`,   sep=" ")
Species_results <- Species_results[c(1,4,5)]
# Remove special characters from the "species" column
Species_results$species <- gsub("[^A-Za-z0-9 ]", "", Species_results$species)

Species_results$species <- sub("unclassified NA", "unclassified", Species_results$species)
Species_results$species <- sub("cannotbeassignedtononviralspecies NA", "cannotbeassignedtononviralspecies", Species_results$species)

Exp_species <- expected_results[,c(3,4,11,13)]

Exp_species <- separate(Exp_species, sample_name, into = c("column1", "column2"), sep = "\\.")
Exp_species$Run <- gsub("Run", "", Exp_species$Run)
Exp_species$sample_name = paste(Exp_species$column1, Exp_species$Run, Exp_species$column2, sep=".")
Exp_species <- Exp_species[,c(1,3,6)]

ncols <- max(stringr::str_count(Exp_species$species, " ")) + 1
colmn <- paste("col",1:ncols)
Exp_species <- tidyr::separate(
  data = Exp_species,
  col = species,
  sep = " ",
  into = colmn,
  remove = FALSE
)
Exp_species$species = paste(Exp_species$`col 1`, Exp_species$`col 2`,   sep=" ")
Exp_species$species <- gsub("[^A-Za-z0-9 ]", "", Exp_species$species)
Exp_species <- Exp_species[,c(1,2,10)]
c <- dcast(Exp_species, species ~ sample_name, value.var = 'abundance', fun.aggregate = sum)
Exp_species <- pivot_longer(c, cols = -species, names_to = "sample_name", values_to = "abundance")

#Check sum
sum_abundance <- aggregate(Exp_species$abundance, by =list(Exp_species$sample_name), FUN = sum)
colnames(sum_abundance)  <- c("sample_name", "total_abundance")
sum_abundance #should be 100

Species_results <- rbind(Exp_species, Species_results)
write.table(Species_results, "Output/Species_results.tsv", sep="\t", quote = FALSE, row.names = F)
