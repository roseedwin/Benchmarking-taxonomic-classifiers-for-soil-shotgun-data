#DATA CLEANING: ABUNDANCE FILES
#EXPECTED RESULTS 
################################################################################################
library(dplyr)
library(tidyr)
#For the 7 runs bring in results and add the run name

#Read Expected abundance files
#Add the name of the Run to the dataframe

# create an empty list to store the data frames
df_list <- list()

# use a for loop to read each CSV file and add a "Run" column with the corresponding run number
for (i in 2:8) {
  filename <- paste0("merged_Run", i, "_abundance.txt")
  run_num <- sprintf("%02d", i)
  df <- read.csv(filename, header=F, sep="")
  df$Run <- run_num
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
AB <- do.call(rbind, df_list)

#Rename the column_names
colnames(AB)[1] <- "sample"
colnames(AB)[2] <- "accession_id"
colnames(AB)[3] <- "abundance"


########################################################################################################
df_species <- read.csv("df_species.csv", header=T)

AB_new <- dplyr::left_join(AB,df_species, by="accession_id")
AB_new['Classifier'] = 'insilico'

AB_new$sample_name = paste(AB_new$sample, AB_new$Run,  AB_new$Classifier, sep=".")

expected_results <- AB_new[,c(3,5,7)]

expected_results$abundance <- expected_results$abundance * 100
colnames(expected_results)[2] <- "species"

write.csv(expected_results,"expected_results.csv")
rm(list= setdiff(ls(), "expected_results"))


################################################################################################
###########################.  TAXONOMY CLASSIFIERS .############################################
################################################################################################
################################################################################################

################################################################################################
#DATA CLEANING: Kaiju Results
################################################################################################

df_list <- list()

# use a for loop to read each file, add a "Run" column with the corresponding run number, and store the data frame in the list
for (i in 2:8) {
  filename <- paste0("merged_Run", i, "_Kaiju.txt")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = F, sep = "\t", dec = ".")
  df$Run <- run_num
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
KJ_results <- do.call(rbind, df_list)

KJ_results <- KJ_results[,-c(2,4,5)]

colnames(KJ_results)[1] <- "sample"
colnames(KJ_results)[2] <- "abundance"
colnames(KJ_results)[3] <- "species"

KJ_results['Classifier'] = 'kaiju'

KJ_results$sample_name = paste(KJ_results$sample, KJ_results$Run,  KJ_results$Classifier, sep=".")
KJ_results <- KJ_results[,-c(1,4,5)]

###unclassified calculations ###
#unclassified <- subset(KJ_results, species == 'unclassified')
#avg_unclassified <-  mean(unclassified$abundance)

#Removed unclassified
KJ_results <-  subset(KJ_results, species != "unclassified")
KJ_sum <- KJ_results %>% 
  group_by(sample_name, species) %>%
  summarize(abundance= sum(abundance)) %>%
  mutate(abundance = (abundance/ sum(abundance)) *100)

#Removed "cannot be assigned to a (non-viral) species"
KJ_results <-  subset(KJ_results, species != "cannot be assigned to a (non-viral) species")
KJ_sum <- KJ_results %>% 
  group_by(sample_name, species) %>%
  summarize(abundance= sum(abundance)) %>%
  mutate(abundance = (abundance/ sum(abundance)) *100)

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

################################################################################################
#DATA CLEANING: Metaphlan4 Results
################################################################################################

# create an empty list to store the data frames
df_list <- list()

# use a for loop to read each file, preprocess it, add a "Run" column with the corresponding run number, and store the data frame in the list
for (i in 2:8) {
  filename <- paste0("merged_Run", i, "_metaphlan4.txt")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = T, sep = "\t", dec = ".")
  df <- df[grep("s__", df$clade_name),]
  df$clade_name <- gsub(".*\\|s__","", df$clade_name)
  df$clade_name <- gsub("_"," ", df$clade_name)
  colnames(df) <- gsub("_R1_mpa_out","", colnames(df))
  df <- df %>% pivot_longer(-c(clade_name), names_to="sample_names", values_to = "abundance")
  df$Run <- run_num
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
MP4_results <- do.call(rbind, df_list)

# print the combined data frame
MP4_results

colnames(MP4_results)[1] <-  "species"

MP4_results['Classifier'] = 'metaphlan4'

MP4_results$sample_name = paste(MP4_results$sample_names, MP4_results$Run,  MP4_results$Classifier, sep=".")
MP4_results <- MP4_results[,c(1,3,6)]

MP4_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", MP4_results$species)

ncols <- max(stringr::str_count(MP4_results$species, " ")) + 1
colmn <- paste("col",1:ncols)


################################################################################################
#DATA CLEANING: Kraken Results
################################################################################################

# create an empty list to store the data frames
df_list <- list()

# use a for loop to read each file, add a "Run" column with the corresponding run number, and store the data frame in the list
for (i in 2:8) {
  filename <- paste0("merged_Run", i, "_Kraken.txt")
  run_num <- sprintf("%02d", i)
  df <- read.delim(filename, header = F, sep = "\t", dec = ".")
  df$Run <- run_num
  df_list[[i]] <- df
}

# combine all the data frames into a single data frame using rbind
KR_results <- do.call(rbind, df_list)

# print the combined data frame
KR_results

KR_results <- KR_results[,-c(3:7)]

colnames(KR_results)[1] <- "sample"
colnames(KR_results)[2] <- "species"
colnames(KR_results)[3] <- "abundance"

KR_results['Classifier'] = 'kraken'

KR_results$sample_name = paste(KR_results$sample, KR_results$Run,  KR_results$Classifier, sep=".")
KR_results <- KR_results[,-c(1,4,5)]

KR_results$abundance <- KR_results$abundance * 100

ncols <- max(stringr::str_count(KR_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

KR_results <-
  tidyr::separate(
    data = KR_results,
    col = species,
    sep = " ",
    into = colmn,
    remove = FALSE
  )

KR_results$species = paste(KR_results$`col 1`, KR_results$`col 2`,   sep=" ")
KR_results <- KR_results[c(1,10,11)]

#Check sum
#sum_abundance <- aggregate(KR_results$abundance, by =list(KR_results$sample_name), FUN = sum)
#colnames(sum_abundance)  <- c("sample_name", "total_abundance")
#sum_abundance


################################################################################################
#DATA CLEANING: Kraken_gtdb Results
################################################################################################

# create an empty list to store the data frames
df_list <- list()

# use a for loop to read each file, add a "Run" column with the corresponding run number, and store the data frame in the list
for (i in 2:8) {
  filename <- paste0("merged_Run", i, "_Kraken_gtdb.txt")
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

KG_results['Classifier'] = 'kraken_gtdb'

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

observed_results <- rbind(KG_results, KJ_results, KR_results, MP4_results )
write.csv(observed_results,"observed_results.csv")
