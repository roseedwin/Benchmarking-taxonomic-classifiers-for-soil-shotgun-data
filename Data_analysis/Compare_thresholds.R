##################################################################################
######### SUBSETTING FOR THRESHOLDS: CHECKING SENSITIVITY, PRECISION, BALANCED ACCURACY AND F1 SCORE AT DIFFERENT THRESHOLDS
##################################################################################

library(tidyr)
library(ggplot2)
library(reshape2)

observed_results <- read.csv("observed_results_2.csv", header=T)
observed_results <- observed_results[,c(2,3,4)]
observed_results$species <- gsub('[^[:alnum:] ]','',observed_results$species)

expected_results <- read.csv("expected_results.csv", header=T)
expected_results <- expected_results[,c(2,3,4)]
expected_results$species <- gsub('[^[:alnum:] ]','',expected_results$species)

data <- rbind(observed_results,expected_results)

data_0.0001 <- data %>% filter_at(vars(2), any_vars(. > 0.0001))
data_0.0005 <- data %>% filter_at(vars(2), any_vars(. > 0.0005))
data_0.001 <- data %>% filter_at(vars(2), any_vars(. > 0.001))
data_0.005 <- data %>% filter_at(vars(2), any_vars(. > 0.005))
data_0.01 <- data %>% filter_at(vars(2), any_vars(. > 0.01))
  
rm(expected_results, observed_results)

####################################################################################

#change data here for the different data frames
c<-dcast(data_0.005,species~sample_name, value.var= 'abundance', fun.aggregate = sum)

######################################
c[is.na(c)] <- 0

samp2 <- c[,-1]
rownames(samp2) <- c[,1]
samp2[samp2 > 0] <- 1

# Create a list of column names for each series
kaiju_cols <- paste0("sample", 1:20, ".02.kaiju")
kraken_cols <- paste0("sample", 1:20, ".02.kraken")
kraken_gtdb_cols <- paste0("sample", 1:20, ".02.kraken_gtdb")
metaphlan_cols <- paste0("sample", 1:20, ".02.metaphlan")
metaphlan4_cols <- paste0("sample", 1:20, ".02.metaphlan4")
insilico_cols <- paste0("sample", 1:20, ".02.insilico")

##################################### 
# Loop over the columns and update new columns for each sample and each series
for (i in 1:20) {
  # Kaiju vs insilico
  kaiju_name <- paste0("sample", i, "_kaiju")
  samp2[[paste0(kaiju_name, "_result")]] <- "FP"
  samp2[[paste0(kaiju_name, "_result")]][samp2[, kaiju_cols[i]] == 1 & samp2[, insilico_cols[i]] == 1] <- "TP"
  samp2[[paste0(kaiju_name, "_result")]][samp2[, kaiju_cols[i]] == 0 & samp2[, insilico_cols[i]] == 1] <- "FN"
  samp2[[paste0(kaiju_name, "_result")]][samp2[, kaiju_cols[i]] == 1 & samp2[, insilico_cols[i]] == 0] <- "FP"
}

for (i in 1:20) {
  # Kraken vs insilico
  kraken_name <- paste0("sample", i, "_kraken")
  samp2[[paste0(kraken_name, "_result")]] <- "FP"
  samp2[[paste0(kraken_name, "_result")]][samp2[, kraken_cols[i]] == 1 & samp2[, insilico_cols[i]] == 1] <- "TP"
  samp2[[paste0(kraken_name, "_result")]][samp2[, kraken_cols[i]] == 0 & samp2[, insilico_cols[i]] == 1] <- "FN"
  samp2[[paste0(kraken_name, "_result")]][samp2[, kraken_cols[i]] == 1 & samp2[, insilico_cols[i]] == 0] <- "FP"
}
for (i in 1:20) {
  # Kraken_gtdb vs insilico
  kraken_gtdb_name <- paste0("sample", i, "_kraken.gtdb")
  samp2[[paste0(kraken_gtdb_name, "_result")]] <- "FP"
  samp2[[paste0(kraken_gtdb_name, "_result")]][samp2[, kraken_gtdb_cols[i]] == 1 & samp2[, insilico_cols[i]] == 1] <- "TP"
  samp2[[paste0(kraken_gtdb_name, "_result")]][samp2[, kraken_gtdb_cols[i]] == 0 & samp2[, insilico_cols[i]] == 1] <- "FN"
  samp2[[paste0(kraken_gtdb_name, "_result")]][samp2[, kraken_gtdb_cols[i]] == 1 & samp2[, insilico_cols[i]] == 0] <- "FP"
}
for (i in 1:20) {
  # metaphlan vs insilico
  metaphlan_name <- paste0("sample", i, "_metaphlan")
  samp2[[paste0(metaphlan_name, "_result")]] <- "FP"
  samp2[[paste0(metaphlan_name, "_result")]][samp2[, metaphlan_cols[i]] == 1 & samp2[, insilico_cols[i]] == 1] <- "TP"
  samp2[[paste0(metaphlan_name, "_result")]][samp2[, metaphlan_cols[i]] == 0 & samp2[, insilico_cols[i]] == 1] <- "FN"
  samp2[[paste0(metaphlan_name, "_result")]][samp2[, metaphlan_cols[i]] == 1 & samp2[, insilico_cols[i]] == 0] <- "FP"
}

for (i in 1:20) {
   metaphlan4 vs insilico
  metaphlan4_name <- paste0("sample", i, "_metaphlan4")
  samp2[[paste0(metaphlan4_name, "_result")]][samp2[, metaphlan4_cols[i]] == 1 & samp2[, insilico_cols[i]] == 1] <- "TP"
  samp2[[paste0(metaphlan4_name, "_result")]] <- "FP"
  samp2[[paste0(metaphlan4_name, "_result")]][samp2[, metaphlan4_cols[i]] == 0 & samp2[, insilico_cols[i]] == 1] <- "FN"
  samp2[[paste0(metaphlan4_name, "_result")]][samp2[, metaphlan4_cols[i]] == 1 & samp2[, insilico_cols[i]] == 0] <- "FP"
}

# View the updated data frame
samp2

result_cols <- grepl("_result", colnames(samp2))
result_data <- samp2[, result_cols]

# Create a list of all the result column names for kaiju, kraken, and kraken_gtdb
kaiju_cols <- grep("_kaiju_result", colnames(result_data), value = TRUE)
kraken_cols <- grep("_kraken_result", colnames(result_data), value = TRUE)
kraken_gtdb_cols <- grep("_kraken.gtdb_result", colnames(result_data), value = TRUE)
metaphlan_cols <- grep("_metaphlan_result", colnames(result_data), value = TRUE)
metaphlan4_cols <- grep("_metaphlan4_result", colnames(result_data), value = TRUE)
result_cols <- c(kaiju_cols, kraken_cols, kraken_gtdb_cols, metaphlan_cols, metaphlan4_cols)


####################################################################################

################# SENSITIVITY #################
# Create an empty data frame to store sensitivity values
sensitivity_df <- data.frame(Column = character(length(result_cols)), Sensitivity = numeric(length(result_cols)))

# Loop over the result columns and calculate sensitivity
for (i in 1:length(result_cols)) {
  col <- result_cols[i]
  
  # Extract the sample and tool names from the column name
  sample_tool <- gsub("_result", "", col)
  sample_tool <- gsub("_result", "", sample_tool)
  
  # Extract the TP and FN values for the current column
  TP <- sum(result_data[, col] == "TP")
  FN <- sum(result_data[, col] == "FN")
  
  # Calculate sensitivity for the current column
  sensitivity <- TP / (TP + FN)
  
  # Add the sensitivity value to the sensitivity data frame
  sensitivity_df[i, "Column"] <- col
  sensitivity_df[i, "Sensitivity"] <- sensitivity
}

# Remove any rows with missing sensitivity values
sensitivity_df <- sensitivity_df[complete.cases(sensitivity_df), ]


################# PRECISION #################
# Create an empty data frame to store precision values
precision_df <- data.frame(Column = character(length(result_cols)), precision = numeric(length(result_cols)))

# Loop over the result columns and calculate precision
for (i in 1:length(result_cols)) {
  col <- result_cols[i]
  
  # Extract the sample and tool names from the column name
  sample_tool <- gsub("_result", "", col)
  sample_tool <- gsub("_result", "", sample_tool)
  
  # Extract the TP and FN values for the current column
  TP <- sum(result_data[, col] == "TP")
  FP <- sum(result_data[, col] == "FP")
  
  # Calculate precision for the current column
  precision <- TP / (TP + FP)
  
  # Add the precision value to the precision data frame
  precision_df[i, "Column"] <- col
  precision_df[i, "precision"] <- precision
}

# Remove any rows with missing precision values
precision_df <- precision_df[complete.cases(precision_df), ]


####################################################################################
#repeat the same for all thresholds
df_0.005 <- merge(sensitivity_df,precision_df, by="Column")
df_0.005$Threshold <- "0.005"
df_0.005$Run <- "002"
####################################################################################

df <- rbind(df_0.01,df_0.001,df_0.005,df_0.0001, df_0.0005, df_all)
df$balanced_accuracy <- (df$Sensitivity + df$precision) /2
df$F1_score <- 2 * (df$Sensitivity * df$precision) / (df$Sensitivity + df$precision)

#Repeated the same for other runs
write.csv(df, "Classification_accuracy_Run2_diff_thresholds.csv")
