#DATA CLEANING: ABUNDANCE FILES

##############################################################################################################################
library(dplyr)
library(tidyr)
#For the 7 runs bring in results and add the run name

#Read Expected abundance files
#Add the name of the Run to the dataframe
AB_NX1 <- read.csv("NX1_abundance.txt", header=F, sep="" )
AB_NX1['Run'] = '01'
AB_NX2 <- read.csv("NX2_abundance.txt", header=F, sep="" )
AB_NX2['Run'] = '02'

#Merge all dataframes together
AB <- rbind(AB_NX1,AB_NX2)

#AB <- AB_NX1
#Rename the column_names
colnames(AB)[1] <- "sample"
colnames(AB)[2] <- "accession_id"
colnames(AB)[3] <- "abundance"

########

df_species <- read.csv("df_species.csv", header=T)

AB_new <- dplyr::left_join(AB,df_species, by="accession_id")
AB_new['Classifier'] = 'insilico'

AB_new$sample_name = paste(AB_new$sample, AB_new$Run,  AB_new$Classifier, sep=".")

expected_results <- AB_new[,c(3,5,7)]

expected_results$abundance <- expected_results$abundance * 100
colnames(expected_results)[2] <- "species"

#write.csv(expected_results,"expected_results.csv")

rm(list= setdiff(ls(), "expected_results"))
###########################.  TAXONOMY CLASSIFIERS .###########################
#DATA CLEANING: Kaiju Results
KJ_NX1 <- read.delim("NX1_Kaiju.txt", header = F, sep = "\t", dec = ".")
KJ_NX2 <- read.delim("NX2_Kaiju.txt", header = F, sep = "\t", dec = ".")

KJ_NX1['Run'] = '01'
KJ_NX2['Run'] = '02'

KJ_results <- rbind(KJ_NX1,KJ_NX2)
KJ_results <- KJ_results[,-c(2,4,5)]

colnames(KJ_results)[1] <- "sample"
colnames(KJ_results)[2] <- "abundance"
colnames(KJ_results)[3] <- "species"

KJ_results['Classifier'] = 'kaiju'

KJ_results$sample_name = paste(KJ_results$sample, KJ_results$Run,  KJ_results$Classifier, sep=".")
KJ_results <- KJ_results[,-c(1,4,5)]

#Couldnt double check.. going forward ..

#Removed unclassified
KJ_results <-  subset(KJ_results, species != "unclassified")
KJ_sum <- KJ_results %>% 
  group_by(sample_name, species) %>%
  summarize(abundance= sum(abundance)) %>%
  mutate(abundance = (abundance/ sum(abundance)) *100)

#KJ_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", KJ_results$species)

ncols <- max(stringr::str_count(KJ_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

KJ_results <-
  tidyr::separate(
    data = KJ_results,
    col = species,
    sep = " ",
    into = colmn,
    remove = FALSE
  )

KJ_results$species = paste(KJ_results$`col 1`, KJ_results$`col 2`,   sep=" ")
KJ_results <- KJ_results[c(1,2,14)]

#KJ_sum$abundance[is.na(KJ_sum$abundance)] <- 0
#subset_abundance <- KJ_sum$abundance[KJ_sum$sample_name] == "sample19.01.kaiju"
#sum_abundance <-  sum(KJ_sum$abundance[KJ_sum$sample_name] == "sample19.01.kaiju")
#Couldnt double check.. going forward ..

####
#DATA CLEANING: Metaphlan Results
###############################
#SPECIES METAPHLAN DATA
MP_NX1 <- read.delim("NX1_metaphlan_species.txt", header = T, sep = "\t", dec = ".")
MP_NX2 <- read.delim("NX2_metaphlan_species.txt", header = T, sep = "\t", dec = ".")

MP_NX1 <- MP_NX1 %>% pivot_longer(-c(sample), names_to="sample_names", values_to = "abundance")
MP_NX2 <- MP_NX2 %>% pivot_longer(-c(sample), names_to="sample_names", values_to = "abundance")

MP_NX1['Run'] = '01'
MP_NX2['Run'] = '02'

MP_results <- rbind(MP_NX1,MP_NX2)
colnames(MP_results)[1] <-  "species"

MP_results <- separate(MP_results, col = "sample_names", into = c("sample", "Rest"), sep = "_", extra = "drop") 
MP_results['Classifier'] = 'metaphlan'

MP_results$sample_name = paste(MP_results$sample, MP_results$Run,  MP_results$Classifier, sep=".")
MP_results <- MP_results[,-c(2,3,5,6)]

MP_results$species <-  gsub("_"," ", MP_results$species)
#MP_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", MP_results$species)

ncols <- max(stringr::str_count(MP_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

MP_results <-
  tidyr::separate(
    data = MP_results,
    col = species,
    sep = " ",
    into = colmn,
    remove = FALSE
  )

MP_results$species = paste(MP_results$`col 1`, MP_results$`col 2`,   sep=" ")
MP_results <- MP_results[c(1,8,9)]

#UNKNOWN METAPHLAN DATA
###############################
#Removed the first line manually
MP_unknown_NX1 <- read.csv("NX1_Metaphlan.tsv", header=T, sep ="" )
MP_unknown_NX2 <- read.csv("NX2_Metaphlan.tsv", header=T, sep ="" )

MP_unknown_NX1 <- MP_unknown_NX1[,-2]
MP_unknown_NX2 <- MP_unknown_NX2[,-2]

MP_unknown_NX1 <- MP_unknown_NX1 %>% pivot_longer(-c(clade_name), names_to="sample_names", values_to = "abundance")
MP_unknown_NX2 <- MP_unknown_NX2 %>% pivot_longer(-c(clade_name), names_to="sample_names", values_to = "abundance")

MP_unknown_NX1['Run'] = 'NX1'
MP_unknown_NX2['Run'] = 'NX2'

MP_unknown_NX1 <- MP_unknown_NX1[MP_unknown_NX1$clade_name %like% "UNKNOWN", ]
MP_unknown_NX2 <- MP_unknown_NX2[MP_unknown_NX2$clade_name %like% "UNKNOWN", ]

MP_unknown <- rbind(MP_unknown_NX1,MP_unknown_NX1)
colnames(MP_unknown)[1] <- "species"

MP_results  <- rbind(MP_unknown,MP_results)
rm(MP_unknown)




#DATA CLEANING: Metaphlan4 Results
###############################
#SPECIES METAPHLAN4 DATA
MP4_NX1 <- read.delim("NX1_metaphlan4.txt", header = T, sep = "\t", dec = ".")
MP4_NX2 <- read.delim("NX2_metaphlan4.txt", header = T, sep = "\t", dec = ".")

MP4_NX1 <-  MP4_NX1[grep("s__", MP4_NX1$clade_name),]
MP4_NX2 <-  MP4_NX2[grep("s__", MP4_NX2$clade_name),]

MP4_NX1$clade_name <-  gsub(".*\\|s__","", MP4_NX1$clade_name)
MP4_NX2$clade_name <-  gsub(".*\\|s__","", MP4_NX2$clade_name)

MP4_NX1$clade_name <-  gsub("_"," ", MP4_NX1$clade_name)
MP4_NX2$clade_name <-  gsub("_"," ", MP4_NX2$clade_name)

colnames(MP4_NX1) <-  gsub("_R1_mpa_out","", colnames(MP4_NX1))
colnames(MP4_NX2) <-  gsub("_R1_mpa_out","", colnames(MP4_NX2))

MP4_NX1 <- MP4_NX1 %>% pivot_longer(-c(clade_name), names_to="sample_names", values_to = "abundance")
MP4_NX2 <- MP4_NX2 %>% pivot_longer(-c(clade_name), names_to="sample_names", values_to = "abundance")

MP4_NX1['Run'] = '01'
MP4_NX2['Run'] = '02'

MP4_results <- rbind(MP4_NX1,MP4_NX2)
colnames(MP4_results)[1] <-  "species"

MP4_results['Classifier'] = 'metaphlan4'

MP4_results$sample_name = paste(MP4_results$sample_names, MP4_results$Run,  MP4_results$Classifier, sep=".")
MP4_results <- MP4_results[,c(1,3,6)]

MP4_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", MP4_results$species)

ncols <- max(stringr::str_count(MP4_results$species, " ")) + 1
colmn <- paste("col",1:ncols)

#DATA CLEANING: Kraken Results
#########################################
KR_NX1 <- read.delim("NX1_Kraken.txt", header = F, sep = "\t", dec = ".")
KR_NX2 <- read.delim("NX2_Kraken.txt", header = F, sep = "\t", dec = ".")

KR_NX1['Run'] = '01'
KR_NX2['Run'] = '02'

KR_results <- rbind(KR_NX1,KR_NX2)
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
KR_results <- KR_results[c(1,9,10)]

#KR_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", KR_results$species)

#DATA CLEANING: Kraken_gtdb Results
KG_NX1 <- read.delim("NX1_Kraken_gtdb.txt", header = F, sep = "\t", dec = ".")
KG_NX2 <- read.delim("NX2_Kraken_gtdb.txt", header = F, sep = "\t", dec = ".")

KG_NX1['Run'] = '01'
KG_NX2['Run'] = '02'

KG_results <- rbind(KG_NX1,KG_NX2)
KG_results <- KG_results[,-c(3:7)]

colnames(KG_results)[1] <- "sample"
colnames(KG_results)[2] <- "species"
colnames(KG_results)[3] <- "abundance"

KG_results['Classifier'] = 'kraken_gtdb'

KG_results$sample_name = paste(KG_results$sample, KG_results$Run,  KG_results$Classifier, sep=".")
KG_results <- KG_results[,-c(1,4,5)]

KG_results$abundance <- KG_results$abundance * 100

#KG_results$species <- gsub("^(\\w+\\s+\\w+).*", "\\1", KG_results$species)

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

observed_results <- rbind(KG_results,KJ_results,KR_results,MP_results )
write.csv(observed_results,"observed_results_2.csv")

#Remove all files with NX in their name
rm(list=ls())
