#Create an empty dataframe to store the combined data frame
combined_df <- data.frame()

#use a for loop to read each CSV file and store it as separate data frame
for (i in 1:10) {
  filename <- paste0("Classification_accuracy_Run", i ,"_diff_thresholds.csv")
  df <- read.csv(filename)
  combined_df <- rbind(combined_df, df)
}

#print the combined data frame
combined_df

combined_df <- combined_df %>%
  separate(Column, c("Sample","Classifier"), sep= "_")

#filter classifiers of no interest
#Example
df_new <-  combined_df[combined_df$Classifier !="metaphlan",]
df_new <-  df_new[df_new$Classifier !="metaphlan4",]

df_new$Threshold <- as.character(df_new$Threshold)
df_new$Threshold <- gsub("1e-04", "0.0001", df_new$Threshold)
df_new$Threshold <- gsub("5e-04", "0.0005", df_new$Threshold)

#filter threshold/cut-off of no interest
#Example
df_new <-  df_new[df_new$Threshold !="0.0005",]
df_new <-  df_new[df_new$Threshold !="0.005",]

####################### GGPLOT the barplots for F1 accuracy for classifiers
classfier_colors <- c("kaiju"= "azure2" ,"kraken" = "#BDCDD6" , "kraken.gtdb" = "#93BFCF" , "metaphlan"= "#6096B4" ,"metaphlan4"= "azure4")

png("Classification_F1score.png", res = 300, height = 1500, width = 3000)
plot <- ggplot(df_new, aes(x= Threshold, y= F1_score, fill= Classifier)) +
  geom_boxplot() +
  ggtitle("F1 score of Classifiers at Species level")+
  xlab("Relative abundance cut-off")+
  ylab("precision")+
  scale_fill_manual(values= classfier_colors)+
  theme_bw()+facet_wrap(~Classifier, strip.position="bottom", labeller = label_both) 
plot +
  theme(panel.spacing = unit(2, "lines"),panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12, angle= 90, vjust=0.2),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), 
        legend.text=element_text(size = 12), legend.title = element_text(size = 12),
        strip.text.x = element_text(size = 12))  
dev.off()
