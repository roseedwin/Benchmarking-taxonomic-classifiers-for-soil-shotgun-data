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

df <- combined_df[,-c(1)]

df <- df %>%
  separate(Column, c("Sample","Classifier"), sep= "_")
  
df1 <-  df[,c(1,2,5,6,8)]
df1 <-df1[df1$Classifier != "metaphlan4",] 
df1 <-df1[df1$Classifier != "metaphlan",] 
df1 <-df1[df1$Threshold != "0.1",]
df1 <-df1[df1$Threshold != "5e-04",]
df1 <-df1[df1$Threshold != "0.005",]

df1R <- subset(df1, Run == 1)
df1R <- df1R[,-c(4)]

df2R <- subset(df1, Run == 2)
df2R <- df2R[,-c(4)]


df1R <-  df1R %>%
  pivot_wider(names_from = Classifier, values_from = F1_score)

  df_long <- df1R %>%
  pivot_longer(cols = c("kaiju", "kraken", "kraken.gtdb"),
               names_to = "Tool", values_to = "F1_score")
  
  df_long_ranked <- df_long %>%
  mutate(rank = dense_rank(F1_score))
classfier_colors <- c("#3A98B9","#87CBB9","#A0D995")
colr <- c("grey28","grey28", "grey28")

ggviolin(df_long_ranked, x = "Threshold", y = "rank", add = "jitter", color = "Tool", fill ="Tool")+
  scale_color_manual(values= colr) +
  scale_fill_manual(values= classfier_colors) +
  stat_compare_means(method ="kruskal.test", label="p.format",
                     label.y= max(df_long_ranked$rank) +1,
                     size =4, vjust =-1) 
 

kruskal.test(F1_score ~ Tool, data = df_long_ranked)
kruskal.test(F1_score ~ Threshold* Tool, data = df_long_ranked)

vec2 <- df_long_ranked %>% 
  group_by(Threshold)  %>% 
  rstatix::dunn_test(F1_score ~ Tool, p.adjust.method = "bonferroni", detailed =T)

vec3 <- df_long_ranked %>% 
  group_by(Tool)  %>% 
  rstatix::dunn_test(F1_score ~ Threshold, p.adjust.method = "bonferroni", detailed =T)

vec <- df_long_ranked %>% 
  unite("factor1", c("Tool", "Threshold"), remove = F) %>%
  rstatix::dunn_test(F1_score ~ factor1, p.adjust.method = "bonferroni", detailed =T) %>% 
  filter(p.adj.signif != "ns") %>% 
  arrange(p.adj)
