#Pre-analysis for randomdiversity.R

df <- Merged_Samples_final[,c(1,2,4)]
df <- df[df$Sample !="blank",]
df <- df[df$Sample !="AH2",]
df <- df[df$Sample !="BL5",]
df <- df[df$Sample !="PA17",]

NX1 <- subset(df, Run == "NX1b")
NX1$raw_reads <- NX1$raw_reads * 2
NX2 <- subset(df, Run == "NX2")
NX2$raw_reads <- NX2$raw_reads * 2
NX3 <- subset(df, Run == "NX3")
NX3$raw_reads <- NX3$raw_reads * 2
NX4 <- subset(df, Run == "NX4")
NX4$raw_reads <- NX4$raw_reads * 2
NX5 <- subset(df, Run == "NX5")
NX5$raw_reads <- NX5$raw_reads * 2
NX6 <- subset(df, Run == "NX6")
NX6$raw_reads <- NX6$raw_reads * 2
NX7 <- subset(df, Run == "NX7")
NX7$raw_reads <- NX7$raw_reads * 2

df_names <- c("NX1","NX2","NX3","NX4","NX5","NX6","NX7")

mean_logs <- c()
sd_logs <- c()
mean_sum <- c()

sum_df <- sum(NX1$raw_reads)
mean_sum <- mean(sum_df)

for (name in df_names) {
  df <- get(name)
  
  sum_df <- sum(df$raw_reads)
  mean_sum <- c(mean_sum, sum_df)
  
  mean_log= log(mean(df$raw_reads))
  mean_logs <- c(mean_logs, mean_log)
  
  sd_log <- log(sd(df$raw_reads))
  sd_logs <- c(sd_logs, sd_log)
  
  cat(paste0("mean log ", name, ": ", mean_log, "\n"))
  cat(paste0("sd log ", name, ": ", sd_log, "\n"))
}

mean_mean_log <- mean(mean_logs)
print(mean_mean_log)

mean_sd_log <- mean(sd_logs)
print(mean_sd_log)

mean_sum <- mean(mean_sum)
print(mean_sum)

#install.packages("fitdistrplus")
library("fitdistrplus")

png("log_norm_trunc_NX7.png", res = 300, height = 1000, width = 1000)
hist(NX7$raw_reads, breaks =10, freq= FALSE, main= "Histogram of NX7")
fit <- fitdist(NX7$raw_reads, "lnorm")
curve(dlnorm(x, meanlog = fit$estimate[1], sdlog = fit$estimate[2]), add = TRUE, col ="red")
dev.off()

shapiro.test(log(NX7$raw_reads))

png("Distribution of reads across the Runs.png", res = 300, height = 700, width = 1200)
plot <- ggplot(df, aes(x= Run, y= raw_reads)) +
  geom_boxplot() +
  ggtitle("Distribution of reads across the Runs")+
  xlab("Runs")+
  ylab("Read depth")+
  theme_bw()
plot
dev.off()

