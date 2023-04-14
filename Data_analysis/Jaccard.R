library(ggplot2)
library(vegan)
library(wesanderson)
library(tibble)
library(devtools)
library(ggfortify)
library(ggalt)
library(ggpubr)
library(PCAtools)
library(reshape2)

observed_results <- read.csv("observed_results.csv", header=T)
observed_results <- observed_results[,c(2,3,4)]
observed_results$species <- gsub('[^[:alnum:] ]','',observed_results$species)

expected_results <- read.csv("expected_results.csv", header=T)
expected_results <- expected_results[,c(2,3,4)]
expected_results$species <- gsub('[^[:alnum:] ]','',expected_results$species)

data <- rbind(observed_results,expected_results)
rm(expected_results, observed_results)

c<-dcast(data,species~sample_name, value.var= 'abundance', fun.aggregate = sum)
c[is.na(c)] <- 0
write.csv(c,"c.csv")


data <- read.csv("c.csv", header= T, row.names = 1)
data <- t(data)

#metadata file
map <- read.csv("map.csv", header = T, row.names = 1)

library(vegan)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)

mycols <- wes_palette("Moonrise3")





vec1 <- map[,1]
#vec2 <- map[,4]

meta.mds.sl=metaMDS(data,distance = "jaccard")
plot(meta.mds.sl$points)

write.csv(meta.mds.sl$points,"mds_spec.csv")

df<-read.csv("mds_spec.csv",header=T)

### Plot the distance matrix with ggplot2 ###
Sp<-ggplot(df,aes(MDS1,MDS2))+  
  scale_fill_manual(values = mycols) +  
  scale_color_manual(values = mycols) +  
  scale_shape_manual(values = c(5, 7, 10, 13, 1, 2, 3, 4, 6, 8, 9, 11)) +  
  theme_bw() +  
  theme(    
    legend.title=element_text(face="bold",size=12.5),
    axis.title.x=element_text(face="bold",size=12.5),
    
    axis.title.y=element_text(face="bold",size=12.5),
    
    plot.title=element_text(face="bold",size=15),
    
    legend.text = element_text(size=12.5),
    
    axis.text.y = element_text(size=10),
    
    axis.text.x = element_text(size=12.5),
    
    axis.title = element_text(size=12.5)) +
  
  theme(panel.background = element_rect(colour="black")) +
  
  ggtitle("Species dissimilarity Kaiju") +
  
  geom_point(aes(color = vec1), size=5) +
  
  stat_ellipse(aes(color=vec1))



png("poster.png", res = 300, height = 3000, width = 3000)



Sp



dev.off()


