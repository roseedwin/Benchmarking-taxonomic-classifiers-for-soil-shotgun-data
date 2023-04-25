library(dplyr)
library(magrittr)

sample <- as.data.frame(matrix(paste0("sample", "", 01:20), ncol = 1, nrow = 20)) %>% 
  dplyr::rename(sample =V1) %>%
  dplyr::mutate(n_genotypes =sample(1000:2500,20, replace=T ))

list <- c(sample[,2])

write.csv(list, "genotypes.txt", row.names=F)
write.csv(sample[,1],"sample.txt", quote=FALSE,row.names=F)

#This creates two files "sample.txt" with 20 rows of sample names, and "genomes.txt" with 20 rows of numbers between 1000-2500

#################################################################
#sequencing data is truncated log normal distribution
################################################################
#script for calculating the meanlog, sdlog and checking the distribution of data is pre_ramdomdiversity.R

# distribution based on input reads in statsdada2 files from various sequencing runs
# clean lognormal distribution reads accross samples and features
# use truncated log normal distribution to prevent most samples centering at 0 reads
library(EnvStats)

rand_vect_cont <- function(N, S, sdlog, min,max) {
  vec <-  data.frame(reads=abs(EnvStats::rlnormTrunc(N, meanlog=17.5, sdlog, min, max)))
  all <- (vec / sum(vec)) * S
  even <- all %>% 
    dplyr::mutate(reads=as.integer(reads)) %>%
    dplyr::mutate(reads=ifelse((reads%%2)==0,reads, reads+1)) 
  even
}

#Total reads will be the sum of the two fastq files, so if you want the max read to be 35 million for fastq per sample, 
#you want the total reads to 75 million reads to account for the forward and reverse read/fastq file

total_reads <-  75607898
vec <- rand_vect_cont(20, S=total_reads,sdlog=16.9, min=10000000, max=100000000)
write.csv(vec, "reads.txt", row.names=F)

################## normal distribution of quality scores for HiSeq runs
library(EnvStats)
quality_scores <- data.frame(high_q = as.integer(rnormTrunc(20, mean = 35, sd = 2, min = 28, max = 36)), 
                             low_q =as.integer(rnormTrunc(20, mean = 28, sd = 2, min = 20, max = 28)))
write.table(quality_scores, "quality_scores.txt",quote=FALSE, col.names=F, row.names=F,sep = "\t")

