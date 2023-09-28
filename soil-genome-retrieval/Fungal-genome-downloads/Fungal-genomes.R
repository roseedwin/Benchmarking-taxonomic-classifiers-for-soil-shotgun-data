if(!require(rentrez)){
  install.packages("rentrez")
  library(rentrez)
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

if(!require(knitr)){
  install.packages("knitr")
  library(knitr)
}

##FUNGAL GENOMES

#Genomic assembly submissions are divided in four categories, 
#      small complete genomes (eg. mitochondria), 
#      large complete genomes (eg. chromosomes), 
#      incomplete genomes (eg. WGSs), 
#      high throughput genome sequences (eg. BACs).
#The names of fungal species from soil was obtained from MIAE database maintained in INRA.
#using REntrez to extract the genomes belonging the list of names in 'Fungal_names' from NCBI database

#Creating a file with the filters to identify the accession numbers
Fungal_names <- read.csv("Fungal_species_MIAE.csv", header =FALSE)
colnames(Fungal_names)[1] <- "Fungi"

taxon_id = Fungal_names$Fungi
###########USING FILTERS ############
sink(file="output_rep_chr.txt")
for (X in taxon_id) {
  tax <-  entrez_search(db ="assembly", term=paste("Complete_genome[ASLV] OR Chromosome[ASLV] OR Representative_genome[RCAT] OR Reference_genome[RCAT] AND ", X, "[ORGN] AND Fungi [FILT]",sep="") )
  
  if(is.character(tax$ids)&& (tax$ids =1)) 
  {out  <-paste0("" ,X, sep="")
  print(out) }
} 
sink(file = NULL)

Rep_chr_genomes <- read.table("output_rep_chr.txt")
#119 fungal hits

############ TO OBTAIN ACCESSION NUMBERS ########################
Rep_chr_genomes <- read.table("output_rep_chr.txt", header=F)
Fungi_genomes <- Rep_chr_genomes[2]
colnames(Fungi_genomes)[1] <- "Fungi"

Fungi_genomes <- Fungi_genomes %>% distinct()

taxon_id = Fungi_genomes$Fungi

sink(file="genome_accession.txt")
for (X in taxon_id) {
  tax <-  entrez_search(db ="assembly", term=paste("Complete_genome[ASLV] OR Chromosome[ASLV] OR Representative_genome[RCAT] OR Reference_genome[RCAT] AND ", X, "[ORGN] AND Fungi [FILT]",sep="") )
  taxize_sum <- entrez_summary(db ="assembly", id =tax$ids )
  
  if (length(taxize_sum)==2 ) 
  {out  <-paste0( X, ":",lapply(taxize_sum, `[[`, 6), sep=";")
  print(out) } 
  
  if (length(taxize_sum)==3) 
  {out  <-paste0( X, ":",lapply(taxize_sum, `[[`, 4), sep=";")
  print(out) } 
  
  if (is.character(taxize_sum$assemblyaccession)==TRUE)
  {out  <- paste0( X, ":",taxize_sum$assemblyaccession, sep=";")
  print(out)}
}
sink(file=NULL)

#Data cleaning in excel
#In the absence of super-awesome R skills, I did data cleaning in excel before bringing it back in a csv format
#Sorry I tried, just couldnt waste more time
genomes <- read.csv("genome_accession.csv", header =F)

colnames(genomes)[1] <- "Fungi"
colnames(genomes)[2] <- "Accession_no"

#double checking data in R along with fetching metadata
accession_no <- genomes$Accession_no

sink(file="metadata_fungi.csv")
for (X in accession_no) {
  tax <-  entrez_search(db ="assembly", term=paste(X, "[ASAC]",sep="")) 
    
  query <- entrez_summary(db ="assembly" , id=paste(tax$ids))  
    
    if(is.character(query$speciesname))  
    {id <-extract_from_esummary(query, c("assemblyaccession","refseq_category", "organism", "speciesname","taxid", "assemblystatus"))
    print(knitr::kable(t(id), row.names=NA,escape=FALSE))
    }
}
sink(file=NULL)

metadata_fungi <- read.csv("metadata_fungi.csv", header=T, sep = "|")
metadata_fungi <- metadata_fungi[- grep("assemblyaccession", metadata_fungi$assemblyaccession),]  
metadata_fungi <- metadata_fungi[- grep(":-----", metadata_fungi$taxid),]  
metadata_fungi <- metadata_fungi[,2:7]

#de-duplicate the results
metadata_fungi_dd <- metadata_fungi[!duplicated(metadata_fungi),]

write.csv(metadata_fungi_dd,"metadata_fungi_final.csv")

Accession_numbers <- metadata_fungi_dd$assemblyaccession
write.csv(Accession_numbers,"Accession_numbers.csv")

##DOWNLOADING FUNGAL GENOMES USING COMMAND LINE TOOLS
#Once obtained accession numbers now we can download using command line tools
for i in $(cat Accession_numbers.csv);do
datasets download genome accession "$i" --filename "$i".zip
done
