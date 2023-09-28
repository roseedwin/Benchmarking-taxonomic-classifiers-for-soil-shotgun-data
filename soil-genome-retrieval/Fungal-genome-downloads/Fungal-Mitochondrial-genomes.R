if(!require(rentrez)){
  install.packages("rentrez")
  library(rentrez)
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
sink(file="output_refseq_complete.txt")
for (X in taxon_id) {
  
  tax <-  entrez_search(db ="nuccore", term=paste( X, "[ORGN] AND 
                                                Fungi [FILT] AND 
                                                complete genome[TITL] AND
                                                Refseq[KYWD]", sep=""))
  if(is.character(tax$ids)&& (tax$ids =1)) 
  {out  <-paste0("" ,X, sep="")
  print(out) }
}  
sink(file = NULL)

RefSeq_complete_genomes <- read.table("output_refseq_complete.txt")
#41 genomes out of 258 names


#################################### TO DOWNLOAD GENOMES #############################################

############ REFSEQ COMPLETE GENOMES ########################
dir.create("Fungi_RefSeq_complete_mitochondrion")

RefSeq_complete_genomes <- RefSeq_complete_genomes[2]
RefSeq_complete_genomes <- RefSeq_complete_genomes %>% distinct()
colnames(RefSeq_complete_genomes)[1] <- "Fungi"

taxon_id = RefSeq_complete_genomes$Fungi

for (X in taxon_id) {
  
  tax <-  entrez_search(db ="nuccore", term=paste( X, "[ORGN] AND 
                                                Fungi [FILT] AND 
                                                complete genome[TITL] AND
                                                Refseq[KYWD]", sep=""))
  
  tax$ids
  
  for (i in tax$ids)
  {
    fasta_file <- entrez_fetch(db="nuccore", id=paste(i), rettype="fasta")
    write(fasta_file, file=paste("Fungi_RefSeq_complete_mitochondrion/",i,"_",X,".fasta", sep=""))
  }
}

#All mitochondrial sequences
#After removing duplicates 39 fungal mitochondrial genomes 


########################### TO OBTAIN METADATA ######################

genomes <- read.csv("mitochondrion_names.csv", header =F)

colnames(genomes)[2] <- "Fungi"
colnames(genomes)[1] <- "Accession"

#double checking data in R along with fetching metadata
accession <- genomes$Accession

sink(file="metadata_fungi_mitochondrion.csv")
for (X in accession) {
  tax <-  entrez_search(db ="nuccore", term=paste(X, "[ACCN]",sep="")) 
  
  query <- entrez_summary(db ="nuccore" , id=paste(tax$ids))  
  
  if(is.character(query$organism)==T)  
  {id <-extract_from_esummary(query, c("title","sourcedb", "organism", "strain","taxid", "genome", "accessionversion"))
  print(knitr::kable(t(id), row.names=NA,escape=FALSE))
  }
}
sink(file=NULL)

metadata <- read.csv("metadata_fungi_mitochondrion.csv", header=T, sep = "|")
metadata <- metadata[- grep("title", metadata$title),]  
metadata <- metadata[- grep(":-----", metadata$taxid),]  
metadata <- metadata[,2:8]

#de-duplicate the results
metadata_dd <- metadata[!duplicated(metadata),]

write.csv(metadata_dd,"metadata_fungi_mitochondrion_final.csv")
