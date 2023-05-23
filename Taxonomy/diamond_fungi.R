#FOR SwissProt UniProtKB PROTEIN FASTA
library(rentrez)

#DOWNLOADING ALL SWISSPROT PROTEIN SEQUENCES FOR MUCOROMYCOTA, ASCOMYCOTA AND BASIDIOMYCOTA
#MUCOROMYCOTA
tax <-  entrez_search(db ="protein", term="Mucoromycota[Organism] OR mucoromycota[All Fields]) AND swissprot[filter]", retmax =144)
ids <-  tax$ids
for (id in ids) {
  fasta_file <- entrez_fetch(db="protein", id= id, rettype="fasta")
  taxize_sum <- entrez_summary(db ="protein", id = id )
  x <- taxize_sum$accessionversion
  write(fasta_file, file=paste("",x,".fasta", sep=""))
}

#ASCOMYCOTA
tax <-  entrez_search(db ="protein", term="Ascomycota[Organism] OR ascomycota[All Fields]) AND swissprot[filter]", retmax =33545)
ids <-  tax$ids
for (id in ids) {
  fasta_file <- entrez_fetch(db="protein", id= id, rettype="fasta")
  taxize_sum <- entrez_summary(db ="protein", id = id )
  x <- taxize_sum$accessionversion
  write(fasta_file, file=paste("",x,".fasta", sep=""))
}

#Repeat till all files are downloaded (all hail ChatGPT)
repeat {
  tryCatch({
    # Your code lines here
    
    accession_version_to_find <- get("x")
    search_results <- entrez_search(db = "protein", term = accession_version_to_find, retmax = 1)
    if (search_results$count > 0) {
      id_found <- search_results$ids[1]
      print(paste("ID found for accession version", accession_version_to_find, "is", id_found))
    } else {
      print(paste("No ID found for accession version", accession_version_to_find))
    }
    
    # Read the ids.csv file
    ids_data <- read.csv("ids.csv", header = TRUE)
    specific_id <- get("id_found")
    start_row <- which(ids_data$x == specific_id)
    ids_subset <- ids_data$x[start_row:length(ids_data$x)]
    
    # Iterate over the remaining IDs
    for (id in ids_subset) {
      fasta_file <- entrez_fetch(db = "protein", id = id, rettype = "fasta")
      taxize_sum <- entrez_summary(db = "protein", id = id)
      x <- taxize_sum$accessionversion
      write(fasta_file, file = paste("", x, ".fasta", sep = ""))
    }
    
    # If no error occurs, break out of the loop
    break
  }, error = function(e) {
    # Error handling code here
    
    # Print the error message
    cat("Error:", conditionMessage(e), "\n")
    
    # Wait for a few seconds before retrying
    Sys.sleep(3)
  })
}

#BASIDIOMYCOTA
tax <-  entrez_search(db ="protein", term="Basidiomycota[Organism] OR basidiomycota[All Fields]) AND swissprot[filter]", retmax =33545)
ids <-  tax$ids
write.csv (ids, "ids.csv")
for (id in ids) {
  fasta_file <- entrez_fetch(db="protein", id= id, rettype="fasta")
  taxize_sum <- entrez_summary(db ="protein", id = id )
  x <- taxize_sum$accessionversion
  write(fasta_file, file=paste("",x,".fasta", sep=""))
}























