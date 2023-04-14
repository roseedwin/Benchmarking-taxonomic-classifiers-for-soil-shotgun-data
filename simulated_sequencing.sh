dir=$(pwd)
mkdir simulated_sequencing 
cd simulated_sequencing 

#############################################################################################################################################################################################   
# set up sample files for simulated data including, number of genotypes, reads and names of samples
module load R/4.0.2
Rscript /path_to../randomdiversity.R
module unload R/4.0.2

tail -n +2 sample.txt > sample2.txt
tail -n +2 reads.txt > reads2.txt ######## creates random reads in a truncated log normal distribution
tail -n +2 genotypes.txt > genotypes2.txt ###### names for each output genome

rm -f sample.txt
rm -f reads.txt
rm -f genotypes.txt
#############################################################################################################################################################################################   
#THESE STEPS TAKE ONLY FEW SECS HENCE RUN TOGETHER

# pick random files from available fasta files for mock sequencing data
# create a folder for each sample  (Sample 1-20 in the simulated sequencing folder)
cat sample2.txt | xargs -L 1 mkdir

# 1. create commands for to pick n number of fasta files for each folder
awk '{print "find /data/Food/analysis/R1111_Soilmetadata/Main_data/insilicoseq/FASTA_files/*.fa -type f | shuf -n "$0" > genotypes_list.txt"}' genotypes2.txt >> genotypes_list.txt 

# 2. Pick one line from file and place in each subdirectory
for d in */; do
(cd $d && cat ../genotypes_list.txt | shuf -n 1  > genotypes_list.sh)
done

#THIS STEP TOOK 2 HRs and 19 MINUTES TO RUN
# 3. For each subdirectory, carry out commands from part 1
for d in */;  
do 
(cd $d
sh genotypes_list.sh

awk -F '\t' '{print"head -n 1 "$1"  | cut -f 2 >> species.txt"}' genotypes_list.txt > species.sh
#extract names of fasta sequences
sh species.sh
 
# move all fasta files in genotypes_list.txt to each sample sudirectory and cat the fasta files into one big fasta file
xargs -a genotypes_list.txt cp -t .
cat *.fa > sample.fasta 

# extract only sequnece name of fasta file - accession number
awk '{print $1}' species.txt > Genomenames.txt

# combine names into one row. 
cat Genomenames.txt | sed 's/>//' > genotypes.txt
cat genotypes.txt | paste -s -d "_" > genotypes2.txt)
done

#THIS STEP TOOK 30 MINUTES TO RUN
# 4. rename sample.fasta after parent directory and move up directory
find "../simulated_sequencing/" -type f -iname 'sample.fasta' -exec sh -c '
    path="${1%/*}"; filename="${1##*/}";
    cp -nv "${1}" ../"simulated_sequencing/${path##*/}.fasta" ' sh_cp {} \;

# 5. rename genotypes.txt file after parent directory and move up directory
find "../simulated_sequencing/" -type f -iname '*genotypes2.txt' -exec sh -c '
   path="${1%/*}"; filename="${1##*/}";
  cp -nv "${1}" ../"simulated_sequencing/${path##*/}${filename}" ' sh_cp {} \;


# 6. Remove subdirectories and content. 
rm -rf ../simulated_sequencing/*/

# merge all genome files (contain accession numbers of genotypes in each fasta file) and include file name 
for y in $(ls sample*genotypes2.txt);do a="$(cut -d 'g' -f 1 <<< "$y")";for x in $(cat "$y" );do echo "$a" "$x"  >> genome_list.txt ;done;done
rm sample*genotypes2.txt
   
