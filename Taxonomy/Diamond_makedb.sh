#!/bin/sh
#SBATCH --error err_makedb.aa
#SBATCH --output out_makedb.aa
#SBATCH --job-name makedb.aa.fa
#SBATCH --cpus-per-task=3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xyz@gmail.com
#SBATCH -p Priority,Background
#SBATCH --ntasks=1
#SBATCH -N 1

#Download the aa fasta file
wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz
tar -xzvf gtdb_proteins_aa_reps.tar.gz

for i in protein_faa_reps/archaea/*.faa.gz;
do
gunzip "$i"
done

#Make diamond database
cat protein_faa_reps/archaea/*.faa > protein_faa_reps/archaea/archaea.faa
mv protein_faa_reps/archaea/archaea.faa protein_faa_reps/
gzip protein_faa_reps/archaea/*.faa

for i in protein_faa_reps/bacteria/*.faa.gz;
do 
gunzip "$i"
done 

#Add to the diamond database
cat protein_faa_reps/bacteria/*.faa > protein_faa_reps/bacteria/bacteria.faa
mv protein_faa_reps/bacteria/bacteria.faa protein_faa_reps/

gzip protein_faa_reps/bacteria/*.faa

cd protein_faa_reps
# Concatenate all .faa files into protein.faa
cat *.faa > protein.faa

# Create Diamond database using protein.faa
module load diamond/2.0.12
diamond --makedb --in protein.faa -d gtdb_diamond
module unload diamond/2.0.12


