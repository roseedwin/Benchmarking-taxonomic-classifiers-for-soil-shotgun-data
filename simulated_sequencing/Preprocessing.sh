#####################################################################################
############################### PRE-PROCESSING ######################################
#####################################################################################

######################## decompressing the files from fastq.gz to fastq  ############
# For all files with .gz; do command "gunzip file.fastq.gz"; done!

for f in *.gz; do gunzip $f; done 

######################## creating files with names #############################
#all i know is that it selects the first 14 letters of the file name and make names.txt
#then create a unduplicated file 

for filename in *.fastq
  do
yolo=$(echo "$filename" | rev | cut -c 14- | rev)
echo $yolo >> names.txt
done
sort -u  names.txt >> names_undup.txt
rm names.txt


for item in $(less names_undup.txt)
do
g="${item%%_L*}"
  echo "${g}">>names.txt
done
uniq names.txt names_uniq.txt
sort names_uniq.txt
rm names.txt

cp names_undup.txt /data/Food/analysis/R1111_Soilmetadata/Main_data/Sample
cp names_uniq.txt /data/Food/analysis/R1111_Soilmetadata/Main_data/Sample

####################### combining lanes ########################################

for filename in $(less names_uniq.txt)
do
cat "$filename"_L001_R1_001.fastq "$filename"_L002_R1_001.fastq "$filename"_L003_R1_001.fastq "$filename"_L004_R1_001.fastq > "$filename"_merged_R1.fastq
cat "$filename"_L001_R2_001.fastq "$filename"_L002_R2_001.fastq "$filename"_L003_R2_001.fastq "$filename"_L004_R2_001.fastq > "$filename"_merged_R2.fastq
done

mkdir individual_lanes

mv *001.fastq individual_lanes/

############################## RENAME FILES #######################################
for i in *_R1.fastq; do
  mv "$i" "$(echo "$i" | sed 's/_S[0-9]\+_merged_/_/')"
done
for i in *_R2.fastq; do
  mv "$i" "$(echo "$i" | sed 's/_S[0-9]\+_merged_/_/')"
done

############################## CREATE SAMPLES.TXT #######################################
#!/bin/bash
# create an array of unique prefixes from filenames
prefixes=($(ls *_1.fastq | sed 's/_1\.fastq//' | sort -u))

# write prefixes to samples.txt
printf '%s\n' "${prefixes[@]}" > samples.txt

############################### CLUMPIFY ######################################

module load bbmap/38.22
mkdir Clumpify/

for i in $(cat samples.txt); do
clumpify.sh in="$i"_1.fastq in2="$i"_2.fastq out=Clumpify/"$i"_1.fastq out2=Clumpify/"$i"_2.fastq
done
module unload bbmap/38.22

############################### TRIMGALORE ######################################

module load trimgalore/0.6.1
mkdir TrimmedFastQ/

for i in $(cat samples.txt); do
trim_galore --paired Clumpify/"$i"_1.fastq Clumpify/"$i"_2.fastq --fastqc -j 8 -q 30 -o "$i"_trimout
module unload trimgalore/0.6.1

    mv "$i"_trimout/"$i"_1_val_1.fq TrimmedFastQ/"$i"_1.fastq
    mv "$i"_trimout/"$i"_2_val_2.fq TrimmedFastQ/"$i"_2.fastq
    rm -r "$i"_trimout

done
############################### BBMERGE ######################################

module load bbmap/38.22
mkdir BBMerge/

for i in $(cat samples.txt); do
bbmerge.sh in=TrimmedFastQ/"$i"_1.fastq in2=TrimmedFastQ/"$i"_2.fastq out=Clumpify/"$i"_1.fastq out2=Clumpify/"$i"_2.fastq
done
module unload bbmap/38.22
