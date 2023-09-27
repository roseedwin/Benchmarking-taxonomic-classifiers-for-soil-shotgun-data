#!/bin/sh
#SBATCH --error err_tax
#SBATCH --output out_tax
#SBATCH --job-name MS_Taxonomy
#SBATCH --mail-user xyz@gmail.com
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH -p Priority,Background
#SBATCH -N 1
#SBATCH --cpus-per-task=10

for j in Run1 Run2 Run3 Run4 Run5 Run6 Run7 Run8 Run9 Run10; do

cd "$j"

mkdir TAXONOMY
################################################### KAIJU ###################################################
mkdir TAXONOMY/KAIJU
module load kaiju/1.7.4

# Record start time
start_time_Kaiju=$(date +%s)

for ((sample=1; sample<=20; sample++)); do
  i="sample${sample}"
kaiju -t /kaijudb/nodes.dmp \
      -f /kaiju_db_nr_euk.fmi \
      -i "$i"_R1.fastq.gz -j "$i"_R2.fastq.gz -z 10 -o TAXONOMY/KAIJU/"$i"_out

kaiju2table -t /kaijudb/nodes.dmp \
            -n /kaijudb/names.dmp \
            -r species -o TAXONOMY/KAIJU/"$i"_kaiju_species.tsv TAXONOMY/KAIJU/"$i"_out
            
kaiju2table -t /kaijudb/nodes.dmp \
            -n /kaijudb/names.dmp \
            -r genus -o TAXONOMY/KAIJU/"$i"_kaiju_genus.tsv TAXONOMY/KAIJU/"$i"_out
            
kaiju2table -t /kaijudb/nodes.dmp \
            -n /kaijudb/names.dmp \
            -r family -o TAXONOMY/KAIJU/"$i"_kaiju_family.tsv TAXONOMY/KAIJU/"$i"_out
            
kaiju2table -t /kaijudb/nodes.dmp \
            -n /kaijudb/names.dmp \
            -r order -o TAXONOMY/KAIJU/"$i"_kaiju_order.tsv TAXONOMY/KAIJU/"$i"_out
            
kaiju2table -t /kaijudb/nodes.dmp \
            -n /kaijudb/names.dmp \
            -r class -o TAXONOMY/KAIJU/"$i"_kaiju_class.tsv TAXONOMY/KAIJU/"$i"_out
            
kaiju2table -t /kaijudb/nodes.dmp \
            -n /kaijudb/names.dmp \
            -r phylum -o TAXONOMY/KAIJU/"$i"_kaiju_phylum.tsv TAXONOMY/KAIJU/"$i"_out
            
done
module unload kaiju/1.7.4

cd TAXONOMY/KAIJU
rename 'kaiju_' '' *tsv
  categories=("order" "species" "genus" "class" "phylum" "family")

for category in "${categories[@]}"; do
	for x in $(ls *"$category".tsv | cut -d "_" -f1); do
    	awk -v var="$x" '{print var,$0}' OFS="\t" "$x"_"$category".tsv >> "$x"_samplename_added.report
    	done
   	for x in $(ls *samplename*); do
   	tail -n +2 $x >> "$category"_kaiju.txt
   	done

	rm *.report
	rm *"$category".tsv

    input_file="${category}_kaiju.txt"
    output_file="${category}_Kaiju.txt"
    
    sed -e 's/TAXONOMY\/KAIJU\///g' -e 's/_out//g' "$input_file" > "$output_file"
done

rm *kaiju.txt
cd ../..

# Record end time
end_time_Kaiju=$(date +%s)

# Calculate execution time
  execution_time_Kaiju=$((end_time_Kaiju - start_time_Kaiju))
  echo "Execution time for $j Kaiju: $execution_time_Kaiju seconds"

################################################### KRAKEN ###################################################
mkdir TAXONOMY/KRAKEN

# Record start time
start_time_Kraken=$(date +%s)

  for ((sample=1; sample<=20; sample++)); do
    i="sample${sample}"
  module load kraken2/2.1.1
  kraken2 --db /path_to_kraken_database --threads 10 \
        --paired "$i"_R1.fastq.gz "$i"_R2.fastq.gz --output TAXONOMY/KRAKEN/"$i"_kraken\
        --report TAXONOMY/KRAKEN/"$i"_report

  module load braken/2.2
  bracken -d /data/databases/kraken2_pluspf_20210517 \
        -i TAXONOMY/KRAKEN/"$i"_report \
        -l 'S' -o TAXONOMY/KRAKEN/"$i"_bracken.species.report

  bracken -d /data/databases/kraken2_pluspf_20210517 \
        -i TAXONOMY/KRAKEN/"$i"_report \
        -l 'G' -o TAXONOMY/KRAKEN/"$i"_bracken.genus.report

  bracken -d /data/databases/kraken2_pluspf_20210517 \
        -i TAXONOMY/KRAKEN/"$i"_report \
        -l 'F' -o TAXONOMY/KRAKEN/"$i"_bracken.family.report
        
  bracken -d /data/databases/kraken2_pluspf_20210517 \
        -i TAXONOMY/KRAKEN/"$i"_report \
        -l 'O' -o TAXONOMY/KRAKEN/"$i"_bracken.order.report

  bracken -d /data/databases/kraken2_pluspf_20210517 \
        -i TAXONOMY/KRAKEN/"$i"_report \
        -l 'C' -o TAXONOMY/KRAKEN/"$i"_bracken.class.report
        
  bracken -d /data/databases/kraken2_pluspf_20210517 \
        -i TAXONOMY/KRAKEN/"$i"_report \
        -l 'P' -o TAXONOMY/KRAKEN/"$i"_bracken.phylum.report

  done

  module unload braken/2.2
  module unload kraken2/2.1.1
  
################################################### tidying up reports

  cd TAXONOMY/KRAKEN
  rename 'bracken.' '' *report

  categories=("order" "species" "genus" "class" "phylum" "family")

  for category in "${categories[@]}"; do
    for x in $(ls *"$category".report | cut -d "_" -f1); do
        awk -v var="$x" '{print var,$0}' OFS="\t" "$x"_"$category".report >> "$x"_samplename_added.tsv
    done
    for x in $(ls *samplename*); do
        tail -n +2 $x >> "$category"_Kraken.txt
    done
    rm *.tsv
    rm *"$category".report
  done

  cd ../..

# Record end time
end_time_Kraken=$(date +%s)

# Calculate execution time
  execution_time_Kraken=$((end_time_Kraken - start_time_Kraken))
  echo "Execution time for $j Kraken: $execution_time_Kraken seconds"

################################################### CUSTOM KRAKEN ###################################################
mkdir TAXONOMY/KRAKEN_gtdb2

# Record start time
start_time_customKraken=$(date +%s)

  for ((sample=1; sample<=20; sample++)); do
    i="sample${sample}"
  module load kraken2/2.1.1
  kraken2 --db /path_to_database/kraken_gtdbtk_db2 --threads 10 \
        --paired "$i"_R1.fastq.gz "$i"_R2.fastq.gz --output TAXONOMY/KRAKEN_gtdb2/"$i"_kraken\
        --report TAXONOMY/KRAKEN_gtdb2/"$i"_report

  module load braken/2.2
  bracken -d /Kraken/kraken_gtdbtk_db2 \
        -i TAXONOMY/KRAKEN_gtdb2/"$i"_report \
        -l 'S' -o TAXONOMY/KRAKEN_gtdb2/"$i"_bracken.species.report

  bracken -d /Kraken/kraken_gtdbtk_db2 \
        -i TAXONOMY/KRAKEN_gtdb2/"$i"_report \
        -l 'G' -o TAXONOMY/KRAKEN_gtdb2/"$i"_bracken.genus.report
        
  bracken -d /Kraken/kraken_gtdbtk_db2 \
        -i TAXONOMY/KRAKEN_gtdb2/"$i"_report \
        -l 'C' -o TAXONOMY/KRAKEN_gtdb2/"$i"_bracken.class.report

  bracken -d /Kraken/kraken_gtdbtk_db2 \
        -i TAXONOMY/KRAKEN_gtdb2/"$i"_report \
        -l 'F' -o TAXONOMY/KRAKEN_gtdb2/"$i"_bracken.family.report

  bracken -d /Kraken/kraken_gtdbtk_db2 \
        -i TAXONOMY/KRAKEN_gtdb2/"$i"_report \
        -l 'P' -o TAXONOMY/KRAKEN_gtdb2/"$i"_bracken.phylum.report

  bracken -d /Kraken/kraken_gtdbtk_db2 \
        -i TAXONOMY/KRAKEN_gtdb2/"$i"_report \
        -l 'O' -o TAXONOMY/KRAKEN_gtdb2/"$i"_bracken.order.report
  done

  module unload braken/2.2
  module unload kraken2/2.1.1

################################################### tidying up reports
  cd TAXONOMY/KRAKEN_gtdb2
  rename 'bracken.' '' *report

  categories=("order" "species" "genus" "class" "phylum" "family")

  for category in "${categories[@]}"; do
    for x in $(ls *"$category".report | cut -d "_" -f1); do
        awk -v var="$x" '{print var,$0}' OFS="\t" "$x"_"$category".report >> "$x"_samplename_added.tsv
    done
    for x in $(ls *samplename*); do
        tail -n +2 $x >> "$category"_Kraken_gtdb.txt
    done
    rm *.tsv
    rm *"$category".report
  done

  cd ../..

# Record end time
end_time_customKraken=$(date +%s)

# Calculate execution time
  execution_time_customKraken=$((end_time_customKraken - start_time_customKraken))
  echo "Execution time for $j custom Kraken: $execution_time_customKraken seconds"

################################################### INTERLEAVED FASTQ ###################################################

  mkdir TAXONOMY/METAPHLAN3
  for ((sample=1; sample<=20; sample++)); do
  i="sample${sample}"
  module load bbmap/38.22
  reformat.sh in="$i"_R1.fastq.gz in="$i"_R2.fastq.gz out=TAXONOMY/METAPHLAN3/"$i"_interleaved.fastq.gz
  gunzip TAXONOMY/METAPHLAN3/*.gz
  module unload bbmap/38.22
  done

################################################### METAPHLAN3 ###################################################
# Record start time
start_time_Metaphlan3=$(date +%s)

  for ((sample=1; sample<=20; sample++)); do
    i="sample${sample}"
  module load humann3/3.0
  metaphlan TAXONOMY/METAPHLAN3/"$i"_interleaved.fastq --input_type fastq \
  --bowtie2db /data/databases/MetaPhlAn3 -x mpa_v30_CHOCOPhlAn_201901 --bowtie2out TAXONOMY/METAPHLAN3/"$i"_metaphlanbowtie2out.txt \
  --nreads 10 --unknown_estimation -o TAXONOMY/METAPHLAN3/"$i"_profile.txt 

  merge_metaphlan_tables.py TAXONOMY/METAPHLAN3/*_profile.txt -o TAXONOMY/METAPHLAN3/Metaphlan3.tsv
  module unload humann3/3.0
  done

# Record end time
end_time_Metaphlan3=$(date +%s)

# Calculate execution time
  execution_time_Metaphlan3=$((end_time_Metaphlan3 - start_time_Metaphlan3))
  echo "Execution time for $j Metaphlan3: $execution_time_Metaphlan3 seconds"

################################################### METAPHLAN4 ###################################################
mkdir TAXONOMY/METAPHLAN4

# Record start time
start_time_Metaphlan4=$(date +%s)

  for ((sample=1; sample<=20; sample++)); do
    i="sample${sample}"
  module load humann3/3.0

  /install/software/restart/py3/humann_3.5/bin/metaphlan TAXONOMY/METAPHLAN3/"$i"_interleaved.fastq \
  --bowtie2db /data/databases/MetaPhlan4 -x mpa_vJan21_CHOCOPhlAnSGB_202103 --input_type fastq \
  --unclassified_estimation --nproc 10 --output_file TAXONOMY/METAPHLAN4/"$i"_profile.txt

  merge_metaphlan_tables.py TAXONOMY/METAPHLAN4/*_profile.txt -o TAXONOMY/METAPHLAN4/Metaphlan4.tsv

  module unload humann3/3.0
  done

# Record end time
end_time_Metaphlan4=$(date +%s)

# Calculate execution time
  execution_time_Metaphlan4=$((end_time_Metaphlan4 - start_time_Metaphlan4))
  echo "Execution time for $j Metaphlan4: $execution_time_Metaphlan4 seconds"


###############################################################################################################

echo "EXECUTION TIME TAKEN FOR DIFFERENT CLASSIFIERS RUN"
echo "Execution time for $j Metaphlan4: $execution_time_Metaphlan4 seconds"
echo "Execution time for $j Metaphlan3: $execution_time_Metaphlan3 seconds"
echo "Execution time for $j custom Kraken: $execution_time_customKraken seconds"
echo "Execution time for $j Kraken: $execution_time_Kraken seconds"
echo "Execution time for $j Kaiju: $execution_time_Kaiju seconds"

cd ..
done
