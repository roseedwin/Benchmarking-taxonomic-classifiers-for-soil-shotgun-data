dir=$(pwd)
cd insilicoseq/TrimmedFastQ

mkdir TAXONOMY

#RUN REGULAR KRAKEN
mkdir TAXONOMY/KRAKEN

echo " ************************ 1. Running Kraken 2.1.1 and Bracken 2.2 ***************************************"
echo "*********************************************************************************************************"

module load kraken2/2.1.1
ls *_R1.fastq | awk -F '_R1.fastq' '{print "kraken2 --db /data/databases/kraken2_pluspf_20210517 --threads 20 --paired "$0" "$1"_R2.fastq --output TAXONOMY/KRAKEN/"$1"_kraken --report TAXONOMY/KRAKEN/"$1"_report"}' > kraken2.sh
sh kraken2.sh

module load braken/2.2
ls *_R1.fastq | awk -F '_R1.fastq' '{print "bracken -d /data/databases/kraken2_pluspf_20210517 -i TAXONOMY/KRAKEN/"$1"_report -o TAXONOMY/KRAKEN/"$1"_bracken.species.report"}' > bracken.sh
sh bracken.sh

module unload braken/2.2
module unload kraken2/2.1.1

find TAXONOMY/KRAKEN -type f ! -name '*species.report' -delete

echo "**************************** Completed Running Kraken and Bracken ***************************************"
echo "*********************************************************************************************************"
echo " *************************** 2. Running Kaiju 1.7.4 *****************************************************"

#RUN KAIJU
mkdir TAXONOMY/KAIJU
module load kaiju/1.7.4

ls *_R1.fastq | awk -F '_R1.fastq' '{print "kaiju -t /data/Food/analysis/R1111_Soilmetadata/Main_data/databases/kaijudb/nodes.dmp -f /data/Food/analysis/R1111_Soilmetadata/Main_data/databases/kaijudb/kaiju_db_nr_euk.fmi -i "$0" -j "$1"_R2.fastq -z 25 -o TAXONOMY/KAIJU/"$0"_out"}' > runkaiju.sh
sh runkaiju.sh
 
ls TAXONOMY/KAIJU/*_out | awk -F '_out' '{print "kaiju-addTaxonNames -t /data/Food/analysis/R1111_Soilmetadata/Main_data/databases/kaijudb/nodes.dmp -n /data/Food/analysis/R1111_Soilmetadata/Main_data/databases/kaijudb/names.dmp -i "$0" -r phylum,class,order,family,genus,species -u -o "$0"_kaiju.txt"}' > annot.sh
sh annot.sh

ls TAXONOMY/KAIJU/*_out | awk -F '_out' '{print "kaiju2table -t /data/Food/analysis/R1111_Soilmetadata/Main_data/databases/kaijudb/nodes.dmp -n /data/Food/analysis/R1111_Soilmetadata/Main_data/databases/kaijudb/names.dmp -r species -o "$0"_kaiju_summary.tsv "$0" "}' > tables.sh
sh tables.sh

find TAXONOMY/KAIJU -type f ! -name '*.tsv' -delete

module unload kaiju/1.7.4

echo "**************************** Completed Running Kaiju ****************************************************"
echo "*********************************************************************************************************"
echo " **************************** 3. Running Custom_Kraken **************************************************"

#RUN CUSTOM KRAKEN
mkdir TAXONOMY/KRAKEN_gtdb

module load kraken2/2.1.1
ls *_R1.fastq | awk -F '_R1.fastq' '{print "kraken2 --db /data/Food/analysis/R1111_Soilmetadata/Main_data/Kraken/kraken_gtdbtk_db/ --threads 25 --paired "$0" "$1"_R2.fastq --output TAXONOMY/KRAKEN_gtdb/"$1"_kraken --report TAXONOMY/KRAKEN_gtdb/"$1"_gtdb_report"}' > kraken2_gtdb.sh
sh kraken2_gtdb.sh

module load braken/2.2
ls *_R1.fastq | awk -F '_R1.fastq' '{print "bracken -d /data/Food/analysis/R1111_Soilmetadata/Main_data/Kraken/kraken_gtdbtk_db/ -i TAXONOMY/KRAKEN_gtdb/"$1"_gtdb_report -o TAXONOMY/KRAKEN_gtdb/"$1"_gtdb_bracken.species.report "}' > bracken_custom.sh
sh bracken_custom.sh

module unload braken/2.2
module unload kraken2/2.1.1

find TAXONOMY/KRAKEN_gtdb -type f ! -name '*species.report' -delete

echo "**************************** Completed Running Custom Kraken ********************************************"
echo "*********************************************************************************************************"
echo "**************************** 4. Running Metaphlan 3.0 ***************************************************"

#RUN METAPHLAN3
mkdir TAXONOMY/METAPHLAN3
module load bbmap/38.22

ls *_R1.fastq | awk -F '_R1.fastq' '{print "reformat.sh in="$0" in="$1"_R2.fastq out=TAXONOMY/METAPHLAN3/"$0"_interleaved.fastq"}' > runreformat.sh
sh runreformat.sh

module unload bbmap/38.22
module load humann3/3.0

ls *_R1.fastq | awk -F '_R1.fastq' '{print "metaphlan TAXONOMY/METAPHLAN3/"$0"_interleaved.fastq --input_type fastq --bowtie2db /data/databases/MetaPhlAn3 -x mpa_v30_CHOCOPhlAn_201901 --bowtie2out TAXONOMY/METAPHLAN/"$0"_metaphlanbowtie2out.txt --nreads 10 --unknown_estimation > TAXONOMY/METAPHLAN3/"$0"_profile.txt "}' > metaphlan.sh
sh metaphlan.sh

module load metaphlan2/3.0
merge_metaphlan_tables.py TAXONOMY/METAPHLAN3/*_profile.txt > TAXONOMY/METAPHLAN/metaphlan.tsv
module unload metaphlan2/3.0

echo "**************************** Completed Running Metaphlan3 ***********************************************"
echo "*********************************************************************************************************"
echo "**************************** 5. Running Metaphlan4 ******************************************************"

cd TAXONOMY/METAPHLAN3

mkdir METAPHLAN4

module load humann3/3.6
ls -d *.fastq | awk -F '.fastq' '{print "/install/software/restart/py3/humann_3.5/bin/metaphlan "$0" --output_file METAPHLAN4/"$1"_mpa_out.txt --input_type fastq --bowtie2db /data/databases/MetaPhlan4 --index mpa_vJan21_CHOCOPhlAnSGB_202103 --nproc 10 --unclassified_estimation"}' > metaphlan.sh
sh metaphlan.sh

merge_metaphlan_tables.py METAPHLAN4/*_mpa_out.txt > METAPHLAN4/metaphlan4.txt

cd ..
mv METAPHLAN3/METAPHLAN4/ .

echo "**************************** Completed Running Metaphlan4 ***********************************************"
echo "*********************************************************************************************************"
echo "************************* Hopefully all went well for you! **********************************************"

mkdir Taxonomy_scripts_run
mv *.sh Taxonomy_scripts_run

cd ../..
mv insilicoseq/TrimmedFastQ/TAXONOMY .
