# Benchmarking-taxonomic-classifiers-for-soil-shotgun-data

This repository contains the final version of scripts used to prepare data for simulation, simulation scripts, alongside the prepared pipelines and classifier scripts used to process the simulated sequencing data for the publication on Benchmarking taxonomic classifiers for accurate profiling of fungal, bacterial and archaeal communities in soil shotgun metagenomic data using in-silico tools. Many of the scripts used are modified versions of scripts previously made by Amy H Fitzpatrick available at https://github.com/ahfitzpa/Benchmarking-bioinformatics-norovirus-amplicons

Classifiers compared include Kaiju, Kraken, and Kraken with custom gtdb-tk database, Metaphlan3 and Metaphlan4.

And custom database used for Kraken is the GTDB-TK database released in 2022 April and Fungal genomes belonging to phylum Basidiomycota, Mucoromycota and Ascomycota as available in the NCBI assembly database on 20th January 2023

**For Creating Soil-Specific Mock Community**

**<ins> Fungi-genome-downloads </ins>**

Using REntrez and command line tools download fungal genomes from the Assembly database of NCBI

_Step 1_: To check how many of the fungal names obtained from MIAE have a Reference genome, Representative genome, or chromosome or contig level assembly.

_Step 2_: After obtaining the names that check-in for the filters used, obtain the Accession numbers belonging to each fungi

_Step 3_: Fetching metadata and double checking the pulls Once the metadata and Accession numbers are obtained use command line tools to obtain genome

_Step 4_: Download genomes using command line tools

Fungal-mitochondrial-genomes: In the absence of a huge turnover of Fungal genomes (89 genomes obtained from Assembly database). The mitochondrial fasta files, in RefSeq from the nuccore database were also downloaded to add to the database (38 mitochondrial genomes) and 6 from the RefSoil database

A total of 133 fungal representatives for the soil Fungal database

**<ins> Bacteria-Archaea-genome-downloads </ins>**

To create a large repository of genome sequences belonging to bacteria and archaea, multiple sources were collated.

RefSoil database: The genome sequences and metadata were obtained from https://github.com/germs-lab/ref_soil/tree/master/script_download_refsoil
Since this contained only sequences of cultivated soil species we downloaded more sequences from the NCBI database to obtain sequences of non-cultured soil bacteria.

NCBI Nuccore database: Acquiring all Bacterial complete genome sequences from the 'Nuccore' database in NCBI from RefSeq database whose isolation source is described as soil.
