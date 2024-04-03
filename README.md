# Benchmarking-taxonomic-classifiers-for-soil-shotgun-data

This repository hosts the final versions of scripts and pipelines used in our publication titled "An in-depth evaluation of metagenomic classifiers for soil microbiomes" It includes data preparation scripts for simulation, simulation scripts, and classifier scripts for processing the simulated sequencing data. Many scripts are modified versions of those originally created by Amy H Fitzpatrick, available at https://github.com/ahfitzpa/Benchmarking-bioinformatics-norovirus-amplicons

Classifiers compared within this study include Kaiju, Kraken, Kraken2 and Bracken with a custom GTDB-TK database, MetaPhlAn3, and MetaPhlAn4. The custom database used for Kraken incorporates the GTDB-TK database released in April 2022 and fungal genomes belonging to the phyla Basidiomycota, Mucoromycota, and Ascomycota, as available in the NCBI Assembly database as of January 20th, 2023.

**For Creating Soil-Specific Mock Community**

**<ins> Fungal-genome-downloads </ins>**

Fungal genomes were downloaded from NCBI's Assembly database using REntrez and command-line tools, based on criteria such as having a reference genome, representative genome, or chromosome/contig level assembly.

_Step 1_: To check how many of the fungal names obtained from MIAE have a Reference genome, Representative genome, or chromosome or contig level assembly.

_Step 2_: After obtaining the names that check-in for the filters used, obtain the Accession numbers belonging to each fungi

_Step 3_: Fetching metadata and double checking the pulls Once the metadata and Accession numbers are obtained use command line tools to obtain genome

_Step 4_: Download genomes using command line tools

Fungal-mitochondrial-genomes: To complement the collection, mitochondrial fasta files from the RefSeq database and the RefSoil database were downloaded, totaling 133 fungal representatives.

**<ins> Bacteria-Archaea-genome-downloads </ins>**

A comprehensive repository of bacterial and archaeal genome sequences was created by collating multiple sources, including the RefSoil database (https://github.com/germs-lab/ref_soil/tree/master/script_download_refsoil) and the NCBI Nuccore database, focusing on sequences from non-cultured soil bacteria with an isolation source described as soil. 

**<ins> Final metadata file </ins>**
A comprehensive metadata file, "SoilGenomeDB.xlsx", includes details of soil-specific bacterial, archaeal, and fungal genomes from NCBI.

**Manuscript**

For detailed insights and the scientific background of our project, please refer to our manuscript. You can access the manuscript through the following link: https://doi.org/10.1186/s40793-024-00561-w


