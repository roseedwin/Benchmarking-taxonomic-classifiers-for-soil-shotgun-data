#DOWNLOADING GENOMES USING NCBI COMMAND LINE TOOLS:

#On the mac terminal or linux terminal
#Instructions to install datasets and dataformats are available on NCBI website for command line tools
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

#Download datasets: 
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'

#Download dataformat: 
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'

#Make them executable:
chmod +x datasets dataformat

#First, create a conda environment: 
conda create -n ncbi_datasets

#Then, activate your new environment: 
conda activate ncbi_datasets

#Finally, install the datasets conda package: 
conda install -c conda-forge ncbi-datasets-cli

#After following instructions run the download_genomes.sh to create a bash script (download_genomes.sh) to loop the Accession numbers 
#I had cleaned the Accession_numbers.csv to only contain list of accession_numbers first (yes, manually)
bash download_genomes.sh

#Finally:
deactivate conda
