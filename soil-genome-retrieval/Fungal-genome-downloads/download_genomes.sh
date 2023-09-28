##DOWNLOADING FUNGAL GENOMES USING COMMAND LINE TOOLS
#Once obtained accession numbers now we can download using command line tools
for i in $(cat Accession_numbers.csv);do
datasets download genome accession "$i" --filename "$i".zip
done
