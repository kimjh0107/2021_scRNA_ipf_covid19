#!/bin/bash
gse=$1
echo $gse

wget --random-wait -r -p -e  robots=off  --span-hosts https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/$gse -O data/raw/$gse
# wget --random-wait -r -p -e  robots=off  --span-hosts https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/ -O raw/GSE135893_genes.tsv.gz
# wget --random-wait -r -p -e  robots=off  --span-hosts https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893_matrix.mtx.gz -O raw/GSE135893_matrix.mtx.gz
# wget --random-wait -r -p -e  robots=off  --span-hosts https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893_IPF_metadata.csv.gz -O raw/GSE135893_IPF_metadata.csv.gz