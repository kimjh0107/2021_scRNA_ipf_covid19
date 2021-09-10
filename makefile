GSE135893: GSE135893_barcodes.tsv.gz GSE135893_genes.tsv.gz

GSE135893_barcodes.tsv.gz: code/download_raw_gse.sh
	$< $@

GSE135893_genes.tsv.gz: code/download_raw_gse.sh
	$< $@

data/processed/IPF_SeuratObject.RDS: code/make_seurat_obj_ipf.R \
									data/GSE135893_IPF_1/GSE135893_seurat/ \
									data/IPF_metadata.csv
	$<
