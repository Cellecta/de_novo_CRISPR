python pipeline_v1.py  \
	./test_cds.fasta \
	/run/media/Dongfang/DataUSB/Cellecta_Project/denovo_crispr/db/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	/mnt/project/Pipeline/de_novo_CRISPR/Homo_sapiens.GRCh38.114.gtf \
	./test_3genes.xls --transcript_anno ./ccds_coding_gene_positions.tsv --number 20000 &
