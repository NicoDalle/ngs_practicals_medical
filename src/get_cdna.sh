cd /home/rstudio/disk

#Télécharge une URL, -O : output
wget -O Hsap_cDNA.fa 'http://ensembl.org/biomart/martservice?query=<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "cdna" />
		<Attribute name = "external_gene_name" />
	</Dataset>
</Query>'
