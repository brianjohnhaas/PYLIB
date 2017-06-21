
test:
	haasbio/_tests/test.haasbio.py
	./haasbio/scripts/gtf_gene_db_indexer.py haasbio/test_data/small.gtf __testing
	./haasbio/scripts/gtf_gene_db_indexer.dumper.py __testing
	./haasbio/scripts/gtf_to_fasta.py --gtf haasbio/test_data/small.gtf --genome haasbio/test_data/small.fa --seq_type cDNA
	./haasbio/scripts/gtf_to_fasta.py --gtf haasbio/test_data/small.gtf --genome haasbio/test_data/small.fa --seq_type CDS


clean:
	rm -f ./__testing.*

