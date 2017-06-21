
test:
	haasbio/_tests/test.haasbio.py
	./haasbio/scripts/gtf_gene_db_indexer.py haasbio/test_data/small.gtf __testing
	./haasbio/scripts/gtf_gene_db_indexer.dumper.py __testing

clean:
	rm -f ./__testing.*

