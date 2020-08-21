ifndef PROFILE
override PROFILE = singularity
endif

test:
	nextflow -C modules/$(MODULE)/nextflow.config run modules/$(MODULE) \
	-entry test \
	-profile $(PROFILE) \
	--outdir output-test-module-$(MODULE)

test_contig_annot:
	nextflow run main.nf \
	-profile $(PROFILE) \
	--assembly "test_data/sample*.fasta" \
	--outdir output-test-ctg-annot

clean:
	rm -rf work output-test-* .nextflow*
