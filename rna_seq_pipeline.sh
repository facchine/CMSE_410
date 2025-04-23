nextflow run nf-core/rnaseq -r 3.18.0 \
-profile docker \
-c /path/to/nextflow.config \
--input /path/to/samplesheet.csv \
--fasta /path/to/reference/genome \
--gtf /path/to/annotation \
--outdir /path/to/results \
--aligner star_salmon \
-with-report
