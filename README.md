# scVarCall

This is the pipeline used for calling MT variants in 10X scRNA


## Usage

to generate the variant calls from mitochondria we run the pipeline one of the
10X bam files.

```
bsub -Is \
    -R'select[mem>30000] rusage[mem=30000]' -M30000 -n 1 -R'span[hosts=1]' \
	go run scVarCall.go \
	-i /lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data/ti/cd/5892STDY8357359/possorted_genome_bam.bam \
	-o ../qc_filtered_scvarcall_out -b ../valid_barcode_list.txt
```


and merging the variant calls together into single objects

```
bsub -Is \
    -R'select[mem>30000] rusage[mem=30000]' -M30000 -n 1 -R'span[hosts=1]' \
    Rscript mergeChunkRds.R ../qc_filtered_scvarcall_out/chunk_*/
```
