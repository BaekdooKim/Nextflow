# ChIP-Seq (paired-end)


* It is a simple ChIP-Seq Nextflow pipeline.
* To run this pipeline, you need to provide...
	1. Paired-End fastq input
	2. Annotation file (one gtf file)
	3. Bowtie2 index files

## Nextflow.config
```
input_paired = "<YOUR-INPUT-PATH>/*_{1,2,R1,R2}*.fastq*"
input_gtf = "<YOUR-ANNOTATION-PATH>/<NAME>.gtf"
index = "<YOUR-BOWTIE2-INDEX-PATH>/<PREFIX>"
```

## Tools options
1. Maximum fragment length (default 500)
2. p-value (default 0.05)
3. Effective genome size (default 2700000000)
4. Band width (default 200)

## Run the pipeline
```
$ git clone https://github.com/BaekdooKim/Nextflow-pipelines.git
$ cd ChIP-Seq_paired-end

## Run the pipeline using default option values
$ nextflow nextflow main.nf

## Otherwise..
$ nextflow nextflow main.nf --max_frag_length=<VALUE> --pvalue=<VALUE> --gsize=<VALUE> --bwidth=<VALUE>
```