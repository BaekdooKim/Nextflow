
input_paired = file(params.input_paired)
index_file = file(params.index)
index_name = index_file.getFileName()
index_dir = index_file.getParent()
index_files = file("${index_dir}/*.bt2")
genomeFile = "${params.index}.fa"
macs2_wrapper = file(params.macs2_wrapper)
annotationFile = params.input_gtf
bowtie2_index = params.index
max_frag_length = params.max_frag_length
pvalue = params.pvalue
gsize = params.gsize
bwidth = params.bwidth

input_names = "${input_paired[0].getName()}, ${input_paired[1].getName()}"
annotation_name = file(annotationFile).getName()
genome_name = file(genomeFile).getName()


Channel
    .fromPath(genomeFile)
    .into {genomes}

Channel
    .from(macs2_wrapper)
    .set{macs2_wrapper_ch}

log.info "\n\n\t ===========================================================\n\t ==  C h I P - S e q ( P a i r - E n d ) P I P E L I N E  ==\n\t ===========================================================\n\n"
log.info " ===================================== I N F O M A T I O N ===================================="
log.info " Input pair               : ${input_names}"
log.info " Annotation               : ${annotation_name}"
log.info " Genome file              : ${genome_name}"
log.info " Maximum fragment length  : ${max_frag_length}"
log.info " P-value                  : ${pvalue}"
log.info " Effective genome size    : ${gsize}"
log.info " Band width               : ${bwidth}"
log.info " =============================================================================================="

input_ext = input_paired[0].getExtension()

process Index_setting {
    executor 'local'
    input:
        file(index_dir)
        val(index_name)
    output:
        file("genomeIndex/*") into genome_index

    script:
        """
        mkdir -p genomeIndex
        cp ${index_dir}/${index_name}.* genomeIndex/
        """
}

if (input_ext == 'gz') {
    process gunzip_1 {
        output:
            file 'input_1.fastq' into input_1
        """
            gunzip -c ${input_paired[0]} > input_1.fastq
        """
    }
    process gunzip_2 {
        output:
            file 'input_2.fastq' into input_2
        """
            gunzip -c ${input_paired[1]} > input_2.fastq
        """
    }
    process Bowtie2 {
        publishDir "$baseDir/output/bowtie", mode: 'copy', overwrite: false
        input:
            file input_1 from input_1
            file input_2 from input_2
            val (index_name)
            set val (file_path), file(index_dir) from genome_index.first()
        output:
            file "bowtie_out.bam" into bowtie_out

        script:
        indexPath = file_path.getParent()
        """
        bowtie2 -p 8 -x $indexPath/${index_name} -1 ${input_1} -2 ${input_2} -X ${max_frag_length} -N "0" -L "22" -i "S,1,1.15" --n-ceil "L,0,0.15" --dpad "15" --gbar "4" --end-to-end --score-min "L,-0.6,-0.6" | samtools view -Su - | samtools sort -o - - > bowtie_out.bam
        """
    }
}
else {
    Channel
        .from(input_paired)
        .set {chip_seq_input}
    process Bowtie2 {
        publishDir "$baseDir/output/bowtie", mode: 'copy', overwrite: false
        input:
            file genome from genomeFile
            file fin from chip_seq_input.collect()
            val (index_name)
            set val (file_path), file(index_dir) from genome_index.first()
        output:
            file "bowtie_out.bam" into bowtie_out

        script:
        indexPath = file_path.getParent()
        """
        bowtie2 -p 8 -x $indexPath/${index_name} -1 ${fin[0]} -2 ${fin[1]} -X ${max_frag_length} -N "0" -L "22" -i "S,1,1.15" --n-ceil "L,0,0.15" --dpad "15" --gbar "4" --end-to-end --score-min "L,-0.6,-0.6" | samtools view -Su - | samtools sort -o - - > bowtie_out.bam
        """
    }
}

bowtie_out.into{bowtie_out_1;bowtie_out_2}

process BamToSam {
    publishDir "$baseDir/output/bam2sam", mode: 'copy', overwrite: false

    input:
        file bowtie_output from bowtie_out

    output:
        file "bam2sam_out.sam" into bam2sam_out_1

    """
    samtools view -o bam2sam_out.sam -h ${bowtie_output}
    """
}

process Macs2 {
    publishDir "${baseDir}/output/macs2", mode: 'copy', overwrite: false

    input:
        set val (bowtie_output), file (files) from bowtie_out_2
        file macs2 from macs2_wrapper_ch
    output:
        file 'macs2_*' into macs2_out
    
    """
    python ${macs2} callpeak ${bowtie_output} --format=BAMPE --name="MACS2.1.0 in Nextflow" --bw=${bwidth} --gsize=${gsize} --pvalue=${pvalue} --mfold 10 30 --keep-dup 1 --output-summits=macs2_result.bed --output-extra-files=macs2_report.html --output-extra-files-path=macs2_report --output-narrowpeaks=macs2_narrowpeaks.interval --output-xls-to-interval=macs2_xls_to_interval.interval
    """
}

workflow.onComplete {
    println(workflow.success ? "\n\t Output path: ./output/macs2\n\t Done! ": "\t Oops .. Something went wrong!")
}
