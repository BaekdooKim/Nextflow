/* 
 * RNA-Seq Tuxedo-v1 pipeline script
 *
 * @author
 *  Baekdoo Kim (baegi7942@gmail.com)
 */


def set_ch_annotation(annoFile) {
    Channel
        .fromPath(annoFile)
        .into {annotation_1; annotation_2; annotation_3; annotation_4}
}

def set_ch_bowtie_index(index_dir, index_name) {
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
}

def set_ch_genome(genomeFile) {
    Channel
        .fromPath(genomeFile)
        .into {genomes}
}

input_paired_1 = file(params.input_paired_1)
input_paired_2 = file(params.input_paired_2)
cuffmergeWrapper = params.cm_wrapper
annotationFile = file(params.input_gtf)

mate_std_dev = params.mate_inner_dist
anchor_length = params.anchor_length
segment_length = params.segment_length


if (params.run_split_input) {
    chunk_size = params.split_sequence_number
}

index_file = file(params.index)
index_name = index_file.getFileName()
index_dir = index_file.getParent()
genomeFile = file("${index_dir}/hg38.fa")
cuffmergeFile = file(cuffmergeWrapper)

index_files = file("${index_dir}/*.bt2")
index_files_count = index_files.size()

input_1_names = "${input_paired_1[0].getName()}, ${input_paired_1[1].getName()}"
input_2_names = "${input_paired_2[0].getName()}, ${input_paired_2[1].getName()}"
annotation_name = file(annotationFile).getName()
annotationPath = file(annotationFile).getParent()
genome_name = genomeFile.getName()

log.info "\n\n\t\t\t =====================================================\n\t\t\t ==  R N A - S e q   T U X E D O   P I P E L I N E  ==\n\t\t\t =====================================================\n\n"
log.info "\t===================================== I N F O M A T I O N ===================================="
log.info "\tInput pair 1                : ${input_1_names}"
log.info "\tInput pair 2                : ${input_2_names}"
log.info "\tAnnotation                  : ${annotation_name}"
log.info "\tGenome file                 : ${genome_name}"
log.info "\tIndex name                  : ${index_name}"
log.info "\tMate inner distance         : ${mate_std_dev}"
log.info "\tAnchor length               : ${anchor_length}"
log.info "\tRead segments length        : ${segment_length}"
log.info "\tSplitting input             : ${params.run_split_input}"
if (params.run_split_input) {
    log.info "\tSequence number to split    : ${chunk_size}"
}
log.info "\t==============================================================================================\n\n"

if (!annotationFile.exists()) {
    if (!params.download_missing_input) {
        exit 1, "- Annotation file not found (${annotationFile})"
    }
    else {
        log.info "[INFO] - Downloading Annotation file..\n"
        proc_anno = ['wget', 'http://146.95.173.35:9988/hg38/Annotation/genes.gtf', '-O', "${annotationFile}" ].execute()
        proc_anno.waitForProcessOutput()
        set_ch_annotation(annotationFile)
    }
} else {
    set_ch_annotation(annotationFile)
}
if (!genomeFile.exists()) {
    if (!params.download_missing_input) {
        exit 1, "- Genome file not found (${genomeFile})"
    }
    else {
        log.info "[INFO] - Downloading Genome file..\n"
        proc_gn = ['wget', 'http://146.95.173.35:9988/hg38/Sequence/Bowtie_2/hg38.fa', '-O', "${genomeFile}" ].execute()
        proc_gn.waitForProcessOutput()
        set_ch_genome(genomeFile)
    }
} else {
    set_ch_genome(genomeFile)
}

if (index_files_count < 6) {
    if (!params.download_missing_input) {
        exit 1, "- Please check your Bowtie2 index files (Required: 6 files but found only ${index_files_count})"
    }
    else {
        if (index_files_count > 0) {
            log.info "[INFO] - Moving old Bowtie2 index files to (${index_dir}_backup_old_indices)\n"
            proc_1 = ['mv', "${index_dir}", "${index_dir}_backup_old_indices"].execute()
            proc_1.waitForProcessOutput()
        }
        log.info "[INFO] - Downloading Bowtie2 index files..\n"
        proc_2 = ["wget", "-r", "-R", "index.html", "-N", "-nd", "-np", "http://146.95.173.35:9988/hg38/Sequence/Bowtie_2/", "-P", "${index_dir}"].execute()
        proc_2.waitForProcessOutput()
        set_ch_bowtie_index(index_dir,index_name)
        
    }
} else {
    set_ch_bowtie_index(index_dir,index_name)
}

Channel
    .fromPath(cuffmergeFile)
    .into {cuffmerge_wrapper}

if (params.run_split_input) {
    Channel
        .from(input_paired_1[0])
        .splitFastq(by: chunk_size, file: true)
        .set {read_set_1_1}

    Channel
        .from(input_paired_1[1])
        .splitFastq(by: chunk_size, file: true)
        .set {read_set_1_2}

    Channel
        .from(input_paired_2[0])
        .splitFastq(by: chunk_size, file: true)
        .set {read_set_2_1}

    Channel
        .from(input_paired_2[1])
        .splitFastq(by: chunk_size, file: true)
        .set {read_set_2_2}
}
else {
    Channel
        .from(input_paired_1[0])
        .set {read_set_1_1}

    Channel
        .from(input_paired_1[1])
        .set {read_set_1_2}

    Channel
        .from(input_paired_2[0])
        .set {read_set_2_1}

    Channel
        .from(input_paired_2[1])
        .set {read_set_2_2}
}

process Tophat_inputPair_1 {
    publishDir "$baseDir/output/splitted_tophat_files_1", mode: 'copy', overwrite: false
    input:
        val (index_name)
        set val (file_path), file(index_dir) from genome_index.first()
        file (readFile_1_1) from read_set_1_1
        file (readFile_1_2) from read_set_1_2

    output:
        file ("accepted_hits_${readFile_1_1}.bam") into tophat_out_1

    script:
        indexPath = file_path.getParent()
        """
        tophat2 -p 8 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} $indexPath/${index_name} ${readFile_1_1} ${readFile_1_2}
        mv tophat_out/accepted_hits.bam accepted_hits_${readFile_1_1}.bam
        """
}

process Tophat_inputPair_2 {
    publishDir "$baseDir/output/splitted_tophat_files_2", mode: 'copy', overwrite: false
    input:
        val (index_name)
        set val (file_path), file(index_dir) from genome_index.first()
        file (readFile_2_1) from read_set_2_1
        file (readFile_2_2) from read_set_2_2

    output:
        file ("accepted_hits_${readFile_2_1}.bam") into tophat_out_2

    script:
        indexPath = file_path.getParent()
        """
        tophat2 -p 8 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} $indexPath/${index_name} ${readFile_2_1} ${readFile_2_2}
        mv tophat_out/accepted_hits.bam accepted_hits_${readFile_2_1}.bam
        """
}


if (params.run_split_input) {
    process Merge_Tophat_output_1 {

        input:
            file filePath_1 from tophat_out_1.toList()

        output:
            file ("merged_bam_1.bam") into tophat_merged_out_1
        
        shell:
            '''
            unsorted=(!{filePath_1})
            for i in $(printf "%s\n" "${unsorted[@]}" | sort)
            do
                echo $i
                cat $i >> merged_bam_1.bam
            done
            '''
    }

    process Merge_Tophat_output_2 {

        input:
            file filePath_2 from tophat_out_2.toList()

        output:
            file ("merged_bam_2.bam") into tophat_merged_out_2
        
        shell:
            '''
            unsorted=(!{filePath_2})
            for i in $(printf "%s\n" "${unsorted[@]}" | sort);
            do
                echo $i
                cat $i >> merged_bam_2.bam
            done
            '''
    }

    tophat_merged_out_1
        .into { merged_tophat_output_1_1; merged_tophat_output_1_2 }
    tophat_merged_out_2
        .into { merged_tophat_output_2_1; merged_tophat_output_2_2 }
}
else {
    tophat_out_1
        .into { merged_tophat_output_1_1; merged_tophat_output_1_2 }
    tophat_out_2
        .into { merged_tophat_output_2_1; merged_tophat_output_2_2 }
}

process BamToSam_1 {
    publishDir "$baseDir/output/bam2sam", mode: 'copy', overwrite: false

    input:
        file tophat_output1 from merged_tophat_output_1_1

    output:
        file "BamToSam_1.sam" into bam2sam_out_1

    """
    mkdir -p ${baseDir}/output/bam2sam
    samtools view -o BamToSam_1.sam -h ${tophat_output1}
    """
}

process BamToSam_2 {
    publishDir "$baseDir/output/bam2sam", mode: 'copy', overwrite: false

    input:
        file tophat_output2 from merged_tophat_output_2_1

    output:
        file "BamToSam_2.sam" into bam2sam_out_2

    """
    mkdir -p ${baseDir}/output/bam2sam
    samtools view -o BamToSam_2.sam -h ${tophat_output2}
    """
}

process Cufflinks_1 {
    publishDir "$baseDir/output/cufflinks", mode: 'copy', overwrite: false

    input:
        file tophat_output1 from merged_tophat_output_1_2
        file anno_file from annotation_1

    output:
        file "transcripts_1.gtf" into cufflinks_gtf_out_1
    """
    cufflinks -q --no-update-check -p 8 -I 300000 -F 0.1 -j 0.15 -G ${anno_file} ${tophat_output1}
    mv transcripts.gtf transcripts_1.gtf
    """
}

process Cufflinks_2 {
    publishDir "$baseDir/output/cufflinks", mode: 'copy', overwrite: false

    input:
        file tophat_output2 from merged_tophat_output_2_2
        file anno_file from annotation_2

    output:
        file "transcripts_2.gtf" into cufflinks_gtf_out_2
    """
    cufflinks -q --no-update-check -p 8 -I 300000 -F 0.1 -j 0.15 -G ${anno_file} ${tophat_output2}
    mv transcripts.gtf transcripts_2.gtf
    """
}

process Cuffmerge {
    publishDir "$baseDir/output/cuffmerge", mode: 'copy', overwrite: false

    input:
        file cufflinks_gtf_1 from cufflinks_gtf_out_1
        file cufflinks_gtf_2 from cufflinks_gtf_out_2
        file anno_file from annotation_3
        file cuffmergeWrapper from cuffmerge_wrapper

    output:
        file "merged_transcripts.gtf" into cuffmerge_gtf_out

    """
    python ${cuffmergeWrapper} -p 8 -g ${anno_file} --min-isoform-fraction="0.05" --merged-transcripts="merged_transcripts.gtf" ${cufflinks_gtf_1} ${cufflinks_gtf_2}
    """
}

process Cuffcompare {
    publishDir "$baseDir/output/cuffcompare", mode: 'copy', overwrite: false

    input:
        file genome from genomes
        file anno_file from annotation_4
        file cuffmerge_gtf from cuffmerge_gtf_out

    output:
        file "cuffcmp.combined.gtf" into cuffcompare_out

    """
    cuffcompare -r ${anno_file} -s ${genome} -e 100 -d 100 ${cuffmerge_gtf}
    """
}

process Cuffdiff {
    publishDir "$baseDir/output/cuffdiff", mode: 'copy', overwrite: false

    input:
        file cuffcompare_gtf from cuffcompare_out
        file bam2sam_sam_1 from bam2sam_out_1
        file bam2sam_sam_2 from bam2sam_out_2

    output:
        file "*.diff" into cuffdiff_out
    """
    cuffdiff --no-update-check --FDR=0.05 -p 8 --min-alignment-count=10 --library-norm-method=geometric --dispersion-method=pooled --labels "," ${cuffcompare_gtf} ${bam2sam_sam_1} ${bam2sam_sam_2}
    """
}


workflow.onComplete {
    println(workflow.success ? "\n\t Final output: ./output/cuffdiff/gene_exp.diff\n\t Done! ": "\t Oops .. something went wrong")
}

