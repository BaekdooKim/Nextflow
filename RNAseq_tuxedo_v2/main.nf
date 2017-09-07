/* 
 * RNA-Seq Tuxedo-v2 pipeline script
 *
 *  Author
 *  - Baekdoo Kim (baegi7942@gmail.com)
 *  - Aug 29, 2017
 */


def set_ch_annotation(annoFile) {
    Channel
        .fromPath(annoFile)
        .into {annotation_1; annotation_2; annotation_3; annotation_4}
}

def set_ch_Hisat2_index(index_dir, index_name) {
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
prepDE = params.prepDE
gen_input_list = params.gen_input_list
annotationFile = file(params.input_gtf)

if (params.run_split_input) {
    chunk_size = params.split_sequence_number
}

index_file = file(params.index)
index_name = index_file.getFileName()
index_dir = index_file.getParent()
genomeFile = file("${index_dir}/hg38.fa")
cuffmergeFile = file(cuffmergeWrapper)
prepDE_File = file(prepDE)
genInputListFile = file(gen_input_list)

ref_interacts_file = file(params.ref_interactions)
replace_char_file = file(params.replace_char)
convert_chars_file = file(params.convert_chars)
remove_begin_file = file(params.remove_beginning)
assem_filtering_file = file(params.assembled_filtering)
join_tool_file = file(params.join_tool)
cut_tool_file = file(params.cut_tool)
table_arithmetic_file = file(params.arithmetic_operations)
paste_wrapper_file = file(params.paste_wrapper)
fixed_value_column_file = file(params.fixed_value_column)
merge_column_file = file(params.merge_column)
easyjoin_file = file(params.easyjoin)
correlation_file = file(params.correlation)
tabular2HTML_file = file(params.tabular2HTML)
gen_header_detail_file = file(params.gen_header_detail)
html_header_file = file(params.html_header)


index_files = file("${index_dir}/*.ht2")
index_files_count = index_files.size()

input_1_names = "${input_paired_1[0].getName()}, ${input_paired_1[1].getName()}"
input_2_names = "${input_paired_2[0].getName()}, ${input_paired_2[1].getName()}"
annotation_name = file(annotationFile).getName()
annotationPath = file(annotationFile).getParent()
genome_name = genomeFile.getName()

min_ins = params.min_ins
max_ins = params.max_ins
read_len = params.read_len

log.info "\n\n\t\t\t ============================================================\n\t\t\t ==  R N A - S e q   T U X E D O - ver.2  P I P E L I N E  ==\n\t\t\t ============================================================\n\n"
log.info "\t===================================== I N F O M A T I O N ===================================="
log.info "\tInput pair 1                : ${input_1_names}"
log.info "\tInput pair 2                : ${input_2_names}"
log.info "\tAnnotation                  : ${annotation_name}"
log.info "\tGenome file                 : ${genome_name}"
log.info "\tIndex name                  : ${index_name}"
log.info "\tMin fragment length         : ${min_ins}"
log.info "\tMax fragment length         : ${max_ins}"
log.info "\tThe average read length     : ${read_len}"
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
        proc_gn = ['wget', 'http://146.95.173.35:9988/hg38/Sequence/Hisat_2/hg38.fa', '-O', "${genomeFile}" ].execute()
        proc_gn.waitForProcessOutput()
        set_ch_genome(genomeFile)
    }
} else {
    set_ch_genome(genomeFile)
}

if (index_files_count != 8) {
    if (!params.download_missing_input) {
        exit 1, "- Please check your Hisat2 index files (Required: 8 hisat index files but found ${index_files_count})"
    }
    else {
        if (index_files_count > 0) {
            log.info "[INFO] - Moving old Hisat2 index files to (${index_dir}_backup_old_indices)\n"
            proc_1 = ['mv', "${index_dir}", "${index_dir}_backup_old_indices"].execute()
            proc_1.waitForProcessOutput()
        }
        log.info "[INFO] - Downloading Hisat2 index files..\n"
        proc_2 = ["wget", "-r", "-R", "index.html", "-N", "-nd", "-np", "http://146.95.173.35:9988/hg38/Sequence/Hisat_2/", "-P", "${index_dir}"].execute()
        proc_2.waitForProcessOutput()
        set_ch_Hisat2_index(index_dir,index_name)
        
    }
} else {
    set_ch_Hisat2_index(index_dir,index_name)
}

Channel
    .fromPath(cuffmergeFile)
    .into {cuffmerge_wrapper}

Channel
    .fromPath(prepDE_File)
    .into {prepDE_script}

Channel
    .fromPath(replace_char_file)
    .into {replace_char_ch}

Channel
    .fromPath(ref_interacts_file)
    .into {ref_interacts_ch}

Channel
    .fromPath(convert_chars_file)
    .into {convert_chars_ch}

Channel
    .fromPath (remove_begin_file)
    .into {remove_begin_file_ch}

Channel
    .fromPath (assem_filtering_file)
    .into {assem_filtering_file_ch}

Channel
    .fromPath (join_tool_file)
    .into {join_tool_ch}

Channel
    .fromPath (cut_tool_file)
    .into {cut_tool_ch}

Channel
    .fromPath (table_arithmetic_file)
    .into {table_arithmetic_ch}

Channel
    .fromPath(paste_wrapper_file)
    .into {paste_wrapper_ch}

Channel
    .fromPath(fixed_value_column_file)
    .into {fixed_value_column_ch}

Channel
    .fromPath(merge_column_file)
    .into {merge_column_ch}

Channel
    .fromPath(easyjoin_file)
    .into {easyjoin_ch}

Channel
    .fromPath(correlation_file)
    .into {correlation_ch}

Channel
    .fromPath(tabular2HTML_file)
    .into {tabular2HTML_ch}

Channel
    .fromPath(gen_header_detail_file)
    .into {gen_header_detail_ch}

Channel
    .fromPath(html_header_file)
    .into {html_header_ch}


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

process Hisat2_inputPair_1 {
    publishDir "$baseDir/output/splitted_hisat2_files_1", mode: 'copy', overwrite: false
    input:
        val (min_ins)
        val (max_ins)
        val (index_name)
        set val (file_path), file(index_dir) from genome_index.first()
        file (readFile_1_1) from read_set_1_1
        file (readFile_1_2) from read_set_1_2

    output:
        file ("hisat_aligned_${readFile_1_1}.sam") into hisat_out_1

    script:
        indexPath = file_path.getParent()
        """
        hisat2 -p 8 --dta --secondary -I ${min_ins} -X ${max_ins} -x $indexPath/${index_name} -1 ${readFile_1_1} -2 ${readFile_1_2} -S hisat_aligned_${readFile_1_1}.sam
        """
}

process Hisat2_inputPair_2 {
    publishDir "$baseDir/output/splitted_hisat2_files_2", mode: 'copy', overwrite: false
    input:
        val (min_ins)
        val (max_ins)
        val (index_name)
        set val (file_path), file(index_dir) from genome_index.first()
        file (readFile_2_1) from read_set_2_1
        file (readFile_2_2) from read_set_2_2

    output:
        file ("hisat_aligned_${readFile_2_1}.sam") into hisat_out_2

    script:
        indexPath = file_path.getParent()
        """
        hisat2 -p 8 --dta --secondary -I ${min_ins} -X ${max_ins} -x $indexPath/${index_name} -1 ${readFile_2_1} -2 ${readFile_2_2} -S hisat_aligned_${readFile_2_1}.sam
        """
}

if (params.run_split_input) {
    process Merge_Hisat2_output_1 {

        input:
            file filePath from hisat_out_1.toList()

        output:
            file ("merged_sam_1.sam") into hisat_merged_sam_1

        shell:
            '''
            unsorted=(!{filePath})
            header=${unsorted[0]}
            (grep ^@ $header; for f in $(printf "%s\n" "${unsorted[@]}" | sort -V); do grep -v ^@ $f; done) > merged_sam_1.sam
            '''
    }


    process Merge_Hisat2_output_2 {

        input:
            file filePath from hisat_out_2.toList()

        output:
            file ("merged_sam_2.sam") into hisat_merged_sam_2

        shell:
            '''
            unsorted=(!{filePath})
            header=${unsorted[0]}
            (grep ^@ $header; for f in $(printf "%s\n" "${unsorted[@]}" | sort -V); do grep -v ^@ $f; done) > merged_sam_2.sam
            '''
    }

    hisat_merged_sam_1
        .into { hisat_merged_sam_1_1; hisat_merged_sam_1_2 }
    hisat_merged_sam_2
        .into { hisat_merged_sam_2_1; hisat_merged_sam_2_2 }
}
else {
    hisat_out_1
        .into { hisat_merged_sam_1_1; hisat_merged_sam_1_2 }
    hisat_out_2
        .into { hisat_merged_sam_2_1; hisat_merged_sam_2_2 }
}

process SamToBam_1 {
    publishDir "$baseDir/output/sam2bam", mode: 'copy', overwrite: false

    input:
        file hisat_output from hisat_merged_sam_1_1

    output:
        file "hisat_aligned_1.bam" into sam2bam_out_1

    """
    samtools sort -@ 8 -o hisat_aligned_1.bam ${hisat_output}
    """
}

process SamToBam_2 {
    publishDir "$baseDir/output/sam2bam", mode: 'copy', overwrite: false

    input:
        file hisat_output from hisat_merged_sam_2_1

    output:
        file "hisat_aligned_2.bam" into sam2bam_out_2

    """
    samtools sort -@ 8 -o hisat_aligned_2.bam ${hisat_output}
    """
}

sam2bam_out_1
    .into { sam2bam_out_1_1; sam2bam_out_1_2 }
sam2bam_out_2
    .into { sam2bam_out_2_1; sam2bam_out_2_2 }

process Stringtie_1 {
    publishDir "$baseDir/output/stringtie", mode: 'copy', overwrite: false

    input:
        file hisat_output from sam2bam_out_1_1
        file anno_file from annotation_1

    output:
        file "stringtie_assembled_1.gtf" into stringtie_gtf_out_1

    """
    stringtie -p 8 -e -G ${anno_file} -o stringtie_assembled_1.gtf ${hisat_output}
    """
}

stringtie_gtf_out_1
    .into { stringtie_gtf_out_1_1; stringtie_gtf_out_1_2 }

process Stringtie_2 {
    publishDir "$baseDir/output/stringtie", mode: 'copy', overwrite: false

    input:
        file hisat_output from sam2bam_out_2_1
        file anno_file from annotation_2

    output:
        file "stringtie_assembled_2.gtf" into stringtie_gtf_out_2

    """
    stringtie -p 8 -e -G ${anno_file} -o stringtie_assembled_2.gtf ${hisat_output}
    """
}

process Generate_GeneCount {
    publishDir "$baseDir/output/gene_count_matrix", mode: 'copy', overwrite: false

    input:
        val (read_len)
        file stringtie_1 from stringtie_gtf_out_1_1
        file stringtie_2 from stringtie_gtf_out_2
        file gen_input_list from genInputListFile
        file prepDE from prepDE_script

    output:
        file "gene_count_matrix.csv" into gene_count_matrix
        file "transcript_count_matrix.csv" into transcript_count_matrix

    """
    sh ${gen_input_list}
    python ${prepDE} -l ${read_len} -i input_list.txt
    """
}

gene_count_matrix
    .into { gene_count_matrix_1; gene_count_matrix_2 }

process DESeq2 {
    publishDir "$baseDir/output/DESeq2", mode:"copy", overwrite: false

    input:
        file gene_count from gene_count_matrix_1
        file transcript_count from transcript_count_matrix

    output:
        file "genes.diff" into deseq2_output

    shell:
    '''
    #!/usr/bin/env Rscript

    library("DESeq2")
    countData <- as.matrix(read.csv("!{gene_count}", row.names="gene_id"))
    sampleInfo = data.frame(  row.names = colnames( countData ), condition = c( "normal", "disease"))
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleInfo, design = ~condition)
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered = res[order(res$padj),]
    write.table(resOrdered, file="genes.diff", sep="\t")
    '''
}

process Generate_Report {
    publishDir "$baseDir/output/Report", mode:"copy", overwrite: false

    input:
        file genes from deseq2_output
        file gene_count from gene_count_matrix_2
        file stringtie_assembled from stringtie_gtf_out_1_2
        file replace_chr from replace_char_ch
        file ref_interactions from ref_interacts_ch
        file conv_chars from convert_chars_ch
        file remove_beginning from remove_begin_file_ch
        file assembled_filtering from assem_filtering_file_ch
        file join_tool from join_tool_ch
        file cut_tool from cut_tool_ch
        file table_arithmetic from table_arithmetic_ch
        file paste_wrapper from paste_wrapper_ch
        file fixed_value_column from fixed_value_column_ch
        file merge_column from merge_column_ch
        file easyjoin from easyjoin_ch
        file correlation from correlation_ch
        file tabular2HTML from tabular2HTML_ch
        file gen_header_detail from gen_header_detail_ch
        file html_header from html_header_ch

    output:
        file "*_Genes.html" into reports
        
    shell:
        '''
        sh !{replace_chr} !{genes} replace_chr.diff
        python !{conv_chars} --strip --condense !{gene_count} C convert_chars_geneCount.tabular
        python !{conv_chars} --strip --condense !{stringtie_assembled} Sc convert_chars_stringtie.tabular
        python !{assembled_filtering} convert_chars_stringtie.tabular filtered_stringtie_assembled.tabular "c3==__sq__transcript__sq__" 15 "str,str,str,int,int,str,str,str,str,str,str,str,str,str,str" 0
        perl !{remove_beginning} replace_chr.diff 1 remove_beginning_genes.tabular
        perl !{remove_beginning} convert_chars_geneCount.tabular 1 remove_beginning_geneCount.tabular
        Rscript !{join_tool} --out="join_output.tabular" --he="0" --jc="1" --sep="ta" --nc="0" --in="remove_beginning_genes.tabular" --in="remove_beginning_geneCount.tabular"
        perl !{cut_tool} filtered_stringtie_assembled.tabular "c9" "T" cut_c9.tabular
        sh !{replace_chr} cut_c9.tabular cut_c9_replaced.tabular
        sed -i "s/gene_id //g" cut_c9_replaced.tabular
        perl !{cut_tool} filtered_stringtie_assembled.tabular "c1,c4,c5" "T" cut_c1_c4_c5.tabular
        perl !{cut_tool} cut_c1_c4_c5.tabular "c2" "T" cut_c1_c4_c5_then_c2.tabular
        perl !{cut_tool} cut_c1_c4_c5.tabular "c3" "T" cut_c1_c4_c5_then_c3.tabular
        perl !{table_arithmetic} cut_c1_c4_c5_then_c3.tabular cut_c1_c4_c5_then_c2.tabular Subtraction table_arithmetic_out.tabular
        perl !{paste_wrapper} cut_c1_c4_c5.tabular table_arithmetic_out.tabular T paste_wrapper_out_1.tabular
        perl !{fixed_value_column} paste_wrapper_out_1.tabular add_value_column_out_1.tabular ":" "no"
        perl !{fixed_value_column} add_value_column_out_1.tabular add_value_column_out_2.tabular "-" "no"
        python !{merge_column} add_value_column_out_2.tabular merged_column_out.tabular "1" "5" 2 6 3
        perl !{cut_tool} merged_column_out.tabular "c4,c7" "T" cut_merged_column_c4_c7.tabular
        perl !{paste_wrapper} cut_c9_replaced.tabular cut_merged_column_c4_c7.tabular T paste_wrapper_out_2.tabular
        perl !{easyjoin} -a 1 -t$"\t"  -e "N/A" -o auto  -1 "1" -2 "1" join_output.tabular paste_wrapper_out_2.tabular > easyjoin_output.tabular
        sort -u -t$"\t" -o unique_lines.tabular easyjoin_output.tabular
        perl !{cut_tool} unique_lines.tabular "c1,c11,c10,c8,c9,c2,c3,c6" "T" cut_c1_c11_c10_c8_c9_c2_c3_c6.tabular
        
        perl !{easyjoin} -a 1 -t$"\t"  -e "N/A" -o auto  -1 "1" -2 "1" cut_c1_c11_c10_c8_c9_c2_c3_c6.tabular !{ref_interactions} > easyjoin_cut_ref.tabular
        
        python !{correlation} easyjoin_cut_ref.tabular cor_easyjoin.tabular 4,5 pearson
        tail --lines 1 cor_easyjoin.tabular > cor_easyjoin_selected_last.tabular
        perl !{cut_tool} cor_easyjoin_selected_last.tabular "c1" "T" cor_ej_sel_cut_c1.tabular
        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Pearson Correlation/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header cor_ej_sel_cut_c1.tabular concatenate_dataset_2.html
        wc -l easyjoin_cut_ref.tabular | awk '{ print $1 }' >> line_word_char_count.tabular
        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Total Number of Differentially Expressed Genes/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header line_word_char_count.tabular concatenate_dataset_0.html
        perl !{fixed_value_column} easyjoin_cut_ref.tabular easyjoin_cut_ref_add_val.tabular "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38XlastVirtModeType=defaultXlastVirtModeExtraState=XvirtModeType=defaultXvirtMode=0XnonVirtPosition=Xposition=" "no"
        python !{merge_column} easyjoin_cut_ref_add_val.tabular easyjoin_cut_ref_add_val_merged.tabular "10" "2"
        perl !{fixed_value_column} easyjoin_cut_ref_add_val_merged.tabular easyjoin_cut_ref_add_val_merged_added.tabular "http://string-db.org/cgi/network.pl?identifier=" "no"
        python !{merge_column} easyjoin_cut_ref_add_val_merged_added.tabular easyjoin_cut_ref_add_val_merged_added_merged.tabular "12" "1"
        perl !{fixed_value_column} easyjoin_cut_ref_add_val_merged_added_merged.tabular easyjoin_cut_ref_add_val_merged_added_merged_added.tabular "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" "no"
        python !{merge_column} easyjoin_cut_ref_add_val_merged_added_merged_added.tabular easyjoin_cut_ref_add_val_merged_added_merged_added_merged.tabular "14" "1"
        perl !{cut_tool} easyjoin_cut_ref_add_val_merged_added_merged_added_merged.tabular "c1,c2,c3,c11,c9,c13,c15,c4,c5,c6,c7,c8" "T" multiple_added_merged_cut.tabular
        python !{assembled_filtering} multiple_added_merged_cut.tabular multiple_added_merged_cut_filtered.tabular "c12__lt__=0.05" 12 "str,str,int,str,str,str,str,int,int,float,float,float" 0
        python !{assembled_filtering} multiple_added_merged_cut_filtered.tabular multiple_added_merged_cut_filtered_not_c11.tabular "c11!=0" 12 "str,str,int,str,str,str,str,int,int,float,float,float" 0
        wc -l multiple_added_merged_cut_filtered_not_c11.tabular | awk '{ print $1 }' >> multiple_added_merged_cut_filtered_not_c11_count.tabular
        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Total Significant Genes/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header multiple_added_merged_cut_filtered_not_c11_count.tabular concatenate_dataset_1.html 1

        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Significant Genes/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header
        perl !{tabular2HTML} moded_header multiple_added_merged_cut_filtered_not_c11.tabular concatenate_dataset_5.html 12

        ( LC_ALL=C  sort   --stable -t$"\t"  -k "11n,11"  ) < multiple_added_merged_cut_filtered_not_c11.tabular > multiple_added_merged_cut_filtered_not_c11_sorted.tabular 
        head --lines 1 multiple_added_merged_cut_filtered_not_c11_sorted.tabular > not_c11_sorted_head1.tabular
        tail --lines 1 multiple_added_merged_cut_filtered_not_c11_sorted.tabular > not_c11_sorted_tail1.tabular 
        
        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Gene With the Lowest Change/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header
        perl !{tabular2HTML} moded_header not_c11_sorted_head1.tabular concatenate_dataset_3.html 12
        
        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Gene With the Highest Change:/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header not_c11_sorted_head1.tabular concatenate_dataset_4.html 12

        cat concatenate_dataset_0.html concatenate_dataset_1.html concatenate_dataset_2.html concatenate_dataset_3.html concatenate_dataset_4.html concatenate_dataset_5.html > Total_Differentially_Expressed_Genes.html 


        python !{assembled_filtering} multiple_added_merged_cut_filtered.tabular multiple_added_merged_cut_filtered_c11.tabular "c11==0" 12 "str,str,int,str,str,str,str,int,int,float,float,float" 0
        wc -l multiple_added_merged_cut_filtered_c11.tabular | awk '{ print $1 }' >> filtered_c11_lwc_count.tabular
        
        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Total Undifferentially Expressed Genes/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header filtered_c11_lwc_count.tabular concatenate_dataset_2_1.html 1

        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Genes/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header multiple_added_merged_cut_filtered_c11.tabular concatenate_dataset_2_2.html 12
        cat concatenate_dataset_2_1.html concatenate_dataset_2_2.html > Undifferentially_Expressed_Genes.html

        python !{assembled_filtering} multiple_added_merged_cut_filtered_not_c11.tabular filtered_c11_gt_0.tabular "c11__gt__=0" 12 "str,str,int,str,str,str,str,int,int,float,float,float" 0
        wc -l filtered_c11_gt_0.tabular | awk '{ print $1 }' >> filtered_c11_gt_0_count.tabular
        python !{correlation} filtered_c11_gt_0.tabular filtered_c11_gt_0_cor.tabular 8,9 pearson 
        ( LC_ALL=C  sort   --stable -t$"\t"  -k "1n,1"  ) < filtered_c11_gt_0.tabular > filtered_c11_gt_0_sorted.tabular

        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Total Upregulated Significant Genes/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header filtered_c11_gt_0_count.tabular concatenate_dataset_3_0.html 1

        tail --lines 1 filtered_c11_gt_0_cor.tabular > filtered_c11_gt_0_cor_tail_1.tabular
        perl !{cut_tool} filtered_c11_gt_0_cor_tail_1.tabular "c1" "T" filtered_c11_gt_0_cor_tail_1_cut.tabular
        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Pearson Correlation/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header filtered_c11_gt_0_cor_tail_1_cut.tabular concatenate_dataset_3_1.html 1

        head --lines 1 filtered_c11_gt_0_sorted.tabular > filtered_c11_gt_0_sorted_head_1.tabular
        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Gene with the Lowest Change/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header filtered_c11_gt_0_sorted_head_1.tabular concatenate_dataset_3_2.html 1

        tail --lines 1 filtered_c11_gt_0_sorted.tabular > filtered_c11_gt_0_sorted_tail_1.tabular
        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Gene with the Highest Change/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header filtered_c11_gt_0_sorted_tail_1.tabular concatenate_dataset_3_3.html 1


        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Genes/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Chromosome Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header filtered_c11_gt_0.tabular concatenate_dataset_3_4.html 1

        cat concatenate_dataset_3_0.html concatenate_dataset_3_1.html concatenate_dataset_3_2.html concatenate_dataset_3_3.html concatenate_dataset_3_4.html > Upregulated_Expressed_Genes.html 
        

        python !{assembled_filtering} multiple_added_merged_cut_filtered_not_c11.tabular filtered_c11_lt_0.tabular "c11__lt__=0" 12 "str,str,int,str,str,str,str,int,int,float,float,float" 0
        wc -l filtered_c11_lt_0.tabular | awk '{ print $1 }' >> filtered_c11_lt_0_count.tabular
        python !{correlation} filtered_c11_lt_0.tabular filtered_c11_lt_0_cor.tabular 8,9 pearson 
        ( LC_ALL=C  sort   --stable -t$"\t"  -k "1n,1"  ) < filtered_c11_lt_0.tabular > filtered_c11_lt_0_sorted.tabular

        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Total Downregulated Significant Genes/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header filtered_c11_lt_0_count.tabular concatenate_dataset_4_0.html 1

        tail --lines 1 filtered_c11_lt_0_cor.tabular > filtered_c11_lt_0_cor_tail.tabular
        perl !{cut_tool} filtered_c11_lt_0_cor_tail.tabular "c1,c2" "T" filtered_c11_lt_0_cor_tail_cut.tabular
        cat !{html_header} | sed -e "s/COL_NUM/2/" -e "s/TABLE_NAME/Pearson Correlation/" -e "s/HEADER_DETAIL//" > moded_header
        perl !{tabular2HTML} moded_header filtered_c11_lt_0_cor_tail_cut.tabular concatenate_dataset_4_1.html 1

        head --lines 1 filtered_c11_lt_0_sorted.tabular > filtered_c11_lt_0_sorted_head_1.tabular
        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Gene with the Lowest Change/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header filtered_c11_lt_0_sorted_head_1.tabular concatenate_dataset_4_2.html 1

        tail --lines 1 filtered_c11_lt_0_sorted.tabular > filtered_c11_lt_0_sorted_tail_1.tabular
        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Gene with the Highest Change/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Location Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header filtered_c11_lt_0_sorted_tail_1.tabular concatenate_dataset_4_3.html 1


        cat !{html_header} | sed -e "s/COL_NUM/12/" -e "s/TABLE_NAME/Significant Genes/" -e "s/HEADER_DETAIL/$(sh gen_header_detail.sh Gene Chromosome Length Explore Gene__Interaction Explore Explore Sample__1__Read__Count Sample__2__Read__Count Mean__Count Fold__Change P-value)/" > moded_header && sed -i "s/__/ /g" moded_header    
        perl !{tabular2HTML} moded_header filtered_c11_lt_0_sorted.tabular concatenate_dataset_4_4.html 1

        cat concatenate_dataset_4_0.html concatenate_dataset_4_1.html concatenate_dataset_4_2.html concatenate_dataset_4_3.html concatenate_dataset_4_4.html > Downregulated_Expressed_Genes.html 


        '''

}



workflow.onComplete {
    println(workflow.success ? "\n\t Final output: ./output/DESeq2/genes.diff\n\t Reports are in '${baseDir}/output/Report/'\n\t Done!\n\n": "\n\t [ERROR] - Oops .. something went wrong")
}


