#!/usr/bin/env nextflow

// start pipeline as: 
// nextflow run scripts/batchUPdetection.nf \
// --family_files families.csv \
// --vcf_folder vcfs/ \
// --outdir_roh roh/ \
// --outdir_isec isec/ \
// --outdir_results tagging_and_visualisation/ \

project_dir = projectDir

// workflow
workflow {
    ch_vcf_files = Channel.fromPath(file(params.family_files)).splitCsv(header:true).map { row-> tuple(row.index_file, row.mother_file, row.father_file) } 
    roh_for_index(ch_vcf_files)
    isec(ch_vcf_files)
    tag_and_visualize()
}

// roh detection and output
process roh_for_index {
    debug true
    //errorStrategy 'ignore'
    publishDir "$params.outdir_roh"

    input:
        tuple val(index_vcf), val(mother_vcf), val(father_vcf)

    output:
        path "${index_vcf}_roh.txt"
        tuple val(index_vcf), val(mother_vcf), val(father_vcf)

    script:
        """
        bcftools roh --AF-dflt 0.4 -G 30 -I ${params.vcf_folder}${index_vcf} | awk '\$1=="RG"{print \$0}' > ${index_vcf}_roh.txt
        """
}

process isec {

    publishDir "$params.outdir_isec"

    input:
        tuple val(index_vcf), val(mother_vcf), val(father_vcf)

    output:
        path "${index_vcf}_isec.txt"

    script:
        if( mother_vcf != "" & father_vcf != "" )
            """
            bcftools isec -n +1 ${params.vcf_folder}${index_vcf} ${params.vcf_folder}${mother_vcf} ${params.vcf_folder}${father_vcf} > ${index_vcf}_isec.txt
            """
        else if( mother_vcf != "" & father_vcf == "" )
            """
            bcftools isec -n +1 ${params.vcf_folder}${index_vcf} ${params.vcf_folder}${mother_vcf} > ${index_vcf}_isec.txt
            """
        else if( mother_vcf == "" & father_vcf != "" )
            """
            bcftools isec -n +1 ${params.vcf_folder}${index_vcf} ${params.vcf_folder}${father_vcf} > ${index_vcf}_isec.txt
            """
        else if( mother_vcf == "" & father_vcf == "" )
            """
            echo ${index_vcf} - single > ${index_vcf}_isec.txt
            """
}

// merge results of isec and roh
process merge_isec_roh {

    input:
        file(roh_file)
        file(isec_file)
    
    output:
        stdout

    script:
        """
        ### add python script here, that takes roh and isec files, counts rohs and calculates ratios and outputs a single results file and plots
        """
}

// assign upd/consanguinity tags to samples and output results and plots
process tag_and_visualize {

    script:
        """
        python ${project_dir}/upd_finder.py \
            ${params.outdir_roh} \
            ${params.outdir_isec} \
            ${params.family_files} \
            ${params.vcf_folder} \
            ${params.outdir_results}
        """
}


