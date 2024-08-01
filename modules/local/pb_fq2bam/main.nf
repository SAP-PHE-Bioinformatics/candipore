process PB_FQ2BAM {

    label 'process_gpu_high'
    tag "SAMPLE: ${meta.id}" 
    publishDir "${params.outdir}/bams/${sample}", mode: 'symlink'
    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2),  path(fasta)
	tuple val(meta3), path(index)

    output:
    tuple val(sample), path("*.bam"), path("*.bai"), emit: bam
    path("*_qc_metrics"), emit: qc_metrics
    path("*_pbrun_fq2bam_log.txt"), emit: metrics_logs
    path("*_duplicate_metrics.txt"), emit: duplicate_metrics
    path "versions.yml"                           , emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    

    """
    cp -L bwa/${fasta.baseName}.* .

        pbrun fq2bam \\
          --ref ${fasta} \\
          $in_fq_command \\
          --read-group-sm ${prefix} \\
          --out-bam ${prefix}_markdup.bam \\
          --out-duplicate-metrics ${prefix}_duplicate_metrics.txt \\
          --out-qc-metrics-dir ${prefix}_qc_metrics \\
          --bwa-options="-M" \\
          --fix-mate \\
          --optical-duplicate-pixel-distance 2500 \\
          --logfile ${prefix}_pbrun_fq2bam_log.txt \\
          $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}