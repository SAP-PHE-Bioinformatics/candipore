process PEPPER_MARGIN_DEEPVARIANT {
    tag "$meta.id"

    if (params.deepvariant_gpu) {
        container 'docker://kishwars/pepper_deepvariant:r0.8-gpu'
        label 'process_gpu_high', 'process_high'
    } else {
        container 'docker://kishwars/pepper_deepvariant:r0.8'
        label 'process_high'
    }

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*vcf.gz")    ,  emit: vcf
    tuple val(meta), path("*vcf.gz.tbi"),  emit: tbi
    path "versions.yml"                 ,  emit: versions


    script:
    def args    = task.ext.args ?: ''
    def gpu     = params.deepvariant_gpu ? "-g --device_ids 0,2" : ""
    prefix      = task.ext.prefix ?: "${meta.id}"
    //def regions = intervals ? "--regions ${intervals}" : ""
    //def gvcf    = params.make_gvcf ? "--gvcf" : ""

    """
    mkdir -p "${prefix}"

    run_pepper_margin_deepvariant call_variant \\
        -b "${input}" \\
        -f "${fasta}" \\
        -o "." \\
        -p "${prefix}" \\
        -t ${task.cpus} \\
        $args \\
        --gvcf \\
        $gpu \\
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepper_margin_deepvariant: \$(run_pepper_margin_deepvariant --version | sed 's/VERSION: //g')
    END_VERSIONS
    """
}