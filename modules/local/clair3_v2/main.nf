process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:0.1.10--hdfd78af_0' :
        'quay.io/biocontainers/clair3:0.1.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta), path(fai)
    path(model)

    output:
    tuple val(meta), path(vcf), emit: vcf
    path (clair3_dir), emit: output_dir
    path (clair3_log), emit: log
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    clair3_dir   = "${meta.id}.clair3"
    clair3_log   = "${meta.id}.clair3.log"

    """
    /usr/local/bin/run_clair3.sh \
    --bam_fn=$input \
    --ref_fn=$fasta \
    --threads=$task.cpus \
    --platform="ont" \
    --model_path="${model}" \
    --output=${clair3_dir} \
    $args

    ln -s ${clair3_dir}/merge_output.vcf.gz ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$( /usr/local/bin/run_clair3.sh --version | sed 's/ /,/' )
    END_VERSIONS
    """
}