
process CLAIR3 {
  tag "$meta.id"
  label 'process_high'
  conda 'bioconda::clair3==1.0.8'
  container 'quay.io/biocontainers/clair3:1.0.8--py39h46983ab_2'

  input:
  tuple val(meta), path(bam), path(index)
  tuple val(meta2), path(ref_fasta), path(ref_fasta_index)
  // optional model_path
  path model_path 

  output:
  tuple val(meta.id), path(vcf), emit: vcf
  path (clair3_dir), emit: output_dir
  path (clair3_log), emit: log
  path "versions.yml" , emit: versions

  script:
  vcf          = "${meta.id}.clair3.vcf.gz"
  clair3_dir   = "${meta.id}.clair3"
  clair3_log   = "${clair3_dir}/run_clair3.log"
  model_suffix = "models/${params.clair3_variant_model}"
  using_conda = (workflow.containerEngine == null || workflow.containerEngine == '')
  """
  CLAIR_BIN_DIR=\$(dirname \$(which run_clair3.sh))
  if [[ "${params.clair3_user_variant_model}" != "" ]] ; then
      MODEL_PATH=${params.clair3_user_variant_model} 
  else
      if [[ ${using_conda} = true ]] ; then
          MODEL_PATH="\$CLAIR_BIN_DIR/${model_suffix}"
      else [[ ${using_conda} = false ]]
          MODEL_PATH="/opt/models/${params.clair3_variant_model}"
      fi
  fi

  run_clair3.sh \\
      --bam_fn=${bam[0]} \\
      --ref_fn=$ref_fasta \\
      --model_path="\$MODEL_PATH"\\
      --threads=${task.cpus} \\
      --platform="ont" \\
      --output=${clair3_dir} \\
      --enable_long_indel \\
      --keep_iupac_bases \\
      --include_all_ctgs \\

  ln -s ${clair3_dir}/merge_output.vcf.gz ${vcf}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      clair3: \$(head -n1 ${clair3_dir}/run_clair3.log | sed 's/^.*CLAIR3 VERSION: v//; s/ .*\$//')
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
