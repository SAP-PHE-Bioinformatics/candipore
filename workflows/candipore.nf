/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_candipore_pipeline'
include { NANOPLOT               } from '../modules/nf-core/nanoplot/main'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { BEDTOOLS_GENOMECOV     } from '../modules/nf-core/bedtools/genomecov/main'
include { CHOPPER                } from '../modules/nf-core/chopper/main'
include { CLAIR3                 } from '../modules/local/clair3/main'
include { PEPPER_MARGIN_DEEPVARIANT } from '../modules/local/pepper_margin_deepvariant/main'
include { TABIX_BGZIPTABIX      } from '../modules/nf-core/tabix/bgziptabix/main'
include { SNIFFLES              } from '../modules/nf-core/sniffles/main'
include { SAMBAMBA_MARKDUP      } from '../modules/nf-core/sambamba/markdup/main'
include { SAMTOOLS_INDEX        } from '../modules/nf-core/samtools/index/main'
include { PARABRICKS_MINIMAP2   } from '../modules/local/parabricks/minimap2/main'
include { PARABRICKS_FQ2BAM     } from '../modules/nf-core/parabricks/fq2bam/main'
include { BWA_INDEX             } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM               } from '../modules/nf-core/bwa/mem/main'
include { PB_FQ2BAM             } from '../modules/local/pb_fq2bam/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CANDIPORE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_input = Channel.empty()
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_snp_vcf = Channel.empty()

    //housekeeping steps 3
    // get reference name from params.reference filename 
    Channel
        .fromPath( params.reference )
        .map{ file -> [[file.baseName] , file ]}
        .set { ch_reference }
    ch_reference.view()
    //Index Reference
    SAMTOOLS_FAIDX (
        ch_reference,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_reference = ch_reference.join(SAMTOOLS_FAIDX.out.fai)


    //Inital Read QC
    NANOPLOT (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt)

    //OPTION TO FILTER READS
    if (params.filter) {

        CHOPPER (
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(CHOPPER.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(CHOPPER.out.multiqc_files)

        ch_input = CHOPPER.out.fastq
    }
    else {
        ch_input = ch_samplesheet
    }

    ch_bam = Channel.empty()

    // Align reads to reference 
    if (params.cuda) {
        ch_single = ch_input.filter{meta, file -> meta.single_end == true}
        ch_paired = ch_input.filter{meta, file -> meta.single_end == false}

        // PARABRICKS_MINIMAP2(
        //         ch_single,
        //         ch_reference.first()
        //     )
        //     ch_versions = ch_versions.mix(PARABRICKS_MINIMAP2.out.versions)
        //     ch_bam = ch_bam.mix(PARABRICKS_MINIMAP2.out.bam)
      
        BWA_INDEX(
            ch_reference.map{meta, ref, index -> [meta, ref]}
        )
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

        BWA_MEM(
            ch_paired,
            BWA_INDEX.out.index.first(),
            ch_reference.map{meta, ref, index ->[meta, ref]}.first(),
            true
        )
        ch_versions = ch_versions.mix(BWA_MEM.out.versions)

        PARABRICKS_FQ2BAM(
            ch_paired,
            ch_reference.map{meta, ref, index -> [meta, ref]}.first(),
            BWA_INDEX.out.index.first()
        )
        ch_versions = ch_versions.mix(PARABRICKS_FQ2BAM.out.versions)
        ch_bam = ch_bam.mix(PARABRICKS_FQ2BAM.out.bam)

        // PB_FQ2BAM(
        //     ch_paired,
        //     ch_reference.map{meta, ref, index -> [meta, ref]}.first(),
        //     BWA_INDEX.out.index.first()
        // )
        // ch_versions = ch_versions.mix(PB_FQ2BAM.out.versions)
        // ch_bam = ch_bam.mix(PB_FQ2BAM.out.bam)
    }
    // else{ 
    MINIMAP2_ALIGN(
        ch_input,
        ch_reference.first(),
        true,
        'bai',
        false,
        true
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_bam = MINIMAP2_ALIGN.out.bam
    // }

    // Mark Duplicates
    SAMBAMBA_MARKDUP(
        ch_bam,
    )
    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)

    SAMTOOLS_INDEX(
        SAMBAMBA_MARKDUP.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Coverage
    BEDTOOLS_GENOMECOV (
        SAMBAMBA_MARKDUP.out.bam,
        [],
        "bedGraph",
        true,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
    // if (!params.cuda) {

    // Variant Calling using Clair3
    CLAIR3 (
        SAMBAMBA_MARKDUP.out.bam.join(SAMTOOLS_INDEX.out.bai),
        ch_reference.first(),
        params.clair3_model
 )
    ch_versions = ch_versions.mix(CLAIR3.out.versions)
    // compress and index vcf
    TABIX_BGZIPTABIX(
        CLAIR3.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    // }
    // else {
    // Variant Calling using Pepper Margin DeepVariant
        PEPPER_MARGIN_DEEPVARIANT(
            (SAMBAMBA_MARKDUP.out.bam.join(SAMTOOLS_INDEX.out.bai)).filter{meta, bam,bai -> meta.single_end == true},
            ch_reference.first()
        )
        ch_versions = ch_versions.mix(PEPPER_MARGIN_DEEPVARIANT.out.versions)
    
    ch_snp_vcf = PEPPER_MARGIN_DEEPVARIANT.out
    // }

    //Structural Variant Calling
    if (params.svcaller) {
        SNIFFLES(
            SAMBAMBA_MARKDUP.out.bam.join(SAMTOOLS_INDEX.out.bai),
            ch_reference.first()
        )
    }


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
