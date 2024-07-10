/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { VCF2MAF }                from '../../../modules/local/vcf2maf'
include { TABIX }                  from '../../../modules/local/tabix'
include { GENOMENEXUS_VCF2MAF as vcf2maf_on_vep; GENOMENEXUS_VCF2MAF as vcf2maf_for_annotation} from '../../../modules/msk/genomenexus/vcf2maf/main'
include { GENOMENEXUS_ANNOTATIONPIPELINE     } from '../../../modules/msk/genomenexus/annotationpipeline/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_vepcheck_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VEPCHECK {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_fasta_ref = Channel.value([ "reference_genome", path(params.fasta) ])
    ref_index_list = []
    for(single_genome_ref in params.fasta_index){
        ref_index_list.add(path(single_genome_ref))
    }
    ch_fasta_fai_ref = Channel.value([ "reference_genome_index",ref_index_list])
    ch_exac_filter = Channel.value(["exac_filter", path(params.exac_filter)])
    ch_exac_filter_index = Channel.value(["exac_filter_index", path(params.exac_filter_index)])

    vep_version = []
    vep_list = []
    vep_cache = []

    for(single_vep_version in params.vep_version){
        vep_version.add(params.single_vep_version)
    }

    for(single_vep in params.vep){
        vep_list.add(path(single_vep))
    }

    for(single_vep_cache in params.vep_cache){
        vep_cache.add(path(single_vep_cache))
    }

    if (vep_version.size != vep_list.size ||  vep_list.size != vep_cache.size) {
        error "Error: vep version, cache and binary lists are not the same size "
    }

    vep_index = 0
    vep_data = []

    while(vep_index < vep_list.size){
        vep_data.append([vep_version[vep_index], vep_list[vep_index]], vep_cache[vep_index])
    }

    vep_data_channel = Channel.value(vep_data)


    ch_versions = Channel.empty()

    TABIX(ch_samplesheet)

    ch_versions = ch_versions.mix(TABIX.out.versions)

    vcf_and_index = join_vcf_with_index(ch_samplesheet,tabix.out.vcf_index)

    combined_vcf_and_vep = vcf_and_index.combine(vep_data_channel)

    VCF2MAF(combined_vcf_and_vep,
            ch_fasta_ref,
            ch_fasta_fai_ref,
            ch_exac_filter,
            ch_exac_filter_index)

    VEP(combined_vcf_and_vep)

    ch_versions = ch_versions.mix(VCF2MAF.out.versions)

    vcf2maf_for_annotation ( vcf_and_index )

    vcf2maf_on_vep ( VEP.out.vcf )

    ch_versions = ch_versions.mix(GENOMENEXUS_VCF2MAF.out.versions)

    GENOMENEXUS_ANNOTATIONPIPELINE( GENOMENEXUS_VCF2MAF.out.maf )

    ch_versions = ch_versions.mix(GENOMENEXUS_ANNOTATIONPIPELINE.out.versions)




    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectpath(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    vep_vcf        = VEP.out.vcf
    maf            = VCF2MAF.out.maf
    nexus_annotated_maf = vcf2maf_for_annotation.out.maf
    vep_maf        = vcf2maf_on_vep.out.maf
}

def join_vcf_with_index(vcf,index) {
        vcf_channel = vcf
            .map{
                new Tuple(it[0].id,it)
                }
        index_channel = index
            .map{
                new Tuple(it[0].id,it)
                }
        mergedWithKey = vcf_channel
            .join(index_channel)
        merged = mergedWithKey
            .map{
                new Tuple(it[1][0],it[1][1],it[2][1])
            }
        return merged

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
