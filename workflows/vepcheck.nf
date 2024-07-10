/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { VCF2MAF }                from '../modules/local/vcf2maf'
include { VEP }                    from '../modules/local/vep'
include { GENOMENEXUS_VCF2MAF as vcf2maf_on_vep; GENOMENEXUS_VCF2MAF as vcf2maf_for_annotation } from '../modules/msk/genomenexus/vcf2maf/main'
include { GENOMENEXUS_ANNOTATIONPIPELINE } from '../modules/msk/genomenexus/annotationpipeline/main'
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

    ch_fasta_ref = Channel.value([ "reference_genome", file(params.fasta) ])
    ref_index_list = []
    for(single_genome_ref in params.fasta_index){
        ref_index_list.add(file(single_genome_ref))
    }
    ch_fasta_fai_ref = Channel.value([ "reference_genome_index",ref_index_list])
    ch_exac_filter = Channel.value(["exac_filter", file(params.exac_filter)])
    ch_exac_filter_index = Channel.value(["exac_filter_index", file(params.exac_filter_index)])

    vep_version = []
    vep_list = []
    vep_cache = []

    for(single_vep_version in params.vep_version){
        vep_version.add(single_vep_version)
    }

    for(single_vep in params.vep){
        vep_list.add(file(single_vep))
    }

    for(single_vep_cache in params.vep_cache){
        vep_cache.add(file(single_vep_cache))
    }

    if (vep_version.size != vep_list.size ||  vep_list.size != vep_cache.size) {
        error "Error: vep version, cache and binary lists are not the same size "
    }

    vep_index = 0
    vep_data_list = []

    while(vep_index < vep_list.size){
        vep_item = [vep_version[vep_index], vep_list[vep_index], vep_cache[vep_index]]
        vep_data_list.add(vep_item)
        vep_index++
    }

    vep_data_channel = Channel.fromList(vep_data_list)


    ch_versions = Channel.empty()


    combined_vcf_and_vep = ch_samplesheet.combine(vep_data_channel)


    VCF2MAF(combined_vcf_and_vep,
            ch_fasta_ref,
            ch_fasta_fai_ref,
            ch_exac_filter,
            ch_exac_filter_index)

    VEP(combined_vcf_and_vep)

    ch_versions = ch_versions.mix(VCF2MAF.out.versions)


    vcf2maf_for_annotation ( ch_samplesheet )

    vcf2maf_on_vep ( VEP.out.vcf )

    ch_versions = ch_versions.mix(vcf2maf_for_annotation.out.versions)

    GENOMENEXUS_ANNOTATIONPIPELINE( vcf2maf_for_annotation.out.maf )

    ch_versions = ch_versions.mix(GENOMENEXUS_ANNOTATIONPIPELINE.out.versions)




    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")

    emit:
    versions       = ch_versions                 // channel: [ file(versions.yml) ]
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
