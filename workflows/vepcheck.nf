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
include { MAF_DIFF } from '../modules/local/maf_diff/main'
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

    maf_diff_input = create_maf_diff_list(VCF2MAF.out.maf, vcf2maf_on_vep.out.maf)

    MAF_DIFF( maf_diff_input )

    ch_versions = ch_versions.mix(MAF_DIFF.out.versions)

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
    maff_diff_html = MAF_DIFF.out.html_output
    maff_diff_tsv  = MAF_DIFF.out.tsv_output

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

def create_maf_diff_list(vcf2maf_maf, vcf2maf_light_maf){
    maf_diff_input = Channel.empty()
    test_all_vcf2maf = vcf2maf_maf
                        .collect(flat: false)
                        .flatMap()
                        .map{
                            index = 0
                            maf_label_list = []
                            maf_file_list = []
                            while((index+1) < it.size()){
                                maf_label_list.add(it[index].type + "_" + it[index].vep_version)
                                maf_file_list.add(it[index+1])
                                index = index + 2
                            }
                            new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)
                        }
    test_all_vcf2maf.view()

    maf_diff_input = maf_diff_input.mix(test_all_vcf2maf)

    test_all_vcf2maf_light = vcf2maf_light_maf
                        .collect(flat: false)
                        .flatMap()
                        .map{
                            index = 0
                            maf_label_list = []
                            maf_file_list = []
                            while((index+1) < it.size()){
                                maf_label_list.add(it[index].type + "_" + it[index].vep_version)
                                maf_file_list.add(it[index+1])
                                index = index + 2
                            }
                            new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)
                        }

    maf_diff_input = maf_diff_input.mix(test_all_vcf2maf_light)

    for(single_vep_version in params.vep_version){
        vep_maf_this_vep = vcf2maf_maf
                            .filter{
                                it[0].vep_version == single_vep_version
                            }
        vep_maf_other_veps = vcf2maf_maf
                            .filter{
                                it[0].vep_version != single_vep_version
                            }
        combined_compare = vep_maf_this_vep.combine(vep_maf_other_veps)
                            .map{
                                maf_label_list = []
                                maf_file_list = []
                                maf_label_list.add(it[0].type + "_" + it[0].vep_version)
                                maf_label_list.add(it[2].type + "_" + it[2].vep_version)
                                maf_file_list.add(it[1])
                                maf_file_list.add(it[3])
                                new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)
                            }
        maf_diff_input = maf_diff_input.mix(combined_compare)

    }

    for(single_vep_version in params.vep_version){
        vep_maf_this_vep = vcf2maf_light_maf
                            .filter{
                                it[0].vep_version == single_vep_version
                            }
        vep_maf_other_veps = vcf2maf_light_maf
                            .filter{
                                it[0].vep_version != single_vep_version
                            }
        combined_compare = vep_maf_this_vep.combine(vep_maf_other_veps)
                            .map{
                                maf_label_list = []
                                maf_file_list = []
                                maf_label_list.add(it[0].type + "_" + it[0].vep_version)
                                maf_label_list.add(it[2].type + "_" + it[2].vep_version)
                                maf_file_list.add(it[1])
                                maf_file_list.add(it[3])
                                new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)
                            }

        maf_diff_input = maf_diff_input.mix(combined_compare)

    }

    vep_id_vcf2maf_maf = vcf2maf_maf
                        .map{
                            new Tuple(it[0].vep_version, it[0], it[1])
                        }
    vep_id_vcf2maf_light_maf = vcf2maf_light_maf
                    .map{
                        new Tuple(it[0].vep_version, it[0], it[1])
                    }
    test_comparisons_vcf2maf = vep_id_vcf2maf_maf.combine(vep_id_vcf2maf_light_maf, by: 0)
                                .map{
                                    maf_label_list = []
                                    maf_file_list = []
                                    maf_label_list.add(it[1].type + "_" + it[1].vep_version)
                                    maf_label_list.add(it[3].type + "_" + it[3].vep_version)
                                    maf_file_list.add(it[2])
                                    maf_file_list.add(it[4])
                                    new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)

                                }
    maf_diff_input = maf_diff_input.mix(test_comparisons_vcf2maf)
    return maf_diff_input
}

def compare_vcf2maf_and_light(maf1, maf2) {
    maf_check_input = Channel.empty()
    maf_file_list1 = []
    maf_label_list1 = []
    maf_file_list2 = []
    maf_label_list2 = []
    maf_file1 = maf1
        .map{
            maf_file_list1.add(it[1])
        }
    maf_label1 = maf1
        .map{
            maf_label_list1.add(it[0].type + "_" + it[0].vep_version)
        }
    maf_file2 = maf2
        .map{
            maf_file_list2.add(it[1])
        }
    maf_label2 = maf2
        .map{
            maf_label_list2.add(it[0].type + "_" + it[0].vep_version)
        }
    index = 0

    while ( index < maf_file_list1.size ){
        label_list = [maf_label_list1[index], maf_label_list2[check_index]]
        file_list = [maf_file_list1[index], maf_file_list2[check_index]]
        maf_check_input.mix(
            new Tuple(['id': label_list.join("_")], file_list, label_list)
        )
        index++
    }

    return maf_check_input

}

def create_maf_check_list(maf) {
    maf_file_list = []
    maf_label_list = []
    maf_check_input = Channel.empty()
    maf_file = maf.collect()
        .map{
            print(it)
            maf_file_list.add(it[1])

        }
    print(maf_file_list)
    maf_label = maf.collect()
        .map{
            maf_label_list.add(it[0].type + "_" + it[0].vep_version)

        }

    all_mafs = Channel.of(
        new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)
    )
    maf_check_input.mix(all_mafs)
    index = 0

    while ( index < maf_file_list.size ){
        check_index  = index + 1
        while(check_index < maf_file_list.size){
            label_list = [maf_label_list[index], maf_label_list[check_index]]
            file_list = [maf_file_list[index], maf_file_list[check_index]]
            maf_check_input.mix(
                new Tuple(['id': label_list.join("_")], file_list, label_list)
            )
            check_index++
        }
        index++
    }

    return maf_check_input

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
