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
include { BCFTOOLS_ANNOTATE } from '../modules/local/bcftools_annotate'
include { BCFTOOLS_NORM } from '../modules/local/bcftools_norm'
include { TABIX } from '../modules/local/tabix'
include { BGZIP_ZIP } from '../modules/local/bgzip_zip'
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

    BCFTOOLS_NORM(ch_samplesheet)

    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    BGZIP_ZIP(BCFTOOLS_NORM.out.vcf)

    ch_versions = ch_versions.mix(BGZIP_ZIP.out.versions)

    TABIX(BGZIP_ZIP.out.gz_vcf)

    ch_versions = ch_versions.mix(TABIX.out.versions)

    vcf_and_index = join_vcf_with_index(BGZIP_ZIP.out.gz_vcf, TABIX.out.vcf_index)

    BCFTOOLS_ANNOTATE(vcf_and_index)

    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)

    combined_vcf_and_vep = BCFTOOLS_ANNOTATE.out.vcf.combine(vep_data_channel)


    VCF2MAF(combined_vcf_and_vep,
            ch_fasta_ref,
            ch_fasta_fai_ref,
            ch_exac_filter,
            ch_exac_filter_index)

    VEP(combined_vcf_and_vep)

    ch_versions = ch_versions.mix(VCF2MAF.out.versions)


    vcf2maf_for_annotation ( BCFTOOLS_ANNOTATE.out.vcf )

    vcf2maf_on_vep ( VEP.out.vcf )

    ch_versions = ch_versions.mix(vcf2maf_for_annotation.out.versions)

    GENOMENEXUS_ANNOTATIONPIPELINE( vcf2maf_for_annotation.out.maf )

    ch_versions = ch_versions.mix(GENOMENEXUS_ANNOTATIONPIPELINE.out.versions)

    maf_diff_input = create_maf_diff_list(VCF2MAF.out.maf, vcf2maf_on_vep.out.maf, GENOMENEXUS_ANNOTATIONPIPELINE.out.annotated_maf)

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
    versions         = ch_versions                 // channel: [ file(versions.yml) ]
    vep_vcf          = VEP.out.vcf
    maf              = VCF2MAF.out.maf
    biallelic_vcf    = BCFTOOLS_NORM.out.vcf
    id_vcf           = BCFTOOLS_ANNOTATE.out.vcf
    annotated_maf    = vcf2maf_for_annotation.out.maf
    vep_maf          = vcf2maf_on_vep.out.maf
    genome_nexus_maf = GENOMENEXUS_ANNOTATIONPIPELINE.out.annotated_maf
    maff_diff_html   = MAF_DIFF.out.html_output
    maff_diff_tsv    = MAF_DIFF.out.tsv_output

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

def create_maf_diff_list(vcf2maf_maf, vcf2maf_light_maf, nexus_maf){
    maf_diff_input = Channel.empty()
    test_vcf2maf = vcf2maf_maf
                    .mix(
                        vcf2maf_light_maf
                    )
                    .mix(
                        nexus_maf
                    )
                    .collect(flat: true)
                    .flatMap{
                        input_list = []
                        vcf2maf_list = []
                        vcf2maf_light_list = []
                        nexus_list = []
                        combo_list = []
                        vep_dict = [:]
                        elem_dict = [:]
                        it_index = 0
                        while(it_index < it.size()){
                            meta = it[it_index]
                            file = it[it_index+1]
                            id = meta.type + "_" + meta.vep_version
                            if(meta.type == 'vcf2maf'){
                                vcf2maf_list.add([id,file])
                            }
                            if(meta.type == 'vcf2maf_light'){
                                vcf2maf_light_list.add([id,file])
                            }
                            if(meta.type == 'genome_nexus'){
                                nexus_list.add([id,file])
                            }
                            elem_dict[id] = file

                            for(single_vep_version in params.vep_version){
                                if(meta.vep_version == single_vep_version){
                                    if(single_vep_version in vep_dict){
                                        vep_dict[single_vep_version].add([id,file])
                                    }
                                    else{
                                        vep_dict[single_vep_version] = [[id,file]]
                                    }

                                }
                            }
                            it_index = it_index + 2
                        }

                        for(single_vep_key in vep_dict.keySet()){
                                events = vep_dict[single_vep_key]
                                events_index = 0
                                label_list = []
                                file_list = []
                                while(events_index < events.size()){
                                    label_list.add(events[events_index][0])
                                    file_list.add(events[events_index][1])
                                    events_index++
                                }
                                input_list.add(new Tuple(['id': label_list.join("_")], file_list, label_list))
                        }

                        all_vcf2maf_index = 0
                        label_list = []
                        file_list = []
                        while(all_vcf2maf_index < vcf2maf_list.size()){
                                events = vcf2maf_list[all_vcf2maf_index]
                                label_list.add(events[0])
                                file_list.add(events[1])
                                all_vcf2maf_index++
                            }
                        input_list.add(new Tuple(['id': label_list.join("_")], file_list, label_list))
                        all_vcf2maf_light_index = 0
                        label_list = []
                        file_list = []
                        while(all_vcf2maf_light_index < vcf2maf_light_list.size()){
                                events = vcf2maf_light_list[all_vcf2maf_light_index]
                                label_list.add(events[0])
                                file_list.add(events[1])
                                all_vcf2maf_light_index++
                            }
                        input_list.add(new Tuple(['id': label_list.join("_")], file_list, label_list))



                        for(single_vep_version in params.vep_version){
                            for(other_vep_version in params.vep_version){
                                if(single_vep_version > other_vep_version){
                                    first_vcf2maf_id = "vcf2maf_" + single_vep_version
                                    second_vcf2maf_id = "vcf2maf_" + other_vep_version
                                    vcf2maf_label_list = [first_vcf2maf_id, second_vcf2maf_id]
                                    vcf2maf_file_list = [elem_dict[first_vcf2maf_id], elem_dict[second_vcf2maf_id]]
                                    first_vcf2maf_light_id = "vcf2maf_light_" + single_vep_version
                                    second_vcf2maf_light_id = "vcf2maf_light_" + other_vep_version
                                    vcf2maf_light_label_list = [first_vcf2maf_light_id, second_vcf2maf_light_id]
                                    vcf2maf_light_file_list = [elem_dict[first_vcf2maf_light_id], elem_dict[second_vcf2maf_light_id]]
                                    input_list.add(new Tuple(['id': vcf2maf_label_list.join("_")], vcf2maf_file_list, vcf2maf_label_list))
                                    input_list.add(new Tuple(['id': vcf2maf_light_label_list.join("_")], vcf2maf_light_file_list, vcf2maf_light_label_list))
                                }
                            }

                        }
                        return input_list




                    }

}


def create_maf_diff_list_old(vcf2maf_maf, vcf2maf_light_maf){
    maf_diff_input = Channel.empty()
    test_all_vcf2maf = vcf2maf_maf
                        .collect(flat: false)
                        .map{
                            it_index = 0
                            maf_label_list = []
                            maf_file_list = []
                            while(it_index < it.size()){
                                label = it[it_index][0].type + "_" + it[it_index][0].vep_version
                                if(!(label in maf_label_list)){
                                    maf_label_list.add(it[it_index][0].type + "_" + it[it_index][0].vep_version)
                                    maf_file_list.add(it[it_index][1])
                                }
                                it_index++
                            }
                            new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)
                        }

    maf_diff_input = maf_diff_input.mix(test_all_vcf2maf)

    test_all_vcf2maf_light = vcf2maf_light_maf
                        .collect(flat: false)
                        .map{
                            it_index2 = 0
                            maf_label_list = []
                            maf_file_list = []
                            while(it_index2 < it.size()){
                                label = it[it_index2][0].type + "_" + it[it_index2][0].vep_version
                                if(!(label in maf_label_list)){
                                    maf_label_list.add(it[it_index2][0].type + "_" + it[it_index2][0].vep_version)
                                    maf_file_list.add(it[it_index2][1])
                                }
                                it_index2++
                            }
                            new Tuple(['id': maf_label_list.join("_")], maf_file_list, maf_label_list)
                        }

    maf_diff_input = maf_diff_input.mix(test_all_vcf2maf_light)

    for(single_vep_version in params.vep_version){
        vep_maf_this_vep = vcf2maf_maf
                            .collect(flat: false)
                            .flatMap()
                            .filter{
                                it[0].vep_version == single_vep_version
                            }
        vep_maf_other_veps = vcf2maf_maf
                            .collect(flat: false)
                            .flatMap()
                            .filter{
                                it[0].vep_version != single_vep_version
                            }
        combined_compare = vep_maf_this_vep.combine(vep_maf_other_veps)
                            .collect(flat: false)
                            .flatMap{
                                compare_list = []
                                it_index3 = 0
                                maf_label_list = []
                                maf_file_list = []
                                while(it_index3 < it.size()){
                                        if(!it){
                                            continue
                                        }
                                        if(!it[it_index3]){
                                            continue
                                        }
                                        println(it)
                                        println(it_index3)
                                        print(it[it_index3])

                                        vep_version1 = it[it_index3][0].vep_version
                                        vep_version2 = it[it_index3][2].vep_version
                                        if (vep_version1 > vep_version2){
                                            label = it[it_index3][0].type + "_" + it[it_index3][0].vep_version + "_" + it[it_index3][2].type + "_" + it[it_index3][2].vep_version
                                        }
                                        else{
                                            label = it[it_index3][2].type + "_" + it[it_index3][2].vep_version + "_" + it[it_index3][0].type + "_" + it[it_index3][0].vep_version
                                        }
                                        if(!(label in maf_label_list)){
                                            compare_list.add(new Tuple(['id': label], [it[it_index3][1], it[it_index3][3]], [it[it_index3][0].type + "_" + it[it_index3][0].vep_version, it[it_index3][2].type + "_" + it[it_index3][2].vep_version]))
                                        }
                                    it_index3++

                                }
                                return compare_list
                            }
        maf_diff_input = maf_diff_input.mix(combined_compare)

    }

    for(single_vep_version in params.vep_version){
        vep_maf_this_vep = vcf2maf_light_maf
                            .collect(flat: false)
                            .flatMap()
                            .filter{
                                it[0].vep_version == single_vep_version
                            }
        vep_maf_other_veps = vcf2maf_light_maf
                            .collect(flat: false)
                            .flatMap()
                            .filter{
                                it[0].vep_version != single_vep_version
                            }
        combined_compare = vep_maf_this_vep.combine(vep_maf_other_veps)
                            .collect(flat: false)
                            .flatMap{
                                compare_list = []
                                it_index4 = 0
                                maf_label_list = []
                                maf_file_list = []
                                while(it_index4 < it.size()){
                                    if(it[it_index4])
                                    {
                                        vep_version1 = it[it_index4][0].vep_version
                                        vep_version2 = it[it_index4][2].vep_version
                                        if (vep_version1 > vep_version2){
                                            label = it[it_index4][0].type + "_" + it[it_index4][0].vep_version + "_" + it[it_index4][2].type + "_" + it[it_index4][2].vep_version
                                        }
                                        else{
                                            label = it[it_index4][2].type + "_" + it[it_index4][2].vep_version + "_" + it[it_index4][0].type + "_" + it[it_index4][0].vep_version
                                        }
                                        if(!(label in maf_label_list)){
                                            compare_list.add(new Tuple(['id': label], [it[it_index4][1], it[it_index4][3]], [it[it_index4][0].type + "_" + it[it_index4][0].vep_version, it[it_index4][2].type + "_" + it[it_index4][2].vep_version]))
                                        }
                                    }

                                    it_index4++
                                }

                                return compare_list
                            }
        maf_diff_input = maf_diff_input.mix(combined_compare)

    }
    vep_id_vcf2maf_maf = vcf2maf_maf
                        .collect(flat: false)
                        .flatMap()
                        .map{
                            new Tuple(it[0].vep_version, it[0], it[1])
                        }
    vep_id_vcf2maf_light_maf = vcf2maf_light_maf
                        .collect(flat: false)
                        .flatMap()
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
            maf_file_list.add(it[1])

        }
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
