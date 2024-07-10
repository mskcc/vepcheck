process VCF2MAF {

    tag "${meta.id}_vep_${vep_version}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/vcf2maf:1.6.17':
        'docker.io/mskcc/vcf2maf:1.6.17' }"

    input:
    tuple val(meta),         path(inputVcf),      path(inputVcfIndex), val(vep_version),  path(vep_path),      path(vep_data)
    tuple val(meta2),        path(fasta)
    tuple val(meta3),        path(fai)
    tuple val(meta4),        path(exac_filter)
    tuple val(meta5),        path(exac_filter_tbi)

    output:
    tuple val(meta), path("*.maf")     , emit: maf
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forks = task.cpus * 2

    """
    perl /usr/bin/vcf2maf/vcf2maf.pl \\
        ${args} \\
        --cache-version ${vep_version} \\
        --vep-data ${vep_data} \\
        --vep-path ${vep_path} \\
        --input-vcf ${inputVcf} \\
        --ref-fasta ${fasta} \\
        --filter-vcf ${exac_filter} \\
        --vep-forks ${forks} \\
        --output-maf ${prefix}.vcf2maf_1.6.17.vep_${vep_version}.maf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | head -n 2 | grep -o '(v.*)' | sed 's/[()]//g')
        vcf2maf: 1.6.17
        VEP: 86
        htslib: 1.9
        samtools: 1.9
        bcftools: 1.9
    END_VERSIONS
    """

}
