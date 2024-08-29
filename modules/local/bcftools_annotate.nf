process BCFTOOLS_ANNOTATE {


    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/htslib:1.9':
        'docker.io/mskcc/htslib:1.9' }"

    input:
    tuple val(meta), path(inputVcf), path(vcf_index)

    output:
    tuple val(meta), path("*.vcf")                 , emit: vcf
    path "versions.yml"                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus * 2
    def set_id_cmd =

    """
    /usr/bin/bcftools annotate \\
        --set-id '%CHROM\\_%POS\\_%REF\\_%ALT{0}\\_%ALT{1}' \\
        --threads ${threads} \\
        --output ${prefix}.annotate.vcf \\
        ${inputVcf}



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.9
        htslib: 1.9
    END_VERSIONS
    """

}
