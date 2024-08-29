process TABIX {


    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/htslib:1.9':
        'docker.io/mskcc/htslib:1.9' }"

    input:
    tuple val(meta), path(inputVcf)

    output:
    tuple val(meta), path("*.tbi")  , emit: vcf_index
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /usr/local/bin/tabix \\
        ${args} \\
        --preset vcf \\
        ${inputVcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.9
        htslib: 1.9
        tabix: 1.9
    END_VERSIONS
    """

}
