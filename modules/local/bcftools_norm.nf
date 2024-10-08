process BCFTOOLS_NORM {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/htslib:1.9':
        'docker.io/mskcc/htslib:1.9' }"

    input:
    tuple val(meta), path(inputVcf)

    output:
    tuple val(meta), path("*.vcf")     , emit: vcf
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus * 2

    """
    /usr/bin/bcftools norm \\
        -N \\
        -m-any \\
        ${inputVcf} \\
        --output ${prefix}.biallelic_split.vcf \\
        --threads ${threads} \\
        ${inputVcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.9
        htslib: 1.9
    END_VERSIONS
    """

}
