process BGZIP_UNZIP {


    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/htslib:1.9':
        'docker.io/mskcc/htslib:1.9' }"

    input:
    tuple val(meta), path(compressed_vcf), path(inputVcfIndex)

    output:
    tuple val(meta), path("*.vcf")  , emit: vcf
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /usr/local/bin/bgzip \\
        -d \\
        ${compressed_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.9
        htslib: 1.9
        tabix: 1.9
    END_VERSIONS
    """

}
