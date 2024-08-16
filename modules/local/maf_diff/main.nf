process MAF_DIFF {


    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/vep_check_maff_diff:1.0.0':
        'docker.io/mskcc/vep_check_maff_diff:1.0.0' }"

    input:
    tuple val(meta), path(mafs, arity: '2..*'), val(labels)

    output:
    tuple val(meta), path("*.html")  , emit: html_output
    tuple val(meta), path("*.tsv")   , emit: tsv_output
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    maf_diff.py \\
        --mafs ${mafs} \\
        --labels ${labels}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maf_diff: 1.0.0
    END_VERSIONS
    """

}
