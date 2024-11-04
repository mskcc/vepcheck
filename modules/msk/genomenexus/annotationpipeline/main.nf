process GENOMENEXUS_ANNOTATIONPIPELINE {
    tag "$meta.id"
    errorStrategy 'ignore'
    label 'process_high'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/veptest_genomenexus_annotation-pipeline:1.0.3':
        'docker.io/mskcc/veptest_genomenexus_annotation-pipeline:1.0.3' }"

    input:
    tuple val(meta), path(input_maf)

    output:
    tuple val(output_meta), path("*.maf"), emit: annotated_maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output_meta = meta.clone()
    output_meta.vep_version = "112"
    output_meta.type = "genome_nexus"

    """
    java -Xms${task.memory.toMega()/4}m \\
        -Xmx${task.memory.toGiga()}g \\
        -jar /genome-nexus-annotation-pipeline/annotationPipeline/target/annotationPipeline.jar \\
        -i mskcc --filename ${input_maf} --output-filename \\
        ${input_maf.baseName}.genomenexus_annotation_1.0.3.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomenexus: 'annotation pipeline version 1.0.3'
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${meta.id}_annotated.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomenexus: 'annotation pipeline version 1.0.3'
    END_VERSIONS
    """
}
