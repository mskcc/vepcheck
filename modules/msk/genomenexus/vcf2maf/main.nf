process GENOMENEXUS_VCF2MAF {
    tag "$meta.id"
    label 'process_high'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/veptest_vcf2maf_lite:latest':
        'docker.io/mskcc/veptest_vcf2maf_lite:latest' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 /vcf2maf-lite/vcf2maf_lite.py -i ${vcf} ${args}
    mv vcf2maf_output/${vcf.baseName}.maf ${vcf.baseName}.vcf2maf_lite_0.0.1.maf
    echo '"${task.process}": vcf2maf_lite.py v.0.0.1' > versions.yml
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    mkdir vcf2maf_output
    touch ${meta.id}.maf
    cp ${meta.id}.maf vcf2maf_output/
    echo '"${task.process}": vcf2maf_lite.py v.0.0.1' > versions.yml
    """
}
