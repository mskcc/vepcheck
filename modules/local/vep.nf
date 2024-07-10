process VEP {

    tag "${meta.id}_vep_${vep_version}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/vcf2maf:1.6.17':
        'docker.io/mskcc/vcf2maf:1.6.17' }"

    input:
    tuple val(meta),         path(inputVcf),      path(inputVcfIndex),   val(vep_version),  path(vep_path),      path(vep_data)

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forks = task.cpus * 2

    """
    perl ${vep_path} \\
        ${args} \\
        --dir ${vep_data} \\
        --fasta ${vep_data}/homo_sapiens/${vep_version}_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \\
        --input_file ${inputVcf} \\
        --output_file ${prefix}.vep_${vep_version}.vcf \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VEP: ${vep_version}
    END_VERSIONS
    """

}
