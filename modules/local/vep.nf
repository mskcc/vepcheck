process VEP {

    tag "${meta.id}_vep_${vep_version}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/vcf2maf:1.6.17':
        'docker.io/mskcc/vcf2maf:1.6.17' }"

    input:
    tuple val(meta),         path(inputVcf),     val(vep_version),  path(vep_path),      path(vep_data)

    output:
    tuple val(output_meta), path("*.vcf")     , emit: vcf
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forks = task.cpus * 2
    def vep_script = "vep"
    output_meta = meta.clone()
    output_meta.vep_version = vep_version
    output_meta.type = "vcf2maf_light"
    if (vep_version.toInteger() < 88){
        vep_script = "variant_effect_predictor.pl"
    }

    """
    perl ${vep_path}/${vep_script} \\
        ${args} \\
        --dir ${vep_data} \\
        --fasta ${vep_data}/homo_sapiens/${vep_version}_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
        --input_file ${inputVcf} \\
        --output_file ${prefix}.vep_${vep_version}.vcf \\
	--fork ${forks} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VEP: ${vep_version}
    END_VERSIONS
    """

}
