process SEQUENZAUTILS_RSEQZ {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-sequenza%3A3.0.0--r42h3342da4_5' :
        'biocontainers/r-sequenza%3A3.0.0--r42h3342da4_5' }"

    input:
    tuple val(meta), path(seqz_bin)

    output:
    tuple val(meta), path("*pdf"), emit: plot
    tuple val(meta), path("*txt"), emit: txt
    tuple val(meta), path("*RData"), emit: data

    when:
    task.ext.when == null || task.ext.when

    script:
    def seqz_in = "${seqz_bin.toString().minus(".gz")}"
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    zcat ${seqz_bin} > $seqz_in
    sequenza.R ${seqz_in} ${meta.id} ${meta.gender}
    """
}
