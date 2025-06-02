process SEQUENZAUTILS_MERGESEQZ {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_retry'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenza-utils:3.0.0--py39h67e14b5_5' :
        'quay.io/biocontainers/sequenza-utils:3.0.0--py39h67e14b5_5' }"
    
    input:
    tuple val(meta), path(seqz)

    output:
    tuple val(meta), path("*concat.seqz.gz"), path("*concat.seqz.gz.tbi"), emit: concat_seqz
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    head -n 1  <( zcat ${seqz[0]} ) > header
	cat header <( zcat $seqz | awk '{if (NR!=1 && \$1 != "chromosome") {print \$0}}' ) | bgzip > \
    ${prefix}_concat.seqz.gz
    tabix -f -s 1 -b 2 -e 2 -S 1 ${prefix}_concat.seqz.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: 3.0.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "
    zcat $seqz | \
    awk '{if (NR!=1 && dollar1 != "chromosome") {print dollar0}}' | bgzip > \
    ${prefix}_concat.seqz.gz
    tabix -f -s 1 -b 2 -e 2 -S 1 ${prefix}_concat.seqz.gz
    "
    touch ${prefix}_concat.seqz.gz
    touch ${prefix}_concat.seqz.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: stub
    END_VERSIONS
    """
}