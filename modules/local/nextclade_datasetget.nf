process NEXTCLADE_DATASETGET {
    tag "$meta.id"
    label 'process_low'

    container = null
    executor = 'local'

    input:
    tuple val(meta), path(dataset)

    output:
    tuple val(meta), path("$prefix") , emit: dataset_2
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dataset}_2"

    """
    nextclade \\
        dataset \\
        get \\
        $args \\
        --name $dataset \\
        --output-dir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
    END_VERSIONS
    """
}
