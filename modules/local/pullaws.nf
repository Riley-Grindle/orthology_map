
process PULL_AWS {
    tag "Pulling AWS reference fastas"
    label "process_low"

    input:
    val refs
    path fasta

    output:
    path fasta

    script:
    """
    aws s3 sync $refs $fasta
    echo $fasta
    """

    stub:
    """
    touch $fasta
    """


}