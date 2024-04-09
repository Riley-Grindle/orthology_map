
process HMMER_HMMSCAN {
    
    tag '${meta.id}'
    beforeScript = 'ulimit -Ss'
    label 'process_high'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4' :
    'rgrindle/transdecoder' }"

    input:
    path(hmmfile)
    tuple val(meta), path(longest_orfs)

    output:
    tuple val(meta), path("*.domtblout"), emit: table
    path "versions.yml"

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def output = "Pfam.domtblout"
    """
    hmmpress $hmmfile

    hmmscan --cpu $task.cpus --domtblout $output $hmmfile $longest_orfs $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
