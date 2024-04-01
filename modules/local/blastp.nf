
process BLASTP {

    label 'process_medium'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4' :
    'rgrindle/transdecoder' }"

    input:
    path(query_fasta)
    path(reference_fasta)    

    output:
    path "all_v_all.tsv"             , emit: tbl
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    makeblastdb -in $reference_fasta \\
        -dbtype prot \\
        -out reference_db \\
        ${args}
    
    blastp -db reference_db \\
           -query $query_fasta \\
           -outfmt 6 \\
           -out all_v_all.tsv \\
           -num_threads ${task.cpus} \\
           ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastp: \$(echo \$(blastp --version) | sed -e "s/blastp //g")
    END_VERSIONS
    """
}
