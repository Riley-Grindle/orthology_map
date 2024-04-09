
process STAGE_OUTS {
    tag "Staging orthology tool outputs..."
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path (ortho_f)
    tuple val(meta), path (ortho_l)
    path egg 
    path tree

    output:
    tuple val(meta), path ("./outs"), emit: all_outs

    script: 
    """
    mkdir outs
    cd outs
    mkdir ortho_l
    mkdir ortho_f
    mkdir egg
    mkdir tree
    cd ../
    cp -r $ortho_f ./outs/ortho_f/; cp -r $ortho_l ./outs/ortho_l/; cp -r $egg ./outs/egg/; cp $tree ./outs/tree/
    """

}
