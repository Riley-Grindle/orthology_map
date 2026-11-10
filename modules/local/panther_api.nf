process PANTHER_API {
    tag "Pulling gene ids from PantherDB"
    label "process_medium"

    // Hard cap on concurrent Float/spot submissions for this process,
    // independent of how many chunks a given run's splitText(by: ...)
    // produces (see workflows/orthologymap.nf) -- keeps this stage to at
    // most 100 machines running at once regardless of species/input size.
    maxForks 100

    conda "conda-forge::python=3.9.5"
    // Switched off the bespoke rgrindle/panther_api_sf image: this process no
    // longer runs curl or anything else that image was built for (see
    // bin/post_panther.sh -- gene lookup is now a local file read, not a live
    // API call), and that image has started resolving (via Float's registry
    // mirror, quay.io/rgrindle/panther_api_sf) to a pull with no python3/python
    // on PATH at all -- confirmed via direct `ls /usr/bin/python* /opt/conda/bin/python*
    // /usr/local/bin/python*` inside a real failing task, all empty, not a
    // transient spot-interruption artifact.
    //
    // Plain `python:3.9-slim` doesn't work either: Float rewrites every bare
    // image reference to quay.io/<name> before submitting (confirmed via a
    // real submission log: `--image quay.io/python:3.9-slim` -> "cannot find
    // image or image is inaccessible quay.io/python"), and quay.io has no
    // mirror of Docker Hub's official `python` images under that bare name.
    // `biocontainers/python:3.9--1` doesn't have this problem because
    // biocontainers actually dual-publishes under the same name to both
    // Docker Hub and quay.io -- same tag already used successfully elsewhere
    // in this pipeline (GUNZIP, PREP_INPUT, GTF_TO_GENELISTS, COLLECT_CHUNKS,
    // STAGE_OUTS all use it and have run to completion via Float already).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path tree_out
    val ref_panther_prefixes

    output:
    path "*orthologs.txt"

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    post_panther.sh "${ref_panther_prefixes}"
    mv orthologs.txt "\$(basename $tree_out)_orthologs.txt"
    # Harmless cache-busting comment, v4: fixed a subfamily-id extraction
    # bug in post_panther.sh (grep -o re-anchoring on a literal "SF"
    # substring embedded in a PAINT-curated protein NAME, e.g. "NSFL1",
    # silently corrupting the lookup key and producing a false "NO MATCH"
    # even when a real reference-species ortholog exists in that exact
    # subfamily). splitText(by: ...) also changed in workflows/orthologymap.nf,
    # which alone forces every chunk to re-execute under -resume -- this
    # comment is just belt-and-suspenders documentation of why.
    """

}
