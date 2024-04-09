
process POST_PROC {
    tag "Processing orthology tool outputs..."
    //label "process_medium"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/post_proc' }"

    input:
    tuple val(meta), path (outs)
    tuple val(meta), path (fastas)
    val fasta_file
    path ensembl_data
    path taxa_db

    output:
    path ("./std_outs/")

    script:
    """
    FILE="\$(basename $fasta_file)"
    PREFIX="\$(echo "\$FILE" | rev | cut -d. -f2- | rev)"
    
    if ! [[ -e "./outs/ortho_f/blank.txt" ]]; then
        Rscript /rscripts/ortho_f_post.R ./outs/ortho_f/*/Orthologues/Orthologues_\${PREFIX}
    fi
        
    ORTHO_L_PREFIX="\$(echo "\$FILE" | cut -d. -f1)"
    Rscript /rscripts/ortho_l_post.R ./outs/ortho_l/ "\${ORTHO_L_PREFIX}"
    Rscript /rscripts/eggnog_post.R ./outs/egg
    Rscript /rscripts/tree_post.R ./outs/tree
    mkdir /headers
    for file in ./odbdata/*; do
        FILE_NAME="\$(basename "\$file")"; filename_no_extension="\${FILE_NAME%.*}_headers.txt"
        grep ">" "\$file" | grep -o '^[^[:space:]]*' > "/headers/"\$filename_no_extension""
    done
    Rscript /rscripts/hexadecimal_correction.R ./odbdata /headers/ ./
    mv ./tree_std.csv input_tree.csv
    ensembl_id_2_gene_symbl.py ./input_tree.csv /sup_data/prefix_2_file.json $ensembl_data /sup_data/prefix_2_species.json
    paste -d"," input_tree.csv gene_symbols_ensembl.txt species_names.txt > tree_std.csv
    search_taxa.sh tree_std.csv $taxa_db
    echo "SPECIES" >> tmp.tsv
    cat taxa.tsv >> tmp.tsv
    sed 's/^/"/' tmp.tsv | sed 's/\$/"/' > taxa_col.tsv
    cut -d, -f 1-3 tree_std.csv > editing_tree.csv
    paste -d"," editing_tree.csv taxa_col.tsv > tree_std.csv 
    cut -d"," -f2 eggnog_std.csv | sed "s/^.//" | sed "s/.\$//" > match_ids.txt
    mv eggnog_std.csv input_eggnog.csv
    string_search.py /sup_data/string_db.json ./match_ids.txt
    sed 's/^/"/' gene_symbols.txt | sed 's/\$/"/' > gene_symbols_01.txt
    sed 's/^/"/' taxa_ids.txt | sed 's/\$/"/' > taxa_ids_01.txt
    paste -d "," input_eggnog.csv gene_symbols_01.txt taxa_ids_01.txt > eggnog_std.csv
    mkdir std_outs
    for file in *std.csv; do
        mv "\$file" ${meta.id}.\${file}
    done 
    mv *std.csv std_outs/ 
    """

}
