#!/bin/bash -ue
docker run --rm --name treegrafter -v /Users/rgrindle/Desktop/mdibl/wd/ortho_finder/tool_testing/panther/nextflow_test/data/:/sample -v /Users/rgrindle/Desktop/mdibl/wd/ortho_finder/tool_testing/panther/nextflow_test/treeGrafter1.01_supplemental:/opt/supplemental ningzhithm/treegrafter:1.01 -f /sample/human_ch2_9genes.fasta -o /output/human_nf.out -cpus 6 -d /opt/supplemental -auto