#!/bin/bash -ue
docker run --rm --name treegrafter -v data:/sample -v outs:/output -v treeGrafter1.01_supplemental:/opt/supplemental ningzhithm/treegrafter:1.01 -f /sample/human_ch2_9genes.fasta -o /output/human_nf.out -cpus 6 -d /opt/supplemental -auto
