#!/bin/bash -ue
docker run --rm --name treegrafter -v data:/sample -v outs:/output -v treeGrafter1.01_supplemental:/opt/supplemental ningzhithm/treegrafter:1.01 -f /sample/human_ch2.fasta -o /output/human_nf.out -d /opt/supplemental -auto
