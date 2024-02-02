#!/bin/bash -ue
cd /usr/local/lib/python3.11/site-packages/
    mkdir data
    cd /
    echo -e "y
y
y" | download_eggnog_data.py --data_dir /usr/local/lib/python3.11/site-packages/data
    emapper.py -i fastas'/'human_chr2_short.fasta -m diamond -o human_test --itype proteins --cpu 6
