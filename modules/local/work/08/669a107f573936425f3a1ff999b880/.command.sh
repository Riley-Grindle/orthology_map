#!/bin/bash -ue
emapper.py -i fastas/axo*.fasta -m diamond -o human_test --itype proteins --cpu 6
