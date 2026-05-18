#!/bin/bash

conda activate phylo

FILENAME="$(basename "$1")"
FILEPATH="$(dirname "$1")"

mkdir -p ${FILEPATH}/260317_${FILENAME}.dnds
cd ${FILEPATH}/260317_${FILENAME}.dnds

scp ${FILEPATH}/${FILENAME}.dnds/aa.maf.phy .
FastTree -nt -gtr aa.maf.phy > tree.nwk

cat > codeml_m1.ctl << EOF
seqfile  = aa.maf.phy
treefile = tree.nwk
outfile  = dnds.res.m0.out

runmode = 0
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
fix_omega = 0
omega = 1
EOF

codeml codeml_m1.ctl > 260317_codeml_m1.log
