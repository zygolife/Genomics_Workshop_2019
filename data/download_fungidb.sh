#!/usr/bin/bash
if [ ! -f genomic/Ncrassa.gff ]; then
curl -o genomic/Ncrassa.gff https://fungidb.org/common/downloads/Current_Release/NcrassaOR74A/gff/data/FungiDB-45_NcrassaOR74A.gff
fi
