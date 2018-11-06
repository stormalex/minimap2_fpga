#!/bin/bash

#./minimap2 -t 4 -W 40 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/pacbio/059-060_128_E01_1_filtered_subreads.fasta.gz > 059.060.128.e01.segfault.sam


#./minimap2 -t 4 -W 20 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/PacBio-HG02818-10files/SRR6056348.fasta > SRR6056348.sam

./minimap2 -t 4 -W 20 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/PacBio-HG02818-10files/SRR6056461.fasta > SRR6056461.sam
