#!/bin/bash

#./minimap2 -t 4 -W 40 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/pacbio/059-060_128_E01_1_filtered_subreads.fasta.gz > 059.060.128.e01.segfault.sam


#./minimap2 -t 4 -W 20 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/PacBio-HG02818-10files/SRR6056348.fasta > SRR6056348.sam

#./minimap2 -t 4 -W 20 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/PacBio-HG02818-10files/SRR6056461.fasta > SRR6056461.sam

./minimap2 -t 8 -W 20 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/pacbio/066_133_A01_1_filtered_subreads.fasta.gz > 066_133_A01_1_filtered_subreads.sam
exit
for ((t_num = 1; t_num <= 8; t_num++))
do
	for ((W_num = 1; W_num <= 30; W_num++))
	do
		echo "======  $t_num   $W_num  ==========="
		./minimap2 -t $t_num -W $W_num --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi /home/data/data_minimap/pacbio/066_133_A01_1_filtered_subreads.fasta.gz > a.sam
	done	
done

