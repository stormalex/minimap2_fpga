#!/usr/bin/env python2
# -*- coding :utf-8 -*-
import sys,os
import shutil as sh
import time
# /home/ljb/minimap2/mk.sh
'''
path="/home/data/data_minimap/PacBio-HG02818-10files"
with open(path+"/read") as f:
	for line in f:
		cmd='time ./minimap2 -t 4 -W 40 --no-kalloc -ax map-pb /home/data/ref/ucsc.hg19.pb.mmi %s/%s >%s.sam ' %(path,line[:-1],line[:-1])
		os.system(cmd)
		print cmd

'''
path="/home/data/data_minimap/pacbio"
with open(path+"/read") as f:
	for line in f:
		cmd='time ./minimap2 -t 4 -W 40 --no-kalloc -ax map-pb /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mappb.mmi %s/%s >%s.sam ' %(path,line[:-1],line[:-1])
		os.system(cmd)
		print cmd

path="/home/data/data_minimap/nanopore"
with open(path+"/read") as f:
	for line in f:
		print line;
		cmd='time ./minimap2 -t 4 -W 40 --no-kalloc -ax map-ont /home/data/ref/GCA_000001405.27_GRCh38.p12_genomic.fna.mapont.mmi %s/%s >%s.sam ' %(path,line[:-1],line[:-1])
		os.system(cmd)
		print cmd


path="/home/data/data_minimap/pacbio"
with open(path+"/read") as f:
        for line in f:
                print line;
                cmd='time ./minimap2 -t 4 -W 40 --no-kalloc -ax map-pb /home/data/ref/ucsc.hg19.pb.mmi %s/%s >%s.sam ' %(path,line[:-1],line[:-1])
		os.system(cmd)
                print cmd


path="/home/data/data_minimap/nanopore"
with open(path+"/read") as f:
        for line in f:
                print line;                
		cmd='time ./minimap2 -t 4 -W 40 --no-kalloc -ax map-ont /home/data/ref/ucsc.hg19.mapont.mmi %s/%s >%s.sam ' %(path,line[:-1],line[:-1])
		os.system(cmd)
                print cmd
