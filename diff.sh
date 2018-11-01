#!/usr/bin/env python2
# -*- coding :utf-8 -*-
import sys,os
import shutil as sh
import time
# /home/ljb/minimap2/mk.sh

path="/home/data/data_minimap/pacbio"
with open(path+"/readdiff") as f:
	for line in f:
		cmd='diff  %s/%s.hg38.sam %s.sam >%s.diff ' %(path,line[:-1],line[:-1],line[:-1])
		os.system(cmd)
		print cmd

'''
path="/home/data/data_minimap/nanopore"
with open(path+"/read") as f:
	for line in f:
		print line;
		cmd='time diff  %s/%s %s.sam >%s.diff ' %(path,line[:-1],line[:-1],line[:-1])
		os.system(cmd)
		print cmd


path="/home/data/data_minimap/pacbio"
with open(path+"/read") as f:
        for line in f:
                print line;
		cmd='time diff  %s/%s %s.sam >%s.diff ' %(path,line[:-1],line[:-1],line[:-1])
		os.system(cmd)
                print cmd


path="/home/data/data_minimap/nanopore"
with open(path+"/read") as f:
        for line in f:
                print line;                
		cmd='time diff  %s/%s %s.sam >%s.diff ' %(path,line[:-1],line[:-1],line[:-1])
		os.system(cmd)
                print cmd
'''
