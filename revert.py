#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys, os
import shutil

def revert(fpath):
    with open(fpath) as f:
        with open(fpath + '.r.txt', 'w') as f2:
            while True:
                line=f.readline()
                if not line:
                    break
                line = line.strip()
                for i in range(64):
                    pos = (63 - i) * 2
                    f2.write(line[pos:pos+2])
                f2.write('\n')

revert('idxb.txt')
revert('idxh.txt')
revert('idxp.txt')
revert('idxv.txt')
os.system('mv idxb.txt.r.txt idxb.txt')
os.system('mv idxh.txt.r.txt idxh.txt')
os.system('mv idxp.txt.r.txt idxp.txt')
os.system('mv idxv.txt.r.txt idxv.txt')

