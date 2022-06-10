#!/usr/bin/env python3

import sys, csv

ogpath = sys.argv[1]
outpath = sys.argv[2]

d = {}

with open(ogpath) as ogfile:
    ogreader = csv.reader(ogfile, delimiter='\t')
    for row in ogreader:
        keys = row[1].split(', ')
        if keys[0] == '':
            continue
        for key in keys:
            d[key] = row[2]

with open(outpath, 'w') as outfile:
    ogwriter = csv.writer(outfile, delimiter='\t')
    for k in sorted(list(d)):
        ogwriter.writerow([k, d[k]])
    
