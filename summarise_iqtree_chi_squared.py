#!/usr/bin/env python3

import sys, os
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

first = True
big_d = {}
failed = []
all_failed = []

with open(infile) as inha:
    for line in inha:
        if line.startswith('Reading alignment file '):
            if first == True:
                first = False
                filename = os.path.basename(line.split(' ')[3])
                little_d = {}
                print("Working on {}...".format(filename))
                continue
            if False in list(little_d.values()):
                failed.append(filename)
            if list(little_d.values()) == [False] * 11:
                all_failed.append(filename)
            big_d[filename] = pd.Series(little_d)
            little_d = {}
            print("Working on {}...".format(filename))
            filename = os.path.basename(line.split(' ')[3])
        if line.startswith('  ') and line.strip() != '':
            if line.strip().startswith("Gap"):
                continue
            data = line.strip().split()
            try:
                if data[3] == 'failed':
                    little_d[data[1]] = False
                elif data[3] == 'passed':
                    little_d[data[1]] = True
                else:
                    print("ERROR not passed/failed")
                    sys.exit(data[3])
            except Exception as e:
                print("ERROR not passed/failed")
                print("Exception: {}".format(e))
                sys.exit(line)
    df = pd.DataFrame(big_d)
    df.to_csv(outfile)
    print(len(failed))
    print(all_failed)

            



