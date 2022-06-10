#!/usr/bin/env python

import sys, re

from Bio import SeqIO

recs = SeqIO.parse(sys.argv[1], 'gb')
tbl = sys.argv[2]
out_tbl = sys.argv[3]
keyd = {}

for rec in recs:
    for feat in rec.features:
        if feat.type == 'CDS':
            keyd[feat.qualifiers['locus_tag'][0]] = feat.qualifiers['protein_id'][0]

with open(tbl) as inh, open(out_tbl, 'w') as outh:
    for i, line in enumerate(inh):
        if line.strip().startswith('protein_id'):
            search_term = r"MAC[JK]_[0-9]{{6}}"
            matches = re.findall(search_term, line)
            for match in matches:
                print("Replaced {} with {} on line {}".format(match, keyd[match], i+1))
                line = line.strip() + '|gb|{}\n'.format()
            outh.write(line)

