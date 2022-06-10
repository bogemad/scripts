#!/usr/bin/env python

import sys
from Bio import SeqIO

infile = sys.argv[1]
outfile = sys.argv[2]

recs = list(SeqIO.parse(infile, 'gb'))

for rec in recs:
    if rec.id == 'api' or rec.id == 'mit':
        continue
    for feat in rec.features:
        for qual in list(feat.qualifiers):
            if not (qual == 'product' or qual == 'locus_tag' or qual == 'transcript_id' or qual == 'protein_id' or qual == 'codon_start' or qual == 'translation'):
                del feat.qualifiers[qual]
        if feat.type == 'CDS' or feat.type == 'mRNA':
            feat.qualifiers['product'] = ["hypothetical protein"]

SeqIO.write(recs, outfile, 'gb')
