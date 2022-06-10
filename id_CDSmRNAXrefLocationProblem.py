#!/usr/bin/env python
from Bio import SeqIO

import sys

gbk = SeqIO.parse(sys.argv[1], 'gb')

for rec in gbk:
    for feat in rec.features:
        if feat == 'mRNA':
            mRNA_start = feat.location.start
            mRNA_end = feat.location.end
        if feat == 'CDS':
            if feat.location.start != mRNA_start or feat.location.end != mRNA_end:
                print(feat.qualifiers['locus_tag'][0])
