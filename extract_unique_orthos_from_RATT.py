#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

recs = list(SeqIO.parse(sys.argv[1], 'gb'))
with open(sys.argv[2]) as uniq_gene_file:
    uniq_genes = uniq_gene_file.read().splitlines()


new_recs = []

for rec in recs:
    new_rec = SeqRecord(rec.seq, id=rec.id, name=rec.name, description=rec.description)
    new_rec.annotations = rec.annotations
    for feat in rec.features:
        if 'locus_tag' in list(feat.qualifiers):
            for ltag in feat.qualifiers['locus_tag']:
                if ltag in uniq_genes:
                    new_rec.features.append(feat)
    print('Adding {} features from chromosome {}'.format(len(new_rec.features), new_rec.id))
    new_recs.append(new_rec)
    
SeqIO.write(new_recs, sys.argv[3], 'gb')
