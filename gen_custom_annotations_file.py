from Bio import SeqIO
import sys, os

gbk = sys.argv[1]
annos = sys.argv[2]

recs = list(SeqIO.parse(gbk, 'gb'))
outfile = open(annos, 'w')

for rec in recs:
    for feat in rec.features:
        if feat.type.endswith('RNA'):
            ltag = feat.qualifiers['locus_tag'][0]
            product = feat.qualifiers['product'][0]
            if product == 'uncharacterized protein':
                product = "hypothetical protein"
            if not product == "hypothetical protein":
                outline = '{}\tproduct\t{}\n'.format(ltag, product)
                outfile.write(outline)

outfile.close()