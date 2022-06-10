#!/usr/bin/env python

import sys
from Bio import SeqIO

def gen_enzyme_db(enzymedat):
    d = {}
    idline = False
    deline = False
    de = ''
    with open(enzymedat) as enz:
        for line in enz:
            if line.startswith('ID'):
                id = line.strip().split('   ')[1]
                idline = True
            elif line.startswith('DE'):
                if de == '':
                    de = line.strip().split('   ')[1][0].lower() + line.strip().split('   ')[1][1:].rstrip('.')
                else:
                    de += line.strip().split('   ')[1].lower().rstrip('.')
                deline = True
            else:
                if idline == True and deline == True:
                    d[id] = de
                    de = ''
                    idline = False
                    deline = False
                elif idline == False and deline == True:
                    print("DE {} has no ID line".format(de))
                    sys.exit()
                elif idline == True and deline == False:
                    print("ID {} has no DE line".format(de))
                    sys.exit()
    return d

def generate_products(d, ltags, gb):
    pd = {}
    for ltag in ltags:
        for rec in gb:
            for feat in rec.features:
                if feat.type == 'CDS':
                    if feat.qualifiers['locus_tag'] == [ltag]:
                        if 'EC_number' in list(feat.qualifiers):
                            EC = feat.qualifiers['EC_number'][0]
                            product = d[EC]
                            pd[ltag] = product
    return pd

def output_file(pd, outfile):
    with open(outfile, 'w') as of:
        for ltag, product in pd.items():
            of.write("{}\tproduct\t{}\n".format(ltag, product))

def main():
    enzdat = sys.argv[1]
    lt_list = [x.strip().split('\t')[3] for x in open(sys.argv[2])]
    gb = list(SeqIO.parse(sys.argv[3], 'gb'))
    outfile = sys.argv[4]
    enzd = gen_enzyme_db(enzdat)
    prod_d = generate_products(enzd, lt_list, gb)
    output_file(prod_d, outfile)


if __name__ == '__main__':
    main()
    