#!/usr/bin/env python

import sys, re, csv
from Bio import SeqIO

def import_key(keyfile):
    keyd = {}
    with open(keyfile) as keyh:
        csvread = csv.reader(keyh, delimiter='\t')
        for row in csvread:
            keyd[row[0]] = row[1]
    return keyd

def get_locus_tag(keyfile):
    with open(keyfile) as keyh:
        return keyh.readline().split('\t')[0].split('_')[0]

def convertfile(infile, keyd, pid_key, locus_tag, outfile):
    with open(infile) as inh, open(outfile, 'w') as outh:
        for i, line in enumerate(inh):
            search_term = r"MAC[JK]_[0-9]{{6}}".format(locus_tag)
            matches = re.findall(search_term, line)
            for match in matches:
                print("Replaced {} with {} on line {}".format(match, keyd[match], i+1))
                line = line.replace(match, keyd[match])
                if line.strip().startswith('protein_id') and match in list(pid_key):
                    line = line.rstrip() + '|gb|{}\n'.format(pid_key[match])
            outh.write(line)

def get_protein_id_key(old_gb):
    keyd = {}
    recs = SeqIO.parse(old_gb, 'gb')
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'CDS':
                keyd[feat.qualifiers['locus_tag'][0]] = feat.qualifiers['protein_id'][0]
    return keyd




def main():
    keyfile = sys.argv[1]
    old_gb = sys.argv[2]
    infile = sys.argv[3]
    outfile = sys.argv[4]
    key = import_key(keyfile)
    pid_key = get_protein_id_key(old_gb)
    ltag = get_locus_tag(keyfile)
    convertfile(infile, key, pid_key, ltag, outfile)

if __name__ == '__main__':
    main()
