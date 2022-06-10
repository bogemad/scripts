#!/usr/bin/env python

import sys, csv
from Bio import SeqIO


def import_file(file):
    d = {}
    with open(file) as handle:
        for line in handle:
            if line.strip() != '':
                data = line.strip('\n').split()
                if len(data) == 2:
                    d[data[0]] = data[1]
                elif len(data) == 1:
                    d[data[0]] = ''
                else:
                    sys.exit("ERROR: Line has more than 2 or less than 1 elements.")
    return d

def id_largest_ltag(gb):
    m = []
    ltag = ''
    recs = SeqIO.parse(gb, 'gb')
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'gene':
                if ltag == '':
                    ltag = feat.qualifiers['locus_tag'][0].split('_')[0]
                n = int(feat.qualifiers['locus_tag'][0].split('_')[1])
                m.append(n)
    return max(m), ltag

def get_original_locus_tags_and_locations(gb):
    ltagsnlocs = {}
    recs = SeqIO.parse(gb, 'genbank')
    for rec in recs:
        ltagsnlocs[rec.id] = []
        for feat in rec.features:
            if feat.type == 'gene':
                if not feat.qualifiers['locus_tag'][0] in list(ltagsnlocs[rec.id]):
                    if len(feat.qualifiers['locus_tag']) != 1:
                        sys.exit('ERROR: Locus tag list too big {}'.format(feat.qualifiers['locus_tag']))
                    ltagsnlocs[rec.id] += feat.qualifiers['locus_tag']
    return ltagsnlocs

def check_if_all_ltags_present(d, gbd):
    outd = {}
    for l in gbd.values():
        for ltag in l:
            if ltag in list(d.keys()):
                outd[ltag] = d[ltag]
            else:
                print("{} is missing from original".format(ltag))
                outd[ltag] = ''
    return outd

def renumber(d, m, ltag_p, apimit):
    outd = {}
    for k, v in d.items():
        if k in apimit:
            continue
        if v == '':
            m += 1
            nv = ltag_p + '_' + str(m).zfill(6)
        else:
            nv = v
        outd[k] = nv
    for k in apimit:
        if d[k] == '':
            m += 1
            nv = ltag_p + '_' + str(m).zfill(6)
        else:
            nv = d[k]
        outd[k] = nv
    return outd

def main():
    infile = sys.argv[1]
    gb = sys.argv[2]
    old_gb = sys.argv[3]
    outfile = sys.argv[4]
    d = import_file(infile)
    gbd = get_original_locus_tags_and_locations(gb)
    m, tag_prefix = id_largest_ltag(old_gb)
    d = check_if_all_ltags_present(d, gbd)
    apimit = gbd['api'] + gbd['mit']
    outd = renumber(d, m, tag_prefix, apimit)
    with open(outfile, 'w') as outh:
        csvw = csv.writer(outh, delimiter='\t')
        for k, v in outd.items():
            csvw.writerow([k, v])

if __name__ == '__main__':
    main()
    