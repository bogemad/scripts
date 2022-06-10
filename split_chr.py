#!/usr/bin/env python3

import sys, os
from Bio import SeqIO
from collections import defaultdict



def id_files(indir):
    gff_l = []
    fasta = ''
    gbk = ''
    for item in os.listdir(indir):
        print(item)
        if os.path.isfile(os.path.join(indir, item)):
            if item.endswith('.scaffolds.fa'):
                fasta = os.path.join(indir, item)
            if item.endswith('.gbk'):
                gbk = os.path.join(indir, item)
            if item.endswith('.gff') or item.endswith('.gff3'):
                gff_l.append(os.path.join(indir, item))
    print("fasta = {}".format(fasta))
    print("gbk = {}".format(gbk))
    print("gffs = {}".format(gff_l))
    return fasta, gbk, gff_l

def split_fasta(indir, fasta):
    recs = SeqIO.parse(fasta, 'fasta')
    for rec in recs:
        out = os.path.join(indir, rec.id)
        os.makedirs(out, exist_ok = True)
        SeqIO.write([rec], os.path.join(out, '{}.fa'.format(rec.id)), 'fasta')

def split_gbk(indir, gbk):
    recs = SeqIO.parse(gbk, 'genbank')
    for rec in recs:
        out = os.path.join(indir, rec.id)
        os.makedirs(out, exist_ok = True)
        SeqIO.write([rec], os.path.join(out, '{}.gbk'.format(rec.id)), 'genbank')

def split_gffs(indir, gff_l):
    for gff in gff_l:
        lines_d = defaultdict(list)
        with open(gff) as gffh:
            for line in gffh:
                if line.startswith('#') or line.strip() == '':
                    continue
                chrom = line.split()[0]
                lines_d[chrom].append(line.strip())
        for chrom in lines_d.keys():
            out = os.path.join(indir, chrom)
            with open(os.path.join(out, os.path.basename(gff)), 'w') as outh:
                outh.write('\n'.join(lines_d[chrom]) + '\n')

def main():
    indir = sys.argv[1]
    print("Splitting files in {}".format(indir))
    fasta, gbk, gff_l = id_files(indir)
    split_fasta(indir, fasta)
    split_gbk(indir, gbk)
    split_gffs(indir, gff_l)

if __name__ == '__main__':
    main()
