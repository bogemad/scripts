#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from collections import defaultdict, OrderedDict


class GeneStruct:
    def __init__(self):
        self.gene = False
        self.rna = False
        self.cds = False
    
    def printall(self):
        print(self.gene)
        print(self.rna)
        print(self.cds)


def find_max_locus_tag(recs):
    ltag_nums = []
    for rec in recs:
        for feat in rec.features:
            if 'locus_tag' in feat.qualifiers:
                ltag = feat.qualifiers['locus_tag'][0]
                ltag_nums.append(int(ltag.split('_')[1]))
                ltag_prefix = ltag.split('_')[0]
    return max(ltag_nums), ltag_prefix

def import_features(recs, max_ltag, ltag_prefix):
    feat_d = defaultdict(dict)
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'gene' or feat.type == 'CDS' or feat.type == 'source' or feat.type.endswith('RNA'):
                if 'locus_tag' in feat.qualifiers:
                    ltag = feat.qualifiers['locus_tag'][0]
                    if not ltag in list(feat_d[rec.id]) and not feat.type == 'source':
                        feat_d[rec.id][ltag] = GeneStruct()
                    if feat.type == 'gene':
                        feat_d[rec.id][ltag].gene = feat
                    elif feat.type.endswith('RNA'):
                        feat_d[rec.id][ltag].rna = feat
                    elif feat.type == 'CDS':
                        feat_d[rec.id][ltag].cds = feat
                    elif feat.type == 'source':
                        feat_d[rec.id]['source'] = feat
                else:
                    if feat.type == 'CDS':
                        max_ltag += 1
                        ltag = ltag_prefix + '_' + str(max_ltag).zfill(6)
                        feat_d[rec.id][ltag] = GeneStruct()
                        feat_d[rec.id][ltag].cds = feat
                        feat_d[rec.id][ltag].cds.qualifiers['locus_tag'] = [ltag]
                        feat_d[rec.id][ltag].rna = SeqFeature(location=feat.location, type='mRNA', strand=feat.strand, qualifiers=OrderedDict([('locus_tag', [ltag]), ('product', ['hypothetical protein'])]))
                        feat_d[rec.id][ltag].gene = SeqFeature(location=FeatureLocation(feat.location.start, feat.location.end, strand=feat.strand), type='gene', strand=feat.strand, qualifiers=OrderedDict([('locus_tag', [ltag])]))                   
                    if feat.type == 'source':
                        feat_d[rec.id]['source'] = feat
    return feat_d

def fix_features(recs, feat_d):
    for rec in recs:
        for ltag in feat_d[rec.id].keys():
            if ltag == 'source':
                continue
            if feat_d[rec.id][ltag].cds != False:
                if feat_d[rec.id][ltag].rna == False or feat_d[rec.id][ltag].gene == False:
                    feat = feat_d[rec.id][ltag].cds
                    feat_d[rec.id][ltag].rna = SeqFeature(location=feat.location, type='mRNA', strand=feat.strand, qualifiers=OrderedDict([('locus_tag', [ltag]), ('product', ['hypothetical protein'])]))
                    feat_d[rec.id][ltag].gene = SeqFeature(location=FeatureLocation(feat.location.start, feat.location.end, strand=feat.strand), type='gene', strand=feat.strand, qualifiers=OrderedDict([('locus_tag', [ltag])]))
                if feat_d[rec.id][ltag].rna.location != feat_d[rec.id][ltag].cds.location:
                    feat_d[rec.id][ltag].rna.location = feat_d[rec.id][ltag].cds.location
                    feat_d[rec.id][ltag].gene.location = FeatureLocation(feat_d[rec.id][ltag].cds.location.start, feat_d[rec.id][ltag].cds.location.end, strand=feat_d[rec.id][ltag].cds.strand)
                if not 'translation' in list(feat_d[rec.id][ltag].cds.qualifiers):
                    try:
                        feat_d[rec.id][ltag].cds.qualifiers['translation'] = [str(feat_d[rec.id][ltag].cds.translate(rec).seq)]
                    except:
                        print(rec.id, feat_d[rec.id][ltag].cds.location)
                if not 'codon_start' in list(feat_d[rec.id][ltag].cds.qualifiers):
                    feat_d[rec.id][ltag].cds.qualifiers['codon_start'] = ['1']
                if not 'protein_id' in list(feat_d[rec.id][ltag].cds.qualifiers):
                    feat_d[rec.id][ltag].cds.qualifiers['protein_id'] = ['{}-T1'.format(feat_d[rec.id][ltag].cds.qualifiers['locus_tag'][0])]
    return feat_d

def reorder_features(recs, feat_d):
    new_recs = []
    for rec in recs:
        rec.features = []
        rec.features.append(feat_d[rec.id].pop('source'))
        for ltag, gs in feat_d[rec.id].items():
            rec.features.append(gs.gene)
            rec.features.append(gs.rna)
            if not gs.cds == False:
                rec.features.append(gs.cds)
        new_recs.append(rec)
    return new_recs

def check_for_empty_SeqFeatures(feat_d):
    ltags = []
    for rec_id, feats in feat_d.items():
        for ltag, gs in feats.items():
            if gs.gene == False or gs.rna == False or gs.cds == False:
                ltags.append(ltag)
    return ltags

def main():
    recs = list(SeqIO.parse(sys.argv[1], 'gb'))
    max_ltag, ltag_prefix = find_max_locus_tag(recs)
    feat_d = import_features(recs, max_ltag, ltag_prefix)
    feat_d = fix_features(recs, feat_d)
    new_recs = reorder_features(recs, feat_d)
    ltags = check_for_empty_SeqFeatures(feat_d)
    SeqIO.write(new_recs, sys.argv[2], 'gb')

if __name__ == '__main__':
    main()
