#!/usr/bin/env python

import sys, os
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
import pandas as pd

def total_predicted_genes(recs):
    genes = []
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'gene':
                genes.append(feat)
    return len(genes)

def percentage_coding(recs):
    all_seq = ''
    coding_seq = ''
    for rec in recs:
        all_seq += str(rec.seq)
        for feat in rec.features:
            if feat.type == 'CDS':
                coding_seq += str(feat.extract(rec.seq))
    return round(float(len(coding_seq))/float(len(all_seq))*100,2)

def total_sequence_length(recs):
    all_seq = ''
    for rec in recs:
        all_seq += str(rec.seq)
    return len(all_seq)

def percentage_gc(recs):
    all_seq = False
    for rec in recs:
        if all_seq == False:
            all_seq = rec.seq
        else:
            all_seq += rec.seq
    return round(GC(all_seq),2)

def longest_gene(recs):
    longest = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'gene':
                longest = max(longest, len(feat))
    return longest

def shortest_gene(recs):
    shortest = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'gene':
                if shortest == 0:
                    shortest = len(feat)
                else:
                    shortest = min(shortest, len(feat))
    return shortest

def total_gene_length(recs):
    l = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'gene':
                l += len(feat)
    return l

def average_gene_length(num, len):
    return round(float(len)/float(num), 0)

def average_CDS_length(recs):
    lenCDS = 0
    numCDS = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'CDS':
                lenCDS += len(feat)
                numCDS += 1
    return round(float(lenCDS)/float(numCDS), 0)

def gene_density(num_genes, len_seq):
    return (round(float(num_genes)/float(len_seq)*1000000, 2))

def percentage_coding_genes_with_introns(recs):
    total = 0
    intron = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'CDS':
                total += 1
                if len(feat.location.parts) > 1:
                    intron += 1
    return (round(float(intron)/float(total)*100, 1))

def total_exon_length(recs):
    l = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                l += len(feat)
    return l

def total_number_of_exons(recs):
    i = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                for loc in feat.location.parts:
                    i += 1
    return i

def longest_exon(recs):
    longest = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                for loc in feat.location.parts:
                    longest = max(longest, len(loc))
    return longest

def shortest_exon(recs):
    shortest = 10000000
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                for loc in feat.location.parts:
                    shortest = min(shortest, len(loc))
    return shortest

def average_exon_length(n, l):
    return round(float(l)/float(n), 1)

def exon_percent_gc(recs):
    exon_seq = False
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                if exon_seq == False:
                    exon_seq = feat.extract(rec.seq)
                else:
                    exon_seq += feat.extract(rec.seq)
    return round(GC(exon_seq),2)

def total_number_of_introns(recs):
    i = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                if type(feat.location) == CompoundLocation:
                    i += (len(feat.location.parts) - 1)
    return i

def convert_cpnd_loc_to_introns(cl, strand):
    introns = []
    last_i = len(cl.parts) - 1
    if last_i == 0:
        return []
    for i, loc in enumerate(cl.parts):
        if strand == 1:
            if i == 0:
                start = loc.end
            elif i == last_i:
                end = loc.start
                introns.append(FeatureLocation(start, end, strand))
            else:
                end = loc.start
                introns.append(FeatureLocation(start, end, strand))
                start = loc.end
        elif strand == -1:
            if i == 0:
                end = loc.start
            elif i == last_i:
                start = loc.end
                introns.append(FeatureLocation(start, end, strand))
            else:
                start = loc.end
                introns.append(FeatureLocation(start, end, strand))
                end = loc.start
        else:
            sys.stderr.write(
                    "No strand indicated. Assuming +\n"
                )
    if len(introns) == 1:
        return introns[0]
    return CompoundLocation(introns)

def longest_intron(recs):
    longest = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                if type(feat.location) == CompoundLocation:
                    locs = convert_cpnd_loc_to_introns(feat.location, feat.strand)
                    if locs == []:
                        continue
                    for loc in locs.parts:
                        longest = max(longest, len(loc))
    return longest

def shortest_intron(recs):
    shortest = 10000000
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                if type(feat.location) == CompoundLocation:
                    for loc in convert_cpnd_loc_to_introns(feat.location, feat.strand).parts:
                      shortest = min(shortest, len(loc))
    return shortest

def total_intron_length(recs):
    l = 0
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                l += len(convert_cpnd_loc_to_introns(feat.location, feat.strand))
    return l

def average_intron_length(n, l):
    return round(float(l)/float(n), 1)

def average_introns_per_gene(ng, ni):
    return round(float(ni)/float(ng), 1)

def intron_percent_gc(recs):
    intron_seq = False
    for rec in recs:
        for feat in rec.features:
            if feat.type.endswith('RNA'):
                locs = convert_cpnd_loc_to_introns(feat.location, feat.strand)
                if locs == []:
                    continue
                if intron_seq == False:
                    intron_seq = locs.extract(rec.seq)
                else:
                    intron_seq += locs.extract(rec.seq)
    return round(GC(intron_seq),2)

def convert_genes_to_intergenes(genes, l):
    intergenes = []
    last_i = len(genes) - 1
    for i, gene in enumerate(genes):
        if i == 0:
            start = 0
            end = gene.location.start
            if end > start:
                intergenes.append(FeatureLocation(start, end, 1))
            start = gene.location.end
        elif i == last_i:
            end = gene.location.start
            if end > start:
                intergenes.append(FeatureLocation(start, end, 1))
            if l > gene.location.end:
                intergenes.append(FeatureLocation(gene.location.end, l, 1))
        else:
            end = gene.location.start
            if end > start:
                intergenes.append(FeatureLocation(start, end, 1))
            start = gene.location.end
    if len(intergenes) == 0:
        return []
    if len(intergenes) == 1:
        return intergenes
    return CompoundLocation(intergenes)

def total_intergenic_regions(recs):
    l = 0
    for rec in recs:
        genes = []
        for feat in rec.features:
            if feat.type == 'gene':
                genes.append(feat)
        intergenes = convert_genes_to_intergenes(genes, len(rec))
        if len(intergenes) == 0:
            continue
        l += len(intergenes.parts)
    return l

def longest_intergenic_region(recs):
    longest = 0
    for rec in recs:
        genes = []
        for feat in rec.features:
            if feat.type == 'gene':
                genes.append(feat)
        intergenes = convert_genes_to_intergenes(genes, len(rec))
        if len(intergenes) == 0:
            continue
        for intergene in intergenes.parts:
            longest = max(len(intergene), longest)
    return longest

def shortest_intergenic_region(recs):
    shortest = 10000000
    for rec in recs:
        genes = []
        for feat in rec.features:
            if feat.type == 'gene':
                genes.append(feat)
        intergenes = convert_genes_to_intergenes(genes, len(rec))
        if len(intergenes) == 0:
            continue
        for intergene in intergenes.parts:
            shortest = min(len(intergene), shortest)
    return shortest

def total_intergenic_length(recs):
    l = 0
    for rec in recs:
        genes = []
        for feat in rec.features:
            if feat.type == 'gene':
                genes.append(feat)
        intergenes = convert_genes_to_intergenes(genes, len(rec))
        if len(intergenes) == 0:
            continue
        for intergene in intergenes.parts:
            l += len(intergene)
    return l

def average_intergenic_length(l, n):
    return round(float(l)/float(n), 1)

def intergenic_percent_gc(recs):
    intergene_seq = False
    for rec in recs:
        genes = []
        for feat in rec.features:
            if feat.type == 'gene':
                genes.append(feat)
        intergenes = convert_genes_to_intergenes(genes, len(rec))
        if len(intergenes) == 0:
            continue
        if intergene_seq == False:
            intergene_seq = intergenes.extract(rec.seq)
        else:
            intergene_seq += intergenes.extract(rec.seq)
    return round(GC(intergene_seq),2)

def main():
    gbs = sys.argv[1:-1]
    output = sys.argv[-1]
    dfd = {}
    for gb in gbs:
        d = {}
        recs = list(SeqIO.parse(gb, 'gb'))
        d["Total predicted genes"] = total_predicted_genes(recs)
        d["Percentage coding"] = percentage_coding(recs)
        d['Total sequence length'] = total_sequence_length(recs)
        d['Percentage GC'] = percentage_gc(recs)
        d['Longest gene'] = longest_gene(recs)
        d['Shortest gene'] = shortest_gene(recs)
        d['Total gene length'] = total_gene_length(recs)
        d['Average gene length'] = average_gene_length(d["Total predicted genes"], d['Total gene length'])
        d['Average gene CDS length'] = average_CDS_length(recs)
        d['Gene density (genes/Mb)'] = gene_density(d["Total predicted genes"], d['Total sequence length'])
        d['Percentage coding genes with introns'] = percentage_coding_genes_with_introns(recs)
        d['Total exon length'] = total_exon_length(recs)
        d['Total number of exons'] = total_number_of_exons(recs)
        d['Longest exon'] = longest_exon(recs)
        d['Shortest exon'] = shortest_exon(recs)
        d['Average exon length'] = average_exon_length(d['Total number of exons'], d['Total exon length'])
        d['Exon percentage GC'] = exon_percent_gc(recs)
        d['Total intron length'] = total_intron_length(recs)
        d['Total number of introns'] = total_number_of_introns(recs)
        d['Longest intron'] = longest_intron(recs)
        d['Shortest intron'] = shortest_intron(recs)
        d['Average intron length'] = average_intron_length(d['Total number of introns'], d['Total intron length'])
        d['Average introns per gene'] = average_introns_per_gene(d["Total predicted genes"], d['Total number of introns'])
        d['Intron percent GC'] = intron_percent_gc(recs)
        d['Total intergenic length'] = total_intergenic_length(recs)
        d['Total intergenic regions'] = total_intergenic_regions(recs)
        d['Longest intergenic region'] = longest_intergenic_region(recs)
        d['Shortest intergenic region'] = shortest_intergenic_region(recs)
        d['Average intergenic length'] = average_intergenic_length(d['Total intergenic length'], d['Total intergenic regions'])
        d['Intergenic percent GC'] = intergenic_percent_gc(recs)
        dfd[os.path.basename(gb)] = pd.Series(d)
    df = pd.DataFrame(dfd)
    df.to_csv(output)
        
if __name__ == '__main__':
    main()
