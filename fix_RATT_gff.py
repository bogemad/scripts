#!/usr/bin/env python

import sys, csv

infile = sys.argv[1]
outfile = sys.argv[2]

class Gene:
    def __init__(self):
        self.gene = []
        self.rna = []
        self.exons = []
        self.cds = []
    
    def printgff(self):
        print(self.gene)
        print(self.rna)
        print(self.exons)
        print(self.cds)
    
    def check_output(self):
        if self.rna == []:
            print("RNA missing for gene:")
            print(self.gene)
            return True
        return False
        
    def output_gene(self, handle):
        if self.check_output() == True:
            return
        self.fix_gene()
        self.fix_rna()
        self.fix_exons()
        self.fix_cds()
        gff.fix_nonsense_genes()
        gff.write_gff(handle)
    
    def fix_gene(self):
        qualifiers = self.gene[8].split(';')
        outqual = []
        for qual in qualifiers:
            if qual.startswith('Name='):
                outqual.append('ID=' + qual.replace('Name=',''))
        self.gene = self.gene[:8] + [';'.join(outqual)]
        
    def fix_gene_from_rna(self, rna_row):
        qualifiers = rna_row[8].split(';')
        outqual = []
        for qual in qualifiers:
            if qual.startswith('note=Parent:gene-'):
                outqual.append('Name=' + qual.replace('note=Parent:gene-',''))
        return (rna_row[:2] + ['gene'] + rna_row[3:8] + [';'.join(outqual)])
        
    def fix_rna(self):
        qualifiers = self.rna[8].split(';')
        outqual = []
        for qual in qualifiers:
            if qual.startswith('note=ID:rna-'):
                outqual.append('ID=' + qual.replace('note=ID:rna-',''))
            if qual.startswith('note=Parent:gene-'):
                outqual.append('Parent=' + qual.replace('note=Parent:gene-',''))
            if qual.startswith('product='):
                outqual.append(qual)
        self.rna = self.rna[:8] + [';'.join(outqual)]
    
    def fix_exons(self):
        outexons = []
        for exon in self.exons:
            qualifiers = exon[8].split(';')
            outqual = []
            for qual in qualifiers:
                if qual.startswith('note=ID:exon-'):
                    outqual.append('ID=' + qual.replace('note=ID:exon-',''))
                if qual.startswith('note=Parent:rna-'):
                    outqual.append('Parent=' + qual.replace('note=Parent:rna-',''))
            outexons.append(exon[:8] + [';'.join(outqual)])
        self.exons = outexons
    
    def fix_cds(self):
        outcds = []
        for cds in self.cds:
            qualifiers = cds[8].split(';')
            outqual = []
            for qual in qualifiers:
                if qual.startswith('note=ID:cds-'):
                    outqual.append('ID=' + qual.replace('note=ID:cds-',''))
                if qual.startswith('note=Parent:rna-'):
                    outqual.append('Parent=' + qual.replace('note=Parent:rna-',''))
            outcds.append(cds[:8] + [';'.join(outqual)])
        self.cds = outcds
    
    def fix_nonsense_genes(self):
        if self.gene[6] == '-':
            # self.gene = self.gene[:3] + self.fix_nonsense_start_end(self.gene[3], self.gene[4]) + self.gene[5:]
            # self.rna = self.rna[:3] + self.fix_nonsense_start_end(self.rna[3], self.rna[4]) + self.rna[5:]
            # print(self.exons)
            # outexons = []
            # for exon in self.exons:
                # outexons.append(exon[:3] + self.fix_nonsense_start_end(exon[3], exon[4]) + exon[5:])
            # self.exons = outexons
            # print(self.exons)
            # outcds = []
            # for cds in self.cds:
                # outcds.append(cds[:3] + self.fix_nonsense_start_end(cds[3], cds[4]) + cds[5:])
            # self.cds = outcds
            self.exons.sort(key=lambda x: int(x[3]), reverse=True)
            self.cds.sort(key=lambda x: int(x[3]), reverse=True)
    
    # def fix_nonsense_start_end(self, s, e):
        # h = max(int(s), int(e))
        # l = min(int(s), int(e))
        # return [str(h), str(l)]
    
    def write_gff(self, handle):
        handle.writerow(gff.gene)
        handle.writerow(gff.rna)
        for exon in gff.exons:
            handle.writerow(exon)
        for cds in gff.cds:
            handle.writerow(cds)


with open(sys.argv[1]) as infile, open(sys.argv[2], 'w') as outfile:
    csv_in = csv.reader(infile, delimiter='\t')
    csv_out = csv.writer(outfile, delimiter='\t')
    gff = Gene()
    no_gene = False
    for i, row in enumerate(csv_in):
        if row[0].startswith('#'):
            csv_out.writerow(row)
            continue
        if row[2] == 'gene' and gff.gene != []:
            gff.output_gene(csv_out)
            gff = Gene()
            gff.gene = row
            no_gene = False
        elif row[2] == 'gene':
            gff.gene = row
            no_gene = False
        elif row[2].endswith('RNA'):
            if gff.rna != []:
                gff.output_gene(csv_out)
                gff = Gene()
                gff.gene = gff.fix_gene_from_rna(row)
            gff.rna = row
        elif row[2] == 'exon':
            if no_gene == True:
                continue
            gff.exons.append(row)
        elif row[2] == 'CDS':
            if no_gene == True:
                continue
            gff.cds.append(row)
        else:
            sys.exit('\t'.join(row))
    gff.output_gene(csv_out)
    