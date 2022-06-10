import sys, csv, os

#run as python merge_RATT_manual_into_funannotate.py manual_annotation_list funannotate_gff RATT_gff manual_gff genious_gff rrna_gff output_gff

from collections import defaultdict
import pprint

class GFF:
    def __init__(self, GFF_file, keyfield):
        with open(GFF_file) as gff:
            self.features = {}
            first_gene = True
            line_list = []
            i = 1
            for line in gff:
                if line.startswith('#') or line.strip() == '':
                    continue
                gffline = self.GFFLine(line)
                if gffline.type == 'gene':
                    if first_gene == True:
                        first_gene = False
                    else:
                        self.features[genekey] = self.GeneStruct(line_list)
                        line_list = []
                    if keyfield == 'increment':
                        genekey = i
                    else:
                        genekey = gffline.attributes[keyfield][0]
                    i += 1
                line_list.append(gffline)
            self.features[genekey] = self.GeneStruct(line_list)
    
    class GeneStruct:
        def __init__(self, line_list):
            self.exons = []
            self.CDS = []
            for line in line_list:
                if line.type == 'gene':
                    self.gene = line
                elif line.type == 'mRNA':
                    self.rna_type = 'mRNA'
                    self.rna = line
                elif line.type == 'tRNA':
                    self.rna_type = 'tRNA'
                    self.rna = line
                elif line.type == 'rRNA':
                    self.rna_type = 'rRNA'
                    self.rna = line
                elif line.type == 'exon':
                    self.exons.append(line)
                elif line.type == 'CDS':
                    self.CDS.append(line)
                else:
                    print("Unrecognised line type: {}".format(line))
                    sys.exit(1)
    
    class GFFLine:
        def __init__(self, line):
            data = line.strip().split('\t')
            self.seqid = data[0]
            self.source = data[1]
            self.type = data[2]
            self.start = int(data[3])
            self.end = int(data[4])
            self.score = data[5]
            self.strand = data[6]
            self.phase = data[7]
            self.attributes = defaultdict(list)
            for attr in data[8].split(';'):
                if not attr.strip() == '':
                    k, v  = attr.split('=')
                    self.attributes[k].append(v)

class EmptyGFF:
    def __init__(self):
        self.features = {}

    def write_line(self, obj, out):
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                                                obj.seqid, 
                                                                obj.source, 
                                                                obj.type, 
                                                                str(obj.start), 
                                                                str(obj.end), 
                                                                obj.score, 
                                                                obj.strand, 
                                                                obj.phase, 
                                                                self.gen_attr_field(obj.attributes)
                                                                ))
    
    def gen_attr_field(self, attr):
        attr_f = ''
        for k, l in attr.items():
            for v in l:
                attr_f += '{}={};'.format(k, v)
        return attr_f
    
    def write(self, outfile):
        with open(outfile, 'w') as outh:
            for feat in self.features.values():
                self.write_line(feat.gene, outh)
                self.write_line(feat.rna, outh)
                for exon in feat.exons:
                    self.write_line(exon, outh)
                for cds in feat.CDS:
                    self.write_line(cds, outh)


def import_change_list(man_ann_file):
    change_list = []
    with open(man_ann_file) as csvh:
        csvr = csv.reader(csvh, delimiter=',')
        for row in csvr:
            tot = []
            mac = []
            for item in row:
                if item.startswith('TOT'):
                    tot.append(item)
                elif item.startswith('MAC'):
                    mac.append(item)
            if len(tot) == 0:
                change_list.append(mac)
            else:
                for item in mac:
                    change_list.append(tot + [item])
    return change_list

def check_ltag_in_change_list(ll, ltag):
    for l in ll:
        if l[-1] == ltag:
            return l
    return []

def change_attr(attrs, old, new):
    new_attrs = []
    for attr in attrs:
        new_attr = attr.replace(old, new)
        new_attrs.append(new_attr)
    return new_attrs

def fix_attributes(gene, ltag, new_ltag):
    out_gene = gene
    out_gene.gene.attributes['ID'] = change_attr(gene.gene.attributes['ID'], ltag, new_ltag)
    out_gene.rna.attributes['ID'] = change_attr(gene.rna.attributes['ID'], ltag, new_ltag)
    out_gene.rna.attributes['Parent'] = change_attr(gene.rna.attributes['Parent'], ltag, new_ltag)
    new_exons = []
    new_CDS = []
    for exon in gene.exons:
        exon.attributes['ID'] = change_attr(exon.attributes['ID'], ltag, new_ltag)
        exon.attributes['Parent'] = change_attr(exon.attributes['Parent'], ltag, new_ltag)
        new_exons.append(exon)
    out_gene.exons = new_exons
    for cds in gene.CDS:
        cds.attributes['ID'] = change_attr(cds.attributes['ID'], ltag, new_ltag)
        cds.attributes['Parent'] = change_attr(cds.attributes['Parent'], ltag, new_ltag)
        new_CDS.append(cds)
    out_gene.CDS = new_CDS
    return out_gene

#flow check if ID in replace_list --> (yep) has it been replaced already (yep) --> add_gene (num) --> remove_gene and add_gene
def fix_ratt_note(note_list, startkey):
    for note in note_list:
        if note.startswith(startkey):
            return [note.replace(startkey, "Annotation sourced from ")]

def fix_ratt_exonsorcds(ecs, new_ltag, strand, ec_type):
    if strand == '-':
        ecs = sorted(ecs, key=lambda x: x.start, reverse=True)
    new_ecs = []
    for i, ec in enumerate(ecs):
        null = ec.attributes.pop('product', None)
        ec.attributes['Parent'] = [new_ltag+'-T1']
        if ec_type == 'cds':
            null = ec.attributes.pop('Name', None)
            ec.attributes['note'] = fix_ratt_note(ec.attributes['note'], 'ID:cds-')
            ec.attributes['ID'] = [new_ltag+'-T1.cds']
        if ec_type == 'exon':
            ec.attributes['ID'] = [new_ltag+'-T1.exon{}'.format(str(i+1))]
            ec.attributes['note'] = fix_ratt_note(ec.attributes['note'], 'ID:exon-')
        new_ecs.append(ec)
    return new_ecs

def fix_geneious_exonsorcds(ecs, new_ltag, strand, ec_type):
    new_ecs = []
    for i, ec in enumerate(ecs):
        ec.attributes = {}
        ec.attributes['Parent'] = [new_ltag+'-T1']
        if ec_type == 'cds':
            ec.attributes['ID'] = [new_ltag+'-T1.cds']
        if ec_type == 'exon':
            ec.attributes['ID'] = [new_ltag+'-T1.exon{}'.format(str(i+1))]
        new_ecs.append(ec)
    return new_ecs

def fix_ratt(gene, ltag, new_ltag):
    out_gene = gene
    #fix gene
    null = out_gene.gene.attributes.pop('note', None)
    null = out_gene.gene.attributes.pop('Name', None)
    out_gene.gene.attributes['ID'] = [new_ltag]
    #fix rna
    null = out_gene.rna.attributes.pop('Name', None)
    out_gene.rna.attributes['note'] = fix_ratt_note(gene.rna.attributes['note'], 'ID:rna-')
    out_gene.rna.attributes['ID'] = [new_ltag+'-T1']
    out_gene.rna.attributes['Parent'] = [new_ltag]
    out_gene.rna.attributes['product'] = ['hypothetical protein']
    #fix exons & cds
    out_gene.exons = fix_ratt_exonsorcds(gene.exons, new_ltag, gene.gene.strand, 'exon')
    out_gene.CDS = fix_ratt_exonsorcds(gene.CDS, new_ltag, gene.gene.strand, 'cds')
    return out_gene

def fix_geneious(gene, new_ltag):
    out_gene = gene
    if gene.rna_type == 'tRNA':
        product = ['tRNA-{}'.format(gene.rna.attributes['isotype'][0])]
    elif gene.rna_type == 'rRNA':
        product = gene.rna.attributes['product']
    else:
        if 'product' in list(gene.rna.attributes):
            product = gene.rna.attributes['product']
        else:
            product = gene.CDS[0].attributes['Name']
    #fix gene
    out_gene.gene.attributes = {}
    out_gene.gene.attributes['ID'] = [new_ltag]
    #fix rna
    out_gene.rna.attributes = {}
    out_gene.rna.attributes['ID'] = [new_ltag+'-T1']
    out_gene.rna.attributes['Parent'] = [new_ltag]
    out_gene.rna.attributes['product'] = product
    #fix exons & cds
    out_gene.exons = fix_geneious_exonsorcds(gene.exons, new_ltag, gene.gene.strand, 'exon')
    out_gene.CDS = fix_geneious_exonsorcds(gene.CDS, new_ltag, gene.gene.strand, 'cds')
    return out_gene

def fix_rrna(gene, new_ltag):
    product = gene.rna.attributes['product']
    out_gene = gene
    out_gene.gene.attributes = {}
    out_gene.gene.attributes['ID'] = [new_ltag]
    out_gene.gene.attributes['Name'] = gene.rna.attributes['Name']
    out_gene.rna.attributes = {}
    out_gene.rna.attributes['ID'] = [new_ltag+'-T1']
    out_gene.rna.attributes['Parent'] = [new_ltag]
    out_gene.rna.attributes['product'] = product
    for i, exon in enumerate(gene.exons):
        exon.attributes = {}
        exon.attributes['ID'] = [new_ltag+'-T1.exon{}'.format(str(i+1))]
        exon.attributes['Parent'] = [new_ltag+'-T1']
    return out_gene

def main():
    man_ann_file = sys.argv[1]
    print("Importing funannotate gff...")
    fun = GFF(sys.argv[2], 'ID')
    print("Importing RATT gff...")
    ratt = GFF(sys.argv[3], 'Name')
    print("Importing manual split gff...")
    man = GFF(sys.argv[4], 'ID')
    gen = GFF(sys.argv[5], 'increment')
    rrna = GFF(sys.argv[6], 'increment')
    outd = EmptyGFF()
    print("Importing change list...")
    cl = import_change_list(man_ann_file)
    # pprint.pprint(cl)
    # sys.exit(0)
    i = 1
    print("Builing merged annotation gff...")
    genes_added = []
    for ltag, gene in fun.features.items():
        ltag_entry = check_ltag_in_change_list(cl, ltag)
        if len(ltag_entry) > 0:
            genes_to_replace = ltag_entry[:-1]
            for rgene_ltag in genes_to_replace:
                if rgene_ltag in genes_added:
                    continue
                if rgene_ltag.startswith('TOT'):
                    new_ltag = ltag.split('_')[0] + '_' + str(i).zfill(6)
                    outd.features[new_ltag] = fix_ratt(ratt.features[rgene_ltag], rgene_ltag, new_ltag)
                    print("{} added as {}.".format(rgene_ltag, new_ltag))
                    i += 1
                    #add gene from RATT (+ fix RATT attributes + reverse exon and CDS order for - strands)
                else:
                    new_ltag = ltag.split('_')[0] + '_' + str(i).zfill(6)
                    outd.features[new_ltag] = fix_attributes(man.features[rgene_ltag], rgene_ltag, new_ltag)
                    print("{} added as {}.".format(rgene_ltag, new_ltag))
                    i += 1
                genes_added.append(rgene_ltag)
        else:
            new_ltag = ltag.split('_')[0] + '_' + str(i).zfill(6)
            outd.features[new_ltag] = fix_attributes(gene, ltag, new_ltag)
            i += 1
    print("Adding non nuclear annotations...")
    for j, gene in gen.features.items():
        new_ltag = ltag.split('_')[0] + '_' + str(i).zfill(6)
        outd.features[new_ltag] = fix_geneious(gene, new_ltag)
        i += 1
    print("Adding rRNA annotations...")
    for j, gene in rrna.features.items():
        new_ltag = ltag.split('_')[0] + '_' + str(i).zfill(6)
        outd.features[new_ltag] = fix_rrna(gene, new_ltag)
        i += 1
    print("Writing merged annotation gff...")
    outd.write(sys.argv[7])
    print("Done.")
    

if __name__ == '__main__':
    main()
