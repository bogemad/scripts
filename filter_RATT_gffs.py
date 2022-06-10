import sys, csv, gffutils, os

#run as python filter_RATT_Gffs.py manual_annotation_list RATT_gff > output_gff

tot_list = []
with open(sys.argv[1]) as csvh:
    csvr = csv.reader(csvh, delimiter=',')
    for row in csvr:
        for item in row:
            if item.startswith('TOT'):
                tot_list.append(item)

db = gffutils.create_db(sys.argv[2], dbfn='test.db', force=True, keep_order=True, merge_strategy="create_unique", sort_attribute_values=True)
db = gffutils.FeatureDB('test.db', keep_order=True)

for gene in db.features_of_type('gene', order_by=('seqid', 'start')):
    for name in gene.attributes['Name']:
        if name in tot_list:
            print(gene)
            for i in db.children(gene, featuretype='mRNA', order_by='start'):
                print(i)
            for i in db.children(gene, featuretype='tRNA', order_by='start'):
                print(i)
            for i in db.children(gene, featuretype='rRNA', order_by='start'):
                print(i)
            for i in db.children(gene, featuretype='exon', order_by='start'):
                print(i)
            for i in db.children(gene, featuretype='CDS', order_by='start'):
                print(i)

os.remove('test.db')