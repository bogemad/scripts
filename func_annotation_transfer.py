#!/usr/bin/env python3

import sys, csv
from Bio import SeqIO
from collections import defaultdict

def gen_ogd(ogs):
    ogd = defaultdict(dict)
    with open(ogs) as ogsh:
        csvr = csv.reader(ogsh, delimiter='\t')
        for i, row in enumerate(csvr):
            if i == 0:
                isolates = row[1:]
            for j, header in enumerate(isolates):
                if not row[j+1].split(', ') == ['']:
                    ogd[row[0]][header] = row[j+1].split(', ')
                else:
                    ogd[row[0]][header] = []
    return ogd

def gen_product_dict(gb, rorq):
    prod_d = {}
    for chrom in gb:
        if chrom.id == 'api' or chrom.id == 'mit':
            continue
        for feat in chrom.features:
            if feat.type == 'CDS':
                if rorq == 'r':
                    prod_d[feat.qualifiers['protein_id'][0]] = feat.qualifiers['product'][0]
                elif rorq == 'q':
                    prod_d[feat.qualifiers['locus_tag'][0]] = feat.qualifiers['product'][0]
                else:
                    sys.exit("rorq is bad")
    return prod_d


def id_func_annots(ref_prod_d, q_prod_d, ogd, query, ref, qgb):
    class products: pass
    products.new_names = {}
    products.too_many_products = {}
    products.no_products = {}
    for chrom in qgb:
        if chrom.id == 'api' or chrom.id == 'mit':
            continue
        for feat in chrom.features:
            if feat.type == 'CDS':
                qltag = feat.qualifiers['locus_tag'][0]
                qproduct = feat.qualifiers['product'][0]
                og = find_query_og(ogd, qltag, query)
                if og != False:
                    ref_ids = ogd[og][ref]
                    if len(ref_ids) > 1:
                        uniq_products = [ ref_prod_d[refid].strip() for refid in ref_ids ]
                        if len(set(uniq_products)) == 1:
                            products.new_names[qltag] = uniq_products[0]
                        else:
                            products.too_many_products[qltag] = [(refid, ref_prod_d[refid]) for refid in ref_ids]
                            products.new_names[qltag] = qproduct
                    elif len(ref_ids) == 1:
                        products.new_names[qltag] = ref_prod_d[ref_ids[0]]
                    elif len(ref_ids) == 0:
                        products.new_names[qltag] = qproduct
                        products.no_products[qltag] = qproduct
                else:
                    products.new_names[qltag] = qproduct
                    products.no_products[qltag] = qproduct
    return products


def find_query_og(ogd, qltag, isolate):
    for og in list(ogd):
        if qltag+'-T1' in ogd[og][isolate]:
            return og
    return False

def generate_new_genbank_file(products, qgb):
    for chrom in qgb:
        if chrom.id == 'api' or chrom.id == 'mit':
            continue
        for feat in chrom.features:
            if feat.type == 'CDS' or feat.type == 'mRNA':
                feat.qualifiers['product'] = [products.new_names[feat.qualifiers['locus_tag'][0]]]
    return qgb

def output_product_details(products, query):
    outfile = '{}.too_many_products.csv'.format(query)
    with open(outfile, 'w') as outh:
        outw = csv.writer(outh)
        for qltag in list(products.too_many_products):
            outw.writerow([qltag, products.too_many_products[qltag][0], products.too_many_products[qltag][1]])
    outfile = '{}.no_products.csv'.format(query)
    with open(outfile, 'w') as outh:
        outw = csv.writer(outh)
        for qltag in list(products.no_products):
            outw.writerow([qltag, products.no_products[qltag]])
    
    
def main():
    ref_gb = list(SeqIO.parse(sys.argv[1], 'gb'))
    q_gbs = [list(SeqIO.parse(sys.argv[2], 'gb')), list(SeqIO.parse(sys.argv[3], 'gb'))]
    ogs = sys.argv[4]
    ref = 'Theileria_orientalis_Shintoku'
    ref_prod_d = gen_product_dict(ref_gb, 'r')
    ogd = gen_ogd(ogs)
    for i, query in enumerate(['Theileria_orientalis_Fish_Creek', 'Theileria_orientalis_Goon_Nure']):
        q_prod_d = gen_product_dict(q_gbs[i], 'q')
        products = id_func_annots(ref_prod_d, q_prod_d, ogd, query, ref, q_gbs[i])
        new_gb = generate_new_genbank_file(products, q_gbs[i])
        SeqIO.write(new_gb, '{}.shintoku_product_names.gbk'.format(query), 'gb')
        output_product_details(products, query)

if __name__ == '__main__':
    main()
