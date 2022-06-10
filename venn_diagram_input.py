#!/usr/bin/env python

import sys, csv, os
from collections import defaultdict
import argparse

def write_gene_names_for_interest_groups(row):
	list = [row[0]]
	for item in row[1:]:
		if item == '':
			continue
		genes = item.split(', ')
		list += genes
	return list

def generate(orthogroups, unassigned, out):
    iso_d = defaultdict(int)

    os.makedirs(out, exist_ok = True)

    interest_groups = [
            ['Theileria_orientalis_Fish_Creek','Theileria_orientalis_Goon_Nure'], 
            ['Theileria_orientalis_Fish_Creek'],
            ['Theileria_orientalis_Goon_Nure'],
            ['Theileria_orientalis_Shintoku'],
            ]

    file_obj_l = []
    csv_l = []
    for i, group in enumerate(interest_groups):
        file_obj_l.append(open(os.path.join(out, "{}.csv".format(".".join(group))), 'w'))
        csv_l.append(csv.writer(file_obj_l[i]))

    for file in (orthogroups, unassigned):
        with open(file) as infile:
            cr = csv.reader(infile, delimiter='\t')
            for i, row in enumerate(cr):
                if i == 0:
                    column_names = row[1:]
                    continue
                iso_list = []
                for j, item in enumerate(row[1:]):
                    if item.strip() == '':
                        continue
                    iso_list.append(column_names[j])
                if len(iso_list) == 2:
                    if 'Theileria_orientalis_Fish_Creek' in iso_list and 'Theileria_orientalis_Goon_Nure' in iso_list:
                        outrow = write_gene_names_for_interest_groups(row)
                        csv_l[0].writerow(outrow)
                    # if 'Theileria_orientalis_Robertson' in iso_list and 'Theileria_orientalis_Shintoku' in iso_list:
                        # outrow = write_gene_names_for_interest_groups(row)
                        # csv_l[1].writerow(outrow)
                if len(iso_list) == 1:
                    if 'Theileria_orientalis_Fish_Creek' in iso_list:
                        outrow = write_gene_names_for_interest_groups(row)
                        csv_l[1].writerow(outrow)
                    if 'Theileria_orientalis_Goon_Nure' in iso_list:
                        outrow = write_gene_names_for_interest_groups(row)
                        csv_l[2].writerow(outrow)
                    # if 'Theileria_orientalis_Robertson' in iso_list:
                        # outrow = write_gene_names_for_interest_groups(row)
                        # csv_l[4].writerow(outrow)
                    if 'Theileria_orientalis_Shintoku' in iso_list:
                        outrow = write_gene_names_for_interest_groups(row)
                        csv_l[3].writerow(outrow)
                iso_d["&".join(iso_list)] += len(write_gene_names_for_interest_groups(row)) - 1

    with open(os.path.join(out, 'theileria.combinations.txt'), 'w') as outfile:
        for key in sorted(iso_d.keys()):
            outfile.write('{}\t{}\n'.format(key, iso_d[key]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store', dest='orthogroups', help='Orthogroups.tsv file from orthofinder')
    parser.add_argument('-u', action='store', dest='unassigned', help='Orthogroups_UnassignedGenes.tsv file from orthofinder')
    parser.add_argument('-o', action='store', dest='out', help='Output directory')
    args = parser.parse_args()
    generate(args.orthogroups, args.unassigned, args.out)
