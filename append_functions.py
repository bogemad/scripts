import sys, csv
from Bio import SeqIO

#run as python append_function.py proteins.fasta protein_id_list target_locus_tag output


def get_protein_ids(filepath):
    pl = []
    file = open(filepath)
    csvr = csv.reader(file, delimiter='\t')
    for row in csvr:
        for val in row:
            if val.startswith(sys.argv[3]):
                pl.append(val)
    file.close()
    return pl


db = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
pl = get_protein_ids(sys.argv[2])

hyp = []
not_hyp = []
for p in pl:
    if ("hypothetical protein" in db[p].description) or ("uncharacterized protein" in db[p].description) or ("predicted protein" in db[p].description):
        hyp.append(db[p].description)
    else:
        not_hyp.append(db[p].description)

print("{} of {} ({:.2f}%) unannotated Shintoku proteins are predicted or hypothetical".format(len(hyp), len(pl), float(len(hyp))/float(len(pl))*100))

with open(sys.argv[4], 'w') as out:
    out.write("{} of {} ({:.2f}%) unannotated Shintoku proteins are predicted or hypothetical\n".format(len(hyp), len(pl), float(len(hyp))/float(len(pl))*100))
    for d in not_hyp:
        out.write('{}\n'.format(d))

