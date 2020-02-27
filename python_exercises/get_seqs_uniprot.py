#!/usr/bin/python3

from Bio import ExPASy
from Bio import SwissProt


# input data
my_prots = ['KCNE1_HUMAN', 'KCNK2_HUMAN', 'KCMB1_HUMAN']

# get data from Uniprot
records = []
for prot in my_prots:
    with ExPASy.get_sprot_raw(prot) as handle:
        record = SwissProt.read(handle)
        records.append(record)  

# Print info and write sequence as fasta
with open('uniprot.fasta', 'w') as fasta_out:
    for rec in records:
        print('{0}\t{1}'.format(rec.entry_name, ' '.join(rec.accessions)))
        line_out = '>{}\n{}\n'.format(rec.entry_name, rec.sequence)
        fasta_out.write(line_out)

