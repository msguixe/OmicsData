#!/usr/bin/env python3

from Bio import Entrez

Entrez.email = "your@mail.here"

with Entrez.esearch(
    db="nucleotide", term="5-hydroxytryptamine receptor"
) as search_handle:
    search_results = Entrez.read(search_handle)

list_of_ids = search_results["IdList"]

with Entrez.efetch(
    db="nucleotide", id=list_of_ids, rettype="gb", retmode="text"
) as fetch_handle:
    data = fetch_handle.read()

with open("5ht.gb", "w") as fasta_out:
    fasta_out.write(data)
