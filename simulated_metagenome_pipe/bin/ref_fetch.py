#!/usr/bin/env python3

from Bio import Entrez
import sys


def fetch_fasta(accession, email):
    """Retrieve FASTA record from NCBI RefSeq database."""

    Entrez.email = email

    # Use the esearch utility to find the GI number associated with the accession number
    handle = Entrez.esearch(db="nucleotide", term=accession)
    record = Entrez.read(handle)
    gi = record["IdList"][0]

    # Use the efetch utility to retrieve the sequence record in FASTA format
    handle = Entrez.efetch(db="nucleotide", id=gi, rettype="fasta", retmode="text")

    return handle


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", type=str, required=True)
    parser.add_argument("--email", type=str, required=True)
    args = parser.parse_args()

    fasta_handle = fetch_fasta(args.accession, args.email)

    print(fasta_handle, file=sys.stdout)


if __name__ == "__main__":
    main()
