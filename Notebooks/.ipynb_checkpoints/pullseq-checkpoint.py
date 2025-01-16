#!/usr/bin/env python3

from Bio import SeqIO
import argparse

def read_headers(headers_include):
    # import list of headers to pull
    headers_list = []

    with open(headers_include) as f:
        for line in f:
            headers_list.append(line.strip())

    return headers_list


def pull_sequences(fasta, headers_list):
    '''this function includes only sequences from the header list'''
    input_seq_iterator = SeqIO.parse(fasta, "fasta")
    filtered_seq_iterator = (record for record in input_seq_iterator if record.id in headers_list)
    return filtered_seq_iterator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str,
                        help='fasta file')
    parser.add_argument('-o', '--output_fasta', type=str,
                        help='output file')
    parser.add_argument('-i', '--headers_include', type=str,
                        help='headers of entries to include')

    args = parser.parse_args()
    fasta = args.fasta
    output_fasta = args.output_fasta
    headers_include = args.headers_include

    # run
    if headers_include:
        headers_list = read_headers(headers_include)
        filtered_seq_iterator = pull_sequences(fasta, headers_list)
        SeqIO.write(filtered_seq_iterator, output_fasta, "fasta")

if __name__ == '__main__':
    main()