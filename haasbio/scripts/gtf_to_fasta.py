#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import sys, os
import pickle
from haasbio import *

import shelve
import logging
import argparse


import collections

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)




def main():

    parser = argparse.ArgumentParser(description="extract sequences for genes, transcripts, and coding regions",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--gtf", dest="gtf", type=str, required=True, help="annotation in gtf format")
    
    parser.add_argument("--genome", dest="genome", type=str, required=True, help="genome fasta file")

    parser.add_argument("--seq_type", dest="seq_type", type=str, choices=('cDNA', 'CDS'), required=True, help="Feature type")

    args = parser.parse_args()


    gtf_reader = GTF_reader()
    logger.info("-Parsing gtf file: {}".format(args.gtf))
    genes_list = gtf_reader.parse_file(args.gtf)

    logger.info("-Parsing genome fasta file: {}".format(args.genome))
    fr = Fasta_reader(args.genome)

    logger.info("-Reorganizing genes according to chromosome")
    chr_to_gene_list = organize_genes_by_chromosome(genes_list)

    seq_dict = fr.get_seq_dict()

    logger.info("-Outputting sequences")
    
    for contig in chr_to_gene_list:

        if contig not in seq_dict:
            logger.error("contig: {} is not found in fasta file, skipping...".format(contig))
            continue
        
        contig_seq = seq_dict[contig]

        gene_list = chr_to_gene_list[contig]

        for gene in gene_list:

            for isoform in gene.get_isoforms():

                acc_val = ";".join([isoform.get_id(), gene.get_id()])

                seq = None

                if args.seq_type == 'cDNA':
                    seq = isoform.reconstruct_cDNA_seq(contig_seq)

                elif args.seq_type == 'CDS' and isoform.is_coding_transcript():
                    seq = isoform.reconstruct_CDS_seq(contig_seq)

                if seq is not None:
                    #print(isoform.to_string())
                    print(">{}\n{}".format(acc_val, seq))
                    


def organize_genes_by_chromosome(genes_list):

    chr_to_gene_list = collections.defaultdict(list)

    for gene in genes_list:
        contig = gene.get_contig()

        contig_gene_list = chr_to_gene_list[contig]
        contig_gene_list.append(gene)


    return chr_to_gene_list





if __name__ == "__main__":
    main()
