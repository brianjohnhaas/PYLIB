#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import sys, os

from haasbio import *


usage = "usage: {} gtf_file.txt\n\n".format(sys.argv[0])
if len(sys.argv) < 2:
    print(usage, file=sys.stderr)
    sys.exit(1)

gtf_file = sys.argv[1]

gtf_reader = GTF_reader()
genes = gtf_reader.parse_file(gtf_file)

for gene in genes:
    collapsed_coords = gene.get_collapsed_exon_coordinates()

    gene_lend = collapsed_coords[0][0]
    gene_rend = collapsed_coords[-1][1]
    
    contig_id = gene.get_contig()
    

    collapsed_coord_text_list = list()
    for coordset in collapsed_coords:
        collapsed_coord_text_list.append("{}-{}".format(coordset[0], coordset[1]))

    print("\t".join([gene.get_id(),
                     contig_id,
                     str(gene_lend), str(gene_rend),
                     ",".join(collapsed_coord_text_list)]))





sys.exit(0)
