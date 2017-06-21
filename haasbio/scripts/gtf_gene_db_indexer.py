#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import sys, os
import shelve
import pickle
import haasbio
from haasbio.gtf.gtf_reader import *

usage = "usage: {} annotations.gtf output_database_prefix\n\n".format(sys.argv[0])
if len(sys.argv) < 3:
    print(usage, file=sys.stderr)
    sys.exit(1)

gtf_input_file = sys.argv[1]
db_output_prefix = sys.argv[2]



gtf_reader = GTF_reader()
genes = gtf_reader.parse_file(gtf_input_file)

d = shelve.open(db_output_prefix)

for gene in genes:
    #print(gene.to_string())
    print("-indexing {}".format(gene.get_id()))
    d[str(gene.get_id())] = pickle.dumps(gene)
    
    # index the isoforms too
    for isoform in gene.get_isoforms():
        d[str(isoform.get_id())] = pickle.dumps(isoform)

d.close()

print("Done indexing db: {}.db".format(db_output_prefix), file=sys.stderr)




