#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import sys, os
import haasbio
from haasbio.gtf.gtf_reader import *

import shelve

gtf_file = os.path.dirname(__file__) + "/haasbio/test_data/small.gtf"

gtf_reader = GTF_reader()
genes = gtf_reader.parse_file(gtf_file)

d = shelve.open("genes")

for gene in genes:
    print(gene.to_string())
    d[gene.get_id()] = gene

d.close()




