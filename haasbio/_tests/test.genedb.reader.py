#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)



import sys, os
import haasbio
from haasbio.gtf.gtf_reader import *

import shelve

d = shelve.open("genes")

gene_ids = d.keys()

for gene_id in gene_ids:
    print(gene_id)

d.close()




