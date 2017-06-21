#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import sys, os
import shelve
import pickle

import haasbio
from haasbio.gtf.gtf_reader import *

usage = "usage: {} output_database_prefix\n\n".format(sys.argv[0])
if len(sys.argv) < 2:
    print(usage, file=sys.stderr)
    sys.exit(1)

db_output_prefix = sys.argv[1]


d = shelve.open(db_output_prefix)

print(pickle.loads(d['ENSG00000157911.5']).to_string())
print(pickle.loads(d['ENST00000514502.1']).to_string())

#for (feature_id, obj) in d.items():
#    print(feature_id, obj)
    
d.close()





