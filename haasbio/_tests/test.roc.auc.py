#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import sys, os
import haasbio
from haasbio.roc.auc import *
import logging


logging.getLogger('').setLevel(logging.DEBUG)

logger = logging.getLogger(__name__)

logger.info("testing auc calc")

roc_data =  os.path.dirname(__file__) + "/../test_data/roc_test.data"

lines = open(roc_data).readlines()

lines.pop()
lines.pop(0)

pts = []
for line in lines:
    vals = line.split("\t")
    sens = float(vals[4])
    oneMinusSpec = float(vals[6])
    pts.append((oneMinusSpec, sens))

pts.insert(0, (0,0))
pts.append((1,1))

print("AUC: {}\n".format(auc(pts)))

