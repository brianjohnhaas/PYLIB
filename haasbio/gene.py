#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

    
class Seq_feature(object):

    def __init__(self, id):
        self.id = id
        self.atts = []
        self.lend = None
        self.rend = None
        self.orient = None

    def set_coords_n_orient(self, lend, rend, orient):
        self.lend = lend
        self.rend = rend
        self.orient = orient

    def set_coords_end5_end3(self, end5, end3):
        if end5 > end3:
            self.orient = '-'
            self.lend = end3
            self.rend = end5
        else:
            self.orient = '+'
            self.lend = end5
            self.rend = end3


    def get_coords(self):
        return(self.lend, self.rend)

    def get_orient(self):
        return(self.orient)

    def get_coords_n_orient(self):
        return(self.lend, self.rend, self.orient)

    def get_end5_end3(self):
        if orient == '-':
            return(self.rend, self.lend)
        else:
            return(self.lend, self.rend)

    def __repr__(self):

        ret = self.id
        if self.lend and self.rend:
            ret += " {}-{}".format(self.lend,self.rend)
        if self.orient:
            ret += "[{}]".format(orient)
        
        return(self.id)


    
class Gene(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)
        self.isoforms = []



class Isoform(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)
        self.exons = []


class Exon(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)
        self.cds_segment = None


