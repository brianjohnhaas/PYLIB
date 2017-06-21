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
        self.atts = {}
        self.contig = None
        self.lend = None
        self.rend = None
        self.strand = None

    def set_coords_n_strand(self, contig, lend, rend, strand):
        self.contig = contig
        self.lend = lend
        self.rend = rend
        self.strand = strand

    def set_coords_end5_end3(self, contig, end5, end3):
        self.contig = contig
        if end5 > end3:
            self.strand = '-'
            self.lend = end3
            self.rend = end5
        else:
            self.strand = '+'
            self.lend = end5
            self.rend = end3


    def get_contig(self):
        return self.contig

    def get_coords(self):
        return(self.lend, self.rend)

    def get_strand(self):
        return self.strand

    def get_coords_n_strand(self):
        return(self.lend, self.rend, self.strand)

    def get_end5_end3(self):
        if strand == '-':
            return(self.rend, self.lend)
        else:
            return(self.lend, self.rend)

    def add_attributes(self, atts):
        for (key, val) in atts.items():
            self.atts[key] = val
    

    def __repr__(self):

        ret = self.id
        if self.lend and self.rend:
            ret += " {}-{}".format(self.lend,self.rend)
        if self.strand:
            ret += "[{}]".format(strand)
        
        return self.id


    def copy_from_seq_feature(self, seq_feature_obj):
        self.id = seq_feature_obj.id
        self.atts = seq_feature_obj.atts
        self.contig = seq_feature_obj.contig
        self.lend = seq_feature_obj.lend
        self.rend = seq_feature_obj.rend
        self.strand = seq_feature_obj.strand


    def contains_seq_feature_coords(self, seq_feature_obj):
        if (self.contig == seq_feature_obj.contig
            and
            self.lend <= seq_feature_obj.lend
            and
            self.rend >= seq_feature_obj.rend):

            return True
        else:
            return False

    
class Gene(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)
        self.isoforms = set()

    def add_isoforms(self, transcript_obj_list):
        for transcript in transcript_obj_list:
            self.isoforms.add(transcript)
    



class Transcript(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)
        self.exons = set()


    def add_exon(self, exon_obj):
        self.exons.add(exon_obj)

    def add_exons(self, exon_obj_list):
        for exon in exon_obj_list:
            self.exons.add(exon)
    


class Exon(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)
        self.cds_segment = None

    def add_cds_segment(self, cds_segment):
        if self.has_cds_segment():
            raise RuntimeError("Error, exon {} already has a cds segment. Cannot add {}".format(self, cds_segment))

        self.cds_segment = cds_segment

    def has_cds_segment(self):
        if self.cds_segment:
            return True
        else:
            return False


class CDS_segment(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)


