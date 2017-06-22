#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging
import Bio.Seq

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



class bcolors:
    # borrowed from: https://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
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
    

    def get_attributes(self):
        return self.atts

    def set_attributes(self, atts):
        self.atts = atts

    def get_id(self):
        return self.id

    def __repr__(self):

        ret = self.id
        if self.lend and self.rend:
            ret += " {}-{}".format(self.lend,self.rend)
        if self.strand:
            ret += "[{}]".format(self.strand)
        
        return self.id


    def get_coords_string(self):
        return("{} {}-{} [{}]".format(self.contig, self.lend, self.rend, self.strand))

    def get_atts_string(self):

        key_val_pairs = list()
        for (key,val) in self.atts.items():
            key_val_pairs.append("{}:{}".format(key,val))

        return("; ".join(key_val_pairs))
    

    def copy_from_seq_feature(self, seq_feature_obj):
        # retain original id
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
    
    def get_isoforms(self):
        return self.isoforms


    def to_string(self):

        ret = str(bcolors.OKBLUE + "Gene:" + 
                  " {} {}\n".format(self.get_id(), self.get_coords_string()) + bcolors.ENDC +
                  "\tAtts: {}\n".format(self.get_atts_string()) )
                  
        for isoform in self.get_isoforms():
            isoform_str = isoform.to_string()
            for line in isoform_str.split("\n"):
                ret += "\t" + line + "\n"

        return ret


    def get_collapsed_exon_coordinates(self):

        coordsets = list()
        for isoform in self.get_isoforms():
            for exon in isoform.get_exons():
                (lend, rend) = exon.get_coords()
                coordsets.append( [lend, rend] )

        coordsets = sorted(coordsets, key=lambda x:x[0])

        collapsed = list()
        collapsed.append(coordsets.pop(0))

        while coordsets:
            prev_rend = collapsed[-1][1]
            next_coordset = coordsets.pop(0)
            if next_coordset[0] <= prev_rend:
                if next_coordset[1] > prev_rend:
                    # extend existing range
                    collapsed[-1][1] = next_coordset[1]
            else:
                # not overlapping, add new segment
                collapsed.append(next_coordset)

        return collapsed


    def refine_gene(self):

        coordsets = list()
        contig_id = None
        strand = None
        for isoform in self.get_isoforms():
            isoform.refine_isoform()
            (lend, rend) = isoform.get_coords()
            coordsets.append( (lend, rend) )
            contig_id = isoform.get_contig()
            strand = isoform.get_strand()

        coordsets = sorted(coordsets, key=lambda x:x[0])

        gene_lend = coordsets[0][0]
        gene_rend = coordsets[-1][1]

        self.contig = contig_id
        self.lend = gene_lend
        self.rend = gene_rend
        self.strand = strand
        

class Transcript(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)
        self.exons = set()


    def get_exons(self):
        exons = self.exons
        exons = sorted(exons, key=lambda e:e.lend)
        return exons

    def add_exon(self, exon_obj):
        self.exons.add(exon_obj)

    def add_exons(self, exon_obj_list):
        for exon in exon_obj_list:
            self.exons.add(exon)
    
    def to_string(self):
        ret = str(bcolors.OKGREEN + "Transcript:" + 
                  " {} {}\n".format(self.get_id(), self.get_coords_string()) + bcolors.ENDC +
                  "\tAtts: {}\n".format(self.get_atts_string()) )

        exons = self.get_exons()
        for exon in exons:
            exon_coord_str = exon.get_coords_string()
            ret += "\t\tExon: {} {}".format(exon.get_id(), exon.get_coords_string()) + "\n"
            if exon.has_cds_segment():
                cds_segment = exon.get_cds_segment()
                ret += "\t\t\tCDS: {} {}".format(cds_segment.get_id(), cds_segment.get_coords_string()) + "\n"
            
                
        return ret


    def is_coding_transcript(self):
        for exon in self.get_exons():
            if exon.has_cds_segment():
                return True

        return False


    def reconstruct_cDNA_seq(self, contig_seq):

        seq = ""

        for exon in self.get_exons():
            (lend, rend) = exon.get_coords()
            exon_seq = contig_seq[lend-1:rend]
            seq += exon_seq


        if seq and self.get_strand() == '-':
            bioseq = Bio.Seq.MutableSeq(seq)
            bioseq.reverse_complement()
            seq = str(bioseq)
            
            
        return seq


    def reconstruct_CDS_seq(self, contig_seq):

        if not self.is_coding_transcript():
            return None
        
        seq = ""

        for exon in self.get_exons():
            if exon.has_cds_segment():
                cds_segment = exon.get_cds_segment()
                (lend, rend) = cds_segment.get_coords()
                cds_seq = contig_seq[lend-1:rend]
                seq += cds_seq
        

        if seq and self.get_strand() == '-':
            bioseq = Bio.Seq.MutableSeq(seq)
            bioseq.reverse_complement()
            seq = str(bioseq)

        return seq
    
    def refine_isoform(self):

        coordsets = list()
        contig_id = None
        strand = None
        for exon in self.get_exons():
            (lend, rend) = exon.get_coords()
            coordsets.append( (lend, rend) )
            contig_id = exon.get_contig()
            strand = exon.get_strand()

        coordsets = sorted(coordsets, key=lambda x:x[0])

        isoform_lend = coordsets[0][0]
        isoform_rend = coordsets[-1][1]

        self.contig = contig_id
        self.lend = isoform_lend
        self.rend = isoform_rend
        self.strand = strand


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


    def get_cds_segment(self):
        return self.cds_segment
    

class CDS_segment(Seq_feature):

    def __init__(self, id):
        super(self.__class__, self).__init__(id)


