#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                                                print_function, unicode_literals)
import os, sys, re
import logging
from collections import defaultdict

from haasbio.gene import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

class GTF_reader(object):

    def __init__(self):
        self.genes = []


    def parse_file(self, gtf_filename):

        
        gene_objects = {}
        transcript_objects = {}
        transcript_to_cds_segments = defaultdict(list)

        with open(gtf_filename) as f:
            for line in f:
                if re.match("#", line):
                    continue
                if not re.search("\w", line):
                    continue
                line = line.rstrip()
                x = line.split("\t")
                if len(x) < 8:
                    raise RuntimeError("gtf record lacks 8 fields: {}".format(line))

                contig_id = x[0]
                source = x[1]
                feat_type = x[2]
                lend = x[3]
                rend = x[4]
                score = x[5]
                strand = x[6]
                frame = x[7]
                info = x[8]

                atts = self._parse_atts_from_info(info)

                atts['source'] = source
                atts['score'] = score
                atts['frame'] = frame
                
                gene_id = atts['gene_id']
                transcript_id = atts['transcript_id']

                feat_id = "::".join([feat_type, gene_id, transcript_id, contig_id, lend, rend, strand])

                seq_feat = Seq_feature(feat_id)
                seq_feat.set_coords_n_strand(contig_id, lend, rend, strand)
                seq_feat.set_attributes(atts)
                
                g = gene_objects.get(gene_id, None)
                if not g:
                    g = gene_objects[gene_id] = Gene(gene_id)

                
                # assign sequence feature attributes
                    
                if feat_type == 'gene':
                    g.copy_from_seq_feature(seq_feat)

                else:

                    # need a transcript object for all other feature operations
                    
                    t = transcript_objects.get(transcript_id, None)
                    if not t:
                        t = Transcript(transcript_id)
                        transcript_objects[transcript_id] = t
                        g.add_isoforms([t])
                    

                    if feat_type == 'transcript':
                        t.copy_from_seq_feature(seq_feat)

                    elif feat_type == 'exon':
                        e = Exon(feat_id)
                        e.copy_from_seq_feature(seq_feat)
                        t.add_exon(e)

                    elif feat_type == 'CDS':
                        # below, we'll attach these CDS features to corresponding exon features
                        transcript_to_cds_segments[transcript_id].append(seq_feat)
                    

        self._attach_CDS_to_exons(transcript_objects, transcript_to_cds_segments)



        return gene_objects.values()


    def _parse_atts_from_info(self, info_txt):
        
        atts = {}
        
        info_txt = info_txt.strip()
        
        vals = info_txt.split(";")
        for valpair in vals:
            valpair = valpair.strip()
            valpair = valpair.replace('"',"")
            if re.search(" ", valpair):
                (key,val) = valpair.split(" ", 2)
                if atts.has_key(val):
                    atts[key] += "\cA{}".format(val)
                else:
                    atts[key] = val
                    
                    
        return atts


    def _attach_CDS_to_exons(self, transcript_objects, transcript_to_cds_segments):

        for (transcript_id, cds_segment_list) in transcript_to_cds_segments.items():

            transcript_obj = transcript_objects[transcript_id]

            exons = transcript_obj.get_exons()

            for cds_segment in cds_segment_list:
                (cds_lend, cds_rend) = cds_segment.get_coords()

                mapped_to_exon = False
                for exon in exons:
                    if exon.contains_seq_feature_coords(cds_segment):
                        exon.add_cds_segment(cds_segment)
                        mapped_to_exon = True
                        break

                if not mapped_to_exon:
                    raise RuntimeError("Error, wasn't able to map cds {} to any exon of transcript {}".format(cds_segment, transcript_obj))

                
                                       
            
