#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


class Fasta_reader(object):

    def __init__(self, fastaFile):
        self.filename = fastaFile
        self.fh = open(fastaFile, 'r')
        self.curr_line = "\n";

    def next_entry(self):
        while (self.curr_line and self.curr_line.find(">") != 0):
            self.curr_line = self.fh.readline()
        if not self.curr_line:
            # end of file reached
            return(None)
        header = self.curr_line
        header = header.rstrip()
        header = header[1:] # remove carat
        accession = header.split()[0]

        sequenceLines = []
        self.curr_line = self.fh.readline()
        while (self.curr_line and self.curr_line.find(">") != 0):
            sequenceLines.append(self.curr_line)
            self.curr_line = self.fh.readline()

        sequence = "".join(sequenceLines)

        fastaSequenceObject = Fasta_sequence(accession, header, sequence)
        return(fastaSequenceObject)


    def iter_entries(self):

        while True:
            fastaSeqObj = self.next_entry()
            if fastaSeqObj is not None:
                yield fastaSeqObj
            else:
                break
            


    def get_seq_dict(self):

        seq_dict = dict()
        for fasta_entry in self.iter_entries():
            acc = fasta_entry.get_accession()
            seq = fasta_entry.get_sequence()

            seq_dict[acc] = seq

        return seq_dict


    def close(self):
        self.fh.close()
        

class Fasta_sequence(object):

    def __init__(self, accession, header, sequence):
        self.accession = accession
        self.header = header
        self.sequence = self._getRawSequence(sequence)


    def get_accession(self):
        return self.accession

    def get_header(self):
        return self.header

    def get_sequence(self):
        return self.sequence

    def _getRawSequence(self, seq):
        seq = re.sub(r"\s", "", seq) # remove all whitespace
        return(seq)
    
    


if __name__ == '__main__':
    usage = "\n\nusage: " + sys.argv[0] + " fastaFile\n\n";
    if (len(sys.argv) < 2):
        raise Exception(usage)

    fasta_filename = sys.argv[1]
    fasta_reader = Fasta_reader(fasta_filename)
    for fastaSeqObj in fasta_reader.iter_entries():
        print("Accession: {}".format(fastaSeqObj.get_accession()))
        print("Header: {}".format(fastaSeqObj.get_header()))
        print("Sequence: {}".format(fastaSeqObj.get_sequence()))
        print("\n\n")
    
    sys.exit(0)

