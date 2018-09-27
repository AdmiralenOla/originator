#!/usr/bin/env python

'''
originator.py

Author: Ola Brynildsrud
Date Started: 18/09/2018

Identifying the origin of replication in bacteria (and plasmids?)
Returns replicon where dnaA starts at coordinate 1 (rotated if in the wrong way)
'''
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from pkg_resources import resource_string, resource_filename
from .__init__ import __version__
import argparse, time, os
__version__ = '0.1b'

def rev_comp(seq):
    """Returns the reverse complement of a given sequence"""
    function = {"A":"T", "T":"A", "C":"G", "G":"C"}
    output = [function[x] for x in seq[::-1]]
    return "".join(output)

def rev_comp_Bio(seq):
	"""Takes SeqIO object and reverse complements"""
	return seq.reverse_complement(id=seq.id, description= "reverse complemented")

def rotate_seq(seq,coord):
	"""Rotates the circular sequence so that coord is number 1. Should be able to handle Seq objects"""
	#Make sure coord is 0-based
	coord -= 1
	return seq[coord:] + seq[:coord]

def read_FASTA(f):
    """Read in a file in FASTA format"""
    print ("File: ", f)
    print ("Reading in FASTA...", end=' ')

    in_file = open(f,'r')
    seqDict = {}
    name = None
    for line in in_file:
        line = line.rstrip()
        if line[0]=='>':
            name = line[1:]
            seqDict[name] = ''
        else:
            seqDict[name] = seqDict[name] + line

    print ("DONE!")
    return seqDict

def find_ori_gene(seq, db):
	"""Find the location of the dnaA gene and return the coordinates"""
	blastn_cline = NcbiblastnCommandline(query=seq, db=db, outfmt=5, out='dnaA.xml')
	print("Running BLAST...")
	blastn_cline()
	blast_record = NCBIXML.parse(open('dnaA.xml'))
	hits = {'title': 'None', 'score': 0.0, 'identities': 0, 'strand':(0,0), 'gaps':0, 'positives':0, 'length':0, 'query_start':0,'frame':(0,0)}
	highscore = 0.0
	for Blast in blast_record:
		for al in Blast.alignments:
			for hsp in al.hsps:
				if hsp.score > highscore:
					hits = {'title': al.title, 'score': hsp.score, 'identities': hsp.identities, 'strand':hsp.strand, 'gaps':hsp.gaps, 'positives':hsp.positives, 'length':hsp.align_length, 'query_start':hsp.query_start,'frame':hsp.frame}
					highscore = hsp.score
					
	return hits

def read_FASTA_Bio(f):
	"""Read a FASTA file in as a SeqIO"""
	return SeqIO.read(f, "fasta")

def writeFile(fasta, outfile):
	"""Writes the rotated FASTA to file"""
	print("Writing to file %s" % outfile)
	with open(outfile,'w') as f:
		SeqIO.write(fasta,f,"fasta")

def main():
    """Main function of originator"""
    start = time.clock()
    parser = argparse.ArgumentParser(description='Originator - find oriC in bacterial genomes and rotate chromosome')
    parser.add_argument("-v", "--version", help="Installed version", action="version",
                        version="%(prog)s " + str(__version__))
    parser.add_argument("--db", help="Database of origin genes. Default: dnaA genes", default=os.path.join(resource_filename(__name__, 'db'),'dnaA.fasta'))
    parser.add_argument("fastaFile", help="Bacterial genome FASTA file")
    parser.add_argument("--outfile", help="Name of file to write originated FASTA to", default='Rotated_replicon.fasta')
    args = parser.parse_args()

    hits = find_ori_gene(args.fastaFile, args.db)
    input_fasta = read_FASTA_Bio(args.fastaFile)
    if hits['title'] == "None":
    	print("Unable to find any trace of a dnaA gene")
    else:
    	if hits['frame'] == (1,1):
    		rotated_fasta = rotate_seq(input_fasta, hits['query_start'])
    		print(hits)
    		writeFile(rotated_fasta, args.outfile)
    	elif hits['frame'] == (1,-1):
    		rotated_fasta = rotate_seq(input_fasta, hits['query_start'] + hits['length'])
    		revseq = rev_comp_Bio(rotated_fasta)
    		print(hits)
    		writeFile(revseq, args.outfile)
    	else:
    		print("Something wrong with frame. Contact Ola")

if __name__ == '__main__':
	main()