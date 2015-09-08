#!/usr/local/env python
'''
A script to test the functionality of megahit.py wrapper
'''

from megahit import Megahit

# some test data on server
reads1 = "/mnt/seq/MDU/READS/2015-12600/2015-12600_R1.fastq.gz"
reads2 = "/mnt/seq/MDU/READS/2015-12600/2015-12600_R2.fastq.gz"

#create a new Megahit instance
# set it up with some nullarbor defaults
assemb = Megahit(k_min = 41,
                    k_max = 101,
                    k_step = 20,
                    min_count = 3,
                    min_contig_len = 500,
                    no_mercy = True)

# print the options
print "Megahit will be run with these options"
print assemb
print "\n"

# run assembly for read set above
print "Running Megahit"
assemb.assemble(seq1 = reads1, seq2 = reads2)
print "\n"

print "Megahit has successfully run!"
