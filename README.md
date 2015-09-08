# python_wrappers

A collection of Python wrapper classes to commonly used bioinformatic tools
at MDU

## List of tools

  1. megahit

# Using the tools

A test script is provided for each tool. As an example, if one wants to
use megahit:

    # load the library
    >>> from megahit import Megahit
    
    # create a new megahit object
    >>assem = Megahit()

    #set some defaults parameters
    >>> assemb = Megahit(k_min = 41,
                        k_max = 101,
                        k_step = 20,
                        min_count = 3,
                        min_contig_len = 500,
                        no_mercy = True)

    #assuming a pair of read files for paired-end data
    >>> reads1 = "/path/to/foo_R1.fq.gz"
    >>> reads2 = "/path/to/foo_R2.fq.gz"

    #output directory
    >>> outdir = "/path/to/foo_assemb"

    #run the assembly
    >>> assemb.assemble(seq1 = reads1, seq2 = reads2, out_dir = outdir)

    # assuming you have another set of interleaved data files
    >>> interlvd_reads = "/path/to/bar_interleave.fq.gz"
    >>> outdir = "/path/to/bar_assemb"
    >>> assemb.assemble(seq12 = interlvd_reads, out_dir = outdir)

    # assuming you have some single end reads
    >>> send_reads = "/path/to/fuzzy_reads.fq.gz"
    >>> outdir = "/path/to/fuzzy_assembly"
    >>> assemb.assemble(read = send_reads, out_dir = outdir)

    # to print current options
    >>> print(assemb)

    #to see all the options
    >>> help(assemb)

# Authors

Anders Goncalves da Silva
