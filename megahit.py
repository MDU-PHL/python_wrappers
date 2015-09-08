'''
A class that wraps around megahit to perform de novo genome assemblies
'''

#import some BioPython sugar to create comand lines
from Bio.Application import _Option, AbstractCommandline, _Switch
import os
import subprocess

class Megahit(AbstractCommandline):
    def __init__(self, cmd = "megahit", **kwargs):
        '''
        Create an instance of Megahit
        '''
        # the class holds paramter values to set megahit options
        # default values set here are those used in nullarbor
        self.parameters = [
        #Input options
            _Option(["-1", "s1"],
            "comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>",
            equate = False),
            _Option(["-2", "s2"],
            "comma-separated list of fasta/q paired-end #2 files, paired with files in <pe1>",
            equate = False),
            _Option(["--12", "s12"],
            "comma-separated list of interleaved fasta/q paired-end files",
            equate = False),
            _Option(["--read", "read"],
            "comma-separated list of fasta/q single-end files",
            equate = False),
        #Basic assembly options
            _Option(["--min-count", "min_count"],
                "minimum multiplicity for filtering (k_min+1)-mers, default 2",
                equate = False),
            _Option(["--k-min", "k_min"],
            "minimum kmer size (<= 127), must be odd number, default 21",
            equate = False),
            _Option(["--k-max", "k_max"],
            "maximum kmer size (<= 127), must be odd number, default 99",
            equate = False),
            _Option(["--k-step","k_step"],
            "increment of kmer size of each iteration (<= 28), must be even number, default 10",
            equate = False),
            _Option(["--k-list", "k_list"],
            '''
            comma-separated list of kmer size (all must be odd, in the range 15-127, increment <= 28);
            override `--k-min', `--k-max' and `--k-step'
            ''',
            equate = False),
        #Advanced assembly options
            _Switch(["--no-mercy", "no_mercy"],
            "do not add mercy kmers"),
            _Switch(["--no-bubble", "no_bubble"],
            "do not merge bubbles"),
            _Option(["--merge-level", "merge_level"],
            "merge complex bubbles of length <= l*kmer_size and similarity >= s, default 20,0.98",
            equate = False),
            _Option(["--prune-level", "prune_level"],
            "strength of local low depth pruning (0-2), default 2"),
            _Option(["--low-local-ratio", "low_local_ration"],
            "ratio threshold to define low local coverage contigs, default 0.2",
            equate = False),
            _Option(["--max-tip-length", "max_tip_length"],
            "remove tips less than this value; default 2*k for iteration of kmer_size=k",
            equate = False),
            _Switch(["--no-local", "no_local"],
            "disable local assembly"),
            _Switch(["-kmin-1pass", "kmin_1pass"],
            "use 1pass mode to build SdBG of k_min"),
        #Preset parameters
            _Option(["--presets", "presets"],
            '''
            verride a group of parameters; possible values:
                meta            '--min-count 2 --k-list 21,41,61,81,99'             (generic metagenomes, default)
                meta-sensitive  '--min-count 2 --k-list 21,31,41,51,61,71,81,91,99' (more sensitive but slower)
                meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'       (large & complex metagenomes, like soil)
                bulk            '--min-count 3 --k-list 31,51,71,91,99 --no-mercy'  (experimental, standard bulk sequencing with >= 30x depth)
                single-cell     '--min-coun
            ''',
            equate = False),
        #Hardware options
            _Option(["--memory", "memory"],
            '''
            max memory in byte to be used in SdBG construction; default 0.9
            (if set between 0-1, fraction of the machine's total memory)
            ''',
            equate = False),
            _Option(["--mem-flag", "mem_flag"],
            '''
            SdBG builder memory mode, default 1
            0: minimum; 1: moderate; others: use all memory specified by '-m/--memory'.
            ''',
            equate = False),
            _Switch(["--use-gpu", "use_gpu"],
            "Use GPU"),
            _Option(["--gpu-mem", "gpu_mem"],
            "GPU memory in byte to be used. Default: auto detect to use up all free GPU memory.",
            equate = False),
            _Option(["--num-cpu-threads", "threads"],
            "number of CPU threads, at least 2. Default: auto detect to use all CPU threads.",
            equate = False),
        #Output options
            _Option(["--out-dir", "out_dir"],
            "output directory, default ./megahit_out",
            equate = False),
            _Option(["--out-prefix", "out_prefix"],
            "output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)",
            equate = False),
            _Option(["--min-contig-length", "min_contig_len"],
            "minimum length of contigs to output, default 200",
            equate = False),
            _Switch(["--keep-tmp-files", "keep_tmp_files"],
            "keep all temporary files"),
        #Other arguments
            _Switch(["--continue", "cont"],
            '''
            continue a MEGAHIT run from its last available check point.
            please set the output directory correctly when using this option.
            '''),
            _Switch(["--help", "help"],
            "Print usage message"),
            _Switch(["--version", "version"],
            "Print version"),
            _Switch(["--verbose","verbose"],
            "Verbose mode")
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        return

    def assemble(self, seq1 = None, seq2 = None, seq12 = None, seq = None, outdir = "./megahit_out"):
        '''
        Run megahit
        '''
        #different flavours of commands depending on type of inputted sequences
        if (seq1 is not None and seq2 is not None):
            '''Sequences are paired-end reads'''
            for f in [seq1,seq2]:
                if not os.path.isfile(f):
                    raise RuntimeError('''
                    Cannot find one or more of the read files!
                    ''')
            self.s1 = seq1
            self.s2 = seq2
        elif (seq12 is not None):
            '''Paired-end interleaved sequences'''
            self.s12 = seq12
        elif(seq is not None):
            '''Single end sequences'''
            self.read = seq
        else:
            raise ValueError('''
            Could not figure out what type of read data megahit is supposed
            to use. One of the following must be specified:
                - seq1 and seq2 for paired-end reads in separate files
                - seq12 for interleaved paired-end reads
                - seq for single-end reads
            ''')
        self.out_dir = outdir
        print(str(self))
        out, err = self()
        stats = err.split("\n")
        stats = [rec for rec in stats if rec[0:10] == "--- [STAT]"]
        stats = stats.split(",")
        print(stats)
        return
