'''
A class that wraps around megahit to perform de novo genome assemblies
'''

import os
import subprocess

class Megahit:

    def __init__(self, cmd_path = "megahit"):
        '''
        Create an instance of Megahit
        '''
        # gives the option of automatically searching for megahit in the
        # path or for the user to provide a unique path to the megahit
        # executable
        # the class initially assumes that megahit is not in the path,
        # and that the user will provide the path. when trying to run
        # for the first time, if a command path is not provided, we search in
        # PATH and throw an exception if nothing found
        self.exists_in_path = False
        self.cmd_path = cmd_path
        # the class holds paramter values to set megahit options
        # default values set here are those used in nullarbor
        self.k_min = 41
        self.k_max = 101
        self.k_step = 20
        self.min_count = 3
        self.min_contig_len = 500
        self.no_mercy = True
        self.n_cpus = 16
        return

    def in_path(self):
        '''
        Check if megahit exists in the path.
        '''
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, "megahit")
            if os.path.isfile(exe_file):
                self.exists_in_path = True
                return
        raise RuntimeError('''
        Could NOT find an executable for megahit. Please
        make sure it is in your path.
        ''')

    def create_cmd(self):
        '''
        Create the megahit command line
        '''
        return

    def run(self, seq1 = None, seq2 = None, seq12 = None, seq = None, outdir = "./megahit_out"):
        '''
        Run megahit
        '''
        #check if there is an executable for megahit
        if self.cmd_path == "megahit" and not self.exists_in_path:
            self.in_path()
        #different flavours of commands depending on type of inputted sequences
        if (seq1 is not None and seq2 is not None):
            '''Sequences are paired-end reads'''
            for f in [seq1,seq2]:
                if not os.path.isfile(f):
                    raise RuntimeError('''
                    Cannot find one or more of the read files!
                    ''')
            cmd = [self.cmd_path, "-1", seq1, "-2", seq2, "--out-dir", outdir]
            print(cmd)
            subprocess.Popen(cmd)
            return
        elif (seq12 is not None):
            '''Paired-end interleaved sequences'''
        elif(seq is not None):
            '''Single end sequences'''
        else:
            raise ValueError('''
            Could not figure out what type of read data is
            is supposed to be used by megahit. One of the following must
            be specified:
                - seq1 and seq2 for paired-end reads in separate files
                - seq12 for interleaved paired-end reads
                - seq for single-end reads
            ''')
        return

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
