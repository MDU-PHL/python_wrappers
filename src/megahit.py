'''
A class that wraps around megahit to perform de novo genome assemblies
'''

class Megahit:

    def __init__(self):
        '''
        Create an instance of Megahit
        '''
        self.exists = False
        # default values used in nullarbor
        self.k_min = 41
        self.k_max = 101
        self.k_step = 20
        self.min_count = 3
        self.min_contig_len = 500
        self.no_mercy = True
        self.n_cpus = 16
        return

    def inpath(self):
        '''
        Check if megahit exists in the path, and can be return.
        '''
        return

    def create_cmd(self):
        '''
        Create the megahit command line
        '''
        return

    def run(self):
        '''
        Run megahit
        '''
        return
