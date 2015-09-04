'''
A class that wraps around megahit to perform de novo genome assemblies
'''

class Megahit:

    def __init__(self, cmd_path = None):
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
        Check if megahit exists in the path, and can be return.
        '''
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, "megahit")
            if os.path.isfile(exe_file):
                self.exists_in_path = True
            else:
                raise RuntimeError('''
                Could find an executable for megahit. Please
                make sure it is in your path.
                ''')
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
