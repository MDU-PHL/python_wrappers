#!/usr/bin/env python

"""
Run CPE_ongoing Nullarbor.pl requests on repeatedly updated list of MDU IDs.
Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm
YYYMMDD_HHMM: 20160513_2107 (Friday the 13th)
"""

import os
import argparse
import sys
from subprocess import call, Popen, PIPE
import shlex
import time
from collections import defaultdict

#set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description="Run this script in the working\
                                 directory (WD)on any MDU isolate set.  This\
                                 script will: i) take the MDU-IDs,\
                                 ii) find the lines pertaining to these IDs in\
                                 the WGS QC file, and iii) split the request\
                                 into jobs based on the species ID in this\
                                 table.  To capture the nullarbor.pl 'make'\
                                 commands in a file called 'screenlog.0', run\
                                 this script using 'screen' with a screenlog\
                                 (toggle the append-to-log-file function\
                                 on/off by pressing keys ctrl+a, H).  After\
                                 running, retrieve the run commands by doing\
                                 'grep nice screenlog.0'.  Alternatively,\
                                 save the terminal output to file and grep\
                                 that.  Requirements: i) In the\
                                 WD, have a RefSeqDir containing all the\
                                 .gbk refseqs, ii) inside this script, ensure\
                                 the variable REFSEQDIR is correct and the\
                                 dictionary key, value pairs in REFSEQ_DICT\
                                 are as desired.")
PARSER.add_argument("-i", "--mdu_read_ids", help="Feed this script a text file\
                    containing your MDU queries, with one MDU-ID per line.",
                    required=True)
PARSER.add_argument("-w", "--wgs_qc", help="Text file containing the WGS\
                    QC. Default '/mnt/seq/MDU/QC/mdu-wgs.csv'", 
                    default="/mnt/seq/MDU/QC/mdu-wgs.csv", required=False)
PARSER.add_argument("-p", "--prepare_directories", 
                    help="Prepare run directories by moving previous core\
                    Nullarbor outputs (e.g., 'core.tab',\
                    'core.aln', 'mlst.tab', 'report/' etc.) or 'input.tab.txt'\
                    files to new (time-stamped) files.  E.g., will move\
                    'core.tab' to 'core_YYYYMMDD-HHMMSS.tab' where\
                    YYYYMMDD-HHMMSS is today's date and time.", 
                    default="n", required=False)
PARSER.add_argument("-n", "--run_nullarbor", help="Run nullarbor.pl to\
                    generate the makefiles in the respective folders.", 
                    default="n", required=False)
ARGS = PARSER.parse_args()

WD = os.getcwd()
REFSEQDIR = "RefSeqsGenbank"
REFSEQ_DICT = {"Acinetobacter_baumannii":"Acinetobacter_baumannii_XH386.gbk",
"Aeromonas_hydrophila":"Aeromonas_hydrophila_116.gbk",
"Enterobacter_aerogenes":"Enterobacter_aerogenes_1019_EAER.gbk",
"Enterobacter_cloacae":"Enterobacter_cloacae_1000_ECLO.gbk",
"Enterobacteriaceae_bacterium_strain_FGI_57":
"Enterobacteriaceae_bacterium_strain_FGI_57.gbk",
"Erwinia_pyrifoliae":"Erwinia_pyrifoliae_DSM_12163.gbk", #### A plant pathogen
"Escherichia_coli":"Escherichia_coli_O157-H7_str._F7410.gbk",
"Klebsiella_oxytoca":"Klebsiella_oxytoca_09-7231.gbk",
"Klebsiella_pneumoniae":"Klebsiella_pneumoniae_ERS530437_PROKKA_ST258_Kp.gbk",
"Listeria_monocytogenes":"Listeria_monocytogenes_07PF0776.gbk",
"Proteus_mirabilis":"Proteus_mirabilis_1114_PMIR.gbk",
"Pseudomonas_aeruginosa":"Pseudomonas_aeruginosa_0C2E.gbk",
"Pseudomonas_putida":"Pseudomonas_putida_1A00316.gbk",
"Serratia_marcescens":"Serratia_marcescens_1145_SMAR.gbk",
"Stenotrophomonas_maltophilia":"Stenotrophomonas_maltophilia_1025_SMAL.gbk"}

NULLARBOR_MOVES = ["core.aln",
"core.full.aln",
"core.nway.tab",
"core.tab",
"core.txt",
"core.vcf",
"denovo.tab",
"distances.tab",
"isolates.txt",
"mlst.tab",
"mlst2.tab",
"ref.fa",
"ref.fa.fai",
"report",
"roary",
"tree.gif",
"tree.newick",
"tree.svg"]

TIMESTR = time.strftime("%Y%m%d-%H%M%S")

def wgs_qc_data(mdu_ids, mdu_wgs_path):
    """
    Return the mdu-wgs.csv metadata for isolates in the mdu-ids list.
    """
    isolate_basename = os.path.splitext(mdu_ids)[0]
    wgs_metadata = []
    not_found = []
    with open(mdu_ids, "r") as input_handle:
        ids = [line.rstrip().split(',')[0] for line in input_handle]
        for id in ids:
            #Capture the screen output as a string for each isolate MDU-ID.
            #Comma added after id so that IDs like xxx_M, xxx_S, xxx-1 etc., 
            #are not captured by grep
            proc = Popen(["grep", id+",", mdu_wgs_path],\
                                    stdout=PIPE)
            output = proc.stdout.read()
            if len(output) == 0:
                print id+" not found in QC file "+mdu_wgs_path+"."
                not_found.append(id)
            else:
                wgs_metadata.append(output.rstrip())
    sorted(wgs_metadata)
    #Write the WGS QC metadata to file.
    with open(isolate_basename+"_mdu-wgs.csv","w") as output_wgs_metadata:
        #Write column header from wgs QC file to metadata.
        output_wgs_metadata.write(open(mdu_wgs_path, 'r').readline())
        #Write wgs QC metadata to file.
        output_wgs_metadata.write("\n".join(wgs_metadata)+"\n")
    with open(isolate_basename+"_mdu-ids_not-found.txt", "w") as ids_not_found:
        ids_not_found.write("\n".join(not_found)+",isolate not found\n")
    print "\n====\nProcessing started at: "+TIMESTR
    print "----\nWGS QC metadata for "+ARGS.mdu_read_ids+" request:\n"
    header = open(mdu_wgs_path, 'r').readline().rstrip()
    print header
    table = "\n".join(wgs_metadata)
    print table
    with open("mdu_wgs_QC_"+isolate_basename+".csv", "w") as wgsqc_out:
        wgsqc_out.write(header+"\n"+table)
    return sorted(wgs_metadata)

def species_set_dict(wgslist):
    """
    Return the dictionary of species : MDU-IDs_list in the metadata.
    """
    #To future-proof, if new species added, need a way to work out the 
    #RefSeq list
    species = [i.split(',')[2] for i in wgslist]
    spp_set = sorted(list(set(species)))
    print "\n----\nSpecies in this ("+ARGS.mdu_read_ids+") request:\n"+\
          ", ".join(spp_set)
#     print REFSEQS
    #Initialise default dictionary: defaultdict(<type 'list'>, {}).
    spp_MDUdict = defaultdict(list)
    for i in wgslist:
        line = i.split(",")
        spp_MDUdict[line[2].replace(" ","_")].append(line[0])
    return dict(spp_MDUdict)

def prepare_dirs(spp_MDUdict):
    """
    For isolates=value by species=key in the spp_MDUdict, prepare run dirs. 
    """
    for key in spp_MDUdict:
        foldername = "sp_bin_"+key
        if os.path.exists(WD+"/"+foldername) == True:
            for i in NULLARBOR_MOVES:
                filename_parts = os.path.splitext(i)
                old = WD+"/"+foldername+"/"+i
                new = WD+"/"+foldername+"/"+filename_parts[0]+"_"+TIMESTR+\
                      filename_parts[1]
                if os.path.exists(old):  
                    call(["mv", old, new])

def mdu_reads(spp_MDUdict):
    """
    Create mdu-reads .tab files for nullarbor.pl --input input.tab
    """
    print "\n----"
    for key in spp_MDUdict:
        foldername = "sp_bin_"+key
        input_tab = WD+"/"+foldername+"/input.tab.txt"
        isos_prev_request = None
        if os.path.exists(WD+"/"+foldername+"/") == False:
            os.mkdir(foldername)
            print "Making job directory '"+foldername+"'."
        if os.path.exists(input_tab) == True:
            isos_prev_request = [i.rstrip().split("\t")[0] for i in \
                                 sorted(open(input_tab).readlines())]
        isos_this_request = sorted(spp_MDUdict[key])
        if isos_prev_request != None:
            if isos_this_request == isos_prev_request:# or listdiff==0:
                print "Zero new '"+key.replace("_"," ")+"' isolates detected "\
                "as compared to previous request."
            if len(isos_this_request) < len(isos_prev_request):
                print "Zero new '"+key.replace("_"," ")+"' isolates detected,"\
                      " and as the isolate set for this request is smaller "\
                      "than the previous set, it must be a subset of the "\
                      "previous request."
            elif len(isos_this_request) > len(isos_prev_request):
                print "New '"+key.replace("_"," ")+"'isolates detected in "\
                "this job request:"
                print ",".join(list(set(isos_this_request) -\
                          set(isos_prev_request)))
        if "y" in ARGS.prepare_directories.lower():
            if os.path.exists(input_tab):
                new_f = os.path.splitext(input_tab)[0]+"_"+TIMESTR+".txt"
                mvcmd = "mv "+input_tab+" "+new_f
                os.system(mvcmd)
                print "'Prepare_directories' requested, so old request moved "\
                      "from: "+input_tab+"\nto: "+new_f+""
        if os.path.exists(input_tab) == False:
            reads_file = input_tab
            #mdu-reads is called once on a list of isolates for every key iter.
            cmd = "mdu-reads "+' '.join(spp_MDUdict[key])+" --skip"
            args = shlex.split(cmd)
            proc = Popen(args, stdout=PIPE, stderr=PIPE)
            output = proc.stdout.read()
            read_locations = output.rstrip().split("\n")
            read_loc_to_write = [] 
            if len(read_locations) >= 3:
                for i in read_locations:
                    read_loc_to_write.append(i)
            #If there are only two isolates, create a dummy third
            #for Nullarbor compatibility
            elif len(read_locations) == 2:
                for i in read_locations:
                    read_loc_to_write.append(i)
                read_loc_to_write.append("DUMMY1_"+read_locations[0])
            #If there is only one isolate, create two dummies
            else:
                for i in read_locations:
                    read_loc_to_write.append(i)
                for j in range(1,3):
                    read_loc_to_write.append("DUMMY"+str(j)+"_"+\
                                             read_locations[0])
            sorted(read_loc_to_write)
            with open(input_tab, "w") as output_handle:
                output_handle.write("\n".join(read_loc_to_write)+"\n")
                print "Written sequence-read paths for this request ("+\
                       ARGS.mdu_read_ids+") to "+input_tab+"\n"

def execute_nullarbor(spp_MDUdict):
    """
    If requested, run Nullarbor, but only after preparing the target directory.
    Nullarbor available at https://github.com/tseemann/nullarbor
    """
    for key in spp_MDUdict:
        foldername = "sp_bin_"+key
        target_path = WD+"/"+foldername+"/"
        extant_files = []
        for i in NULLARBOR_MOVES:
            if os.path.exists(target_path+i):
                extant_files.append(i)
        if len(extant_files) > 0:
            print "\n----\nAttempted to process isolates for "+key+", but the"\
                  " following files were detected in "+target_path+":"
            print ', '.join(extant_files)
            print "Before running nullarbor.pl, prepare your target "\
                  "directories using '-p yes'.\n\n----\nExiting now.\n====\n"
            sys.exit()
        else:
            cmd = "nullarbor.pl --accurate --force --name "+key+\
                  "_CPEongoing --ref "+WD+"/"+REFSEQDIR+"/"+REFSEQ_DICT[key]+\
                  " --input "+WD+"/sp_bin_"+key+"/input.tab.txt"+\
                  " --outdir "+"sp_bin_"+key
            print "\n----\nFor '"+key.replace("_"," ")+"', will run the "\
                  "following nullarbor command:\n"+cmd
            os.system(cmd)

#Execute the functions in series.
wgs_list = wgs_qc_data(ARGS.mdu_read_ids, ARGS.wgs_qc)
spp_MDUID_dict = species_set_dict(wgs_list)
if "y" in ARGS.prepare_directories.lower():
    prepare_dirs(spp_MDUID_dict)
mdu_reads(spp_MDUID_dict)
if "y" in ARGS.run_nullarbor.lower():
    execute_nullarbor(spp_MDUID_dict)
print "\n====\nCPE_nullarbor_ongoing_requests.py run completed\n====\n====\n"
