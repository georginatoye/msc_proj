#!/bin/bash


# ---------------------------------------------
# LOAD JSON OBJECTS & EXTRACT CDR-H3 FIELDS 
# ---------------------------------------------

# This script takes as input a folder containing compressed JSON (.json.gz)
# files downloaded from the Observed Antibody Space (OAS) repository. First,
# the data for analysis is copied into folder /d/projects/u/tg001/JSON/data
# ready to run the script.
#
# After determining filepaths of files held by the defined directory, this
# script ungzips any compressed JSON (.json.gz) files. Once decompressed, each
# JSON file represents one sample data unit i.e. one set of pre-processed
# Ig-seq sequences, and its associated sample metadata.
#
# The script can also be used to extract the same data from regular JSON 
# files, and could be adapted for other databases using .json/.json.gz format.
#
# The aim was to extract all sequence_ids & associated CDR-H3 sequences from
# each sample. This was done by loading each file of JSON objects. The
# first JSON object contains the sample metadata, while subsequent JSON
# objects hold the sequencing data.
#
# Once the script has run, data from each JSON file will have been streamed into
# the directory /d/projects/u/tg001/JSON/programs/5_Results/ as separate folders:

# (1) Paths: 	A single file listing all filepaths
# (2) Files:	Numbered files with each file's contents
# (3) Parsed:	Numbered files containing parsed sequencing data from each file
# (4) Meta:	  A single file listing metadata parsed from each file
# (5) CDRH3s:	Numbered files each containing a Python dictionary of the sample
#                CDRH3s i.e. (sequence_id {key}/CDRH3 AAs {values})
# (6) FASTA:	Numbered FASTA files of CDR-H3 sequences from each file, stored
#               on the server, ready for entry into a MEME Suite application
# -----------------------------------------------------------------------------



import os
import gzip
import json
import pprint
import re
import random 

seqset = 0



# Define function to find filepaths of files in a directory
def list_filepaths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))

# Define function to open and save contents of JSON files
def ungzip_file(src):

    # Compressed JSON files are ungzipped
    if f.endswith('.gz'):
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as u:
            print('\nFile_#', seqset, ' ungzipped JSON file '\
              'contents:\n', f, sep = '', file=u)
        for line in gzip.open(src,'rb'):
            if len(line) > 3:           
                with open('/d/projects/u/tg001/JSON/programs/'\
                  '5_Results/Files/file_#' + str(seqset) + '.txt','a') as u:
                    print('\n', line, '\n------------', file=u)
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as u:
            print('\nend_File_#', seqset, \
              '\n\n-------------\n', sep='', file=u)
           
    # Other JSON files are opened and contents saved
    if not f.endswith('.gz') and f.endswith('.json'):
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as n:
            print('\nfile_#', seqset, ' JSON file '\
              'contents:\n', f, sep = '', file=n)
        for line in open(src,'rb'):
            if len(line) > 3:
                with open('/d/projects/u/tg001/JSON/programs/'\
                  '5_Results/Files/file_#' + str(seqset) + '.txt','a') as n:
                    print('\n', line, '\n------------', file=n)
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as n:
            print('\nend_file_#', seqset, \
              '\n\n-------------\n', sep='', file=n)
               
# Define function to parse JSON files
def parse_file(src):
    global seq_data
    
    # Compressed JSON files are ungzipped before parsing
    if f.endswith('.gz'):
        meta_line = True
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Meta/meta.txt', 'a') as d:
            print('\nFile_#', seqset, ' Sample metadata: \n'\
              , f, sep='', file=d)
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Parsed/parsed_#' + str(seqset) + '.txt', 'a') as s:
            print('File_#', seqset, ' Parsed_data: \n'\
              , f, sep='', file=s)

        # Parse first JSON object in each file - metadata
        for line in gzip.open(src,'rb'):
            if meta_line == True and len(line) > 3:      
                metadata = json.loads(line)
                meta_line = False
                with open('/d/projects/u/tg001/JSON/programs/'\
                  '5_Results/Meta/meta.txt', 'a') as m:
                    print("\nMetadata = \n", end = '', file=m)
                    pprint.pprint(metadata, stream=m)
                    continue

            # Parse subsequent JSON objects - sequences & data
            if len(line) > 3: 
                seq_data = json.loads(line)
                pprint.pprint(seq_data)
                with open('/d/projects/u/tg001/JSON/programs/5_Results'\
                  '/Parsed/parsed_#' + str(seqset) + '.txt', 'a') as b:
                    print("\nSequence_data = \n", end = '', file=b)
                    pprint.pprint(seq_data, stream=b)
                 
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Meta/meta.txt', 'a') as i:                
            print('\nend_file_#',seqset,'\n\n-------------\n',sep='',file=i)
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Parsed/parsed_#' + str(seqset) + '.txt', 'a') as i:                
            print('\nend_file_#',seqset,'\n\n-------------\n',sep='',file=i)

    # Other JSON files
    if not f.endswith('.gz'):
        meta_line = True
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Meta/meta.txt', 'a') as d:
            print('\nFile_#', seqset, ' Sample metadata: \n'\
              , f, sep='', file=d)
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Parsed/parsed_#' + str(seqset) + '.txt', 'a') as s:
            print('File_#', seqset, ' Parsed_data: \n'\
              , f, sep='', file=s)

        # Parse first JSON object in each file - metadata
        for line in open(src,'rb'):
            if meta_line == True and len(line) > 3:      
                metadata = json.loads(line)
                meta_line = False
                with open('/d/projects/u/tg001/JSON/programs/'\
                  '5_Results/Meta/meta.txt', 'a') as m:
                    print("\nMetadata = \n", end = '', file=m)
                    pprint.pprint(metadata,stream=m)
                    continue
                  
            # Parse subsequent JSON objects - sequences & data
            if len(line) > 3: 
                seq_data = json.loads(line)
                with open('/d/projects/u/tg001/JSON/programs/'\
                  '5_Results/Parsed/parsed_#'+str(seqset) + '.txt', 'a') as b:
                    print("\nSequence_data = \n", end = '', file=b)
                    pprint.pprint(seq_data, stream=b)
                    
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Meta/meta.txt', 'a') as i:                
            print('\nend_file_#',seqset,'\n\n-------------\n',sep='',file=i)
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Parsed/parsed_#' + str(seqset) + '.txt', 'a') as i:                      
            print('\nend_file_#', seqset,'\n\n-------------\n', sep='',file=i) 

# Define function to extract cdrh3 sequences and save them in fasta format
def cdrh3_regex(src):

    # Compressed JSON files
    if f.endswith('.gz'):

        # Use of regex to select sample id and cdrh3 sequence
        with open('/d/projects/u/tg001/JSON/programs/'\
         '5_Results/Files/file_#' + str(seqset) +'.txt', 'rt') as file:
         
            for line in file:           
                cdrh3_flag = re.compile(r'((\",[\s]+\"cdr3\":[\s]+\")([A-Z]+)'\
                   '\",[\s]+\"original_name\":[\s]+\"(.+?)\",[\s]+\"errors)')
                cdrh3_found = re.search(cdrh3_flag, line)
                if cdrh3_found != None:
                    samp_name = cdrh3_found.group(4)
                    cdrh3_seq = cdrh3_found.group(3)

                    # cdrh3 sequences are recorded consecutively
                    samp_names.append(samp_name)  
                    cdrh3_seqs.append(cdrh3_seq)
      
                    # short cdrh3 sequences are removed (meme minimum 8 aa)  
                    if len(cdrh3_seq) > 7:
                        meme_samp = samp_name
                        meme_cdrh3 = cdrh3_seq                     
                        meme_samps.append(meme_samp)  
                        meme_cdrh3s.append(meme_cdrh3)

            # Small sequence sets are removed (minimum of five sequences)
            if len(meme_cdrh3s) > 4 and len(meme_cdrh3s) <= 40000:
                length = len(meme_cdrh3s)
                print('\nNumber of sequences = ', length)         

                # Cumulative list of samples for meme
                meme_sub1 = dict(zip(meme_samps, meme_cdrh3s))
                with open('/d/projects/u/tg001/JSON/programs/'\
                 '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                   print('\nSUBMIT TO MEME:\n\nDictionary (n =',length,')\n',\
                    meme_sub1, '\n\nAs FASTA file.\n', file = l)
 
                # FASTA file prepared for submission to meme
                print('\nAs FASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub1.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/JSON/programs/5_Results'\
                     '/FASTA/fasta_json#'+str(seqset)+'.fa','a') as s:
                      print('\n','>',meme_samp,'\n',meme_cdrh3,sep='',file=s)
         
            elif len(meme_cdrh3s) > 40000:
                length = len(meme_cdrh3s)
                print('\nSince number of sequences was', length, '\nthe '\
                 'dictionary is reduced: 40000 randomly selected sequences')

                # Cumulative list of samples for meme
                meme_dict = dict(zip(meme_samps, meme_cdrh3s))
                meme_sub2 = dict(random.sample(meme_dict.items(), 40000))
                with open('/d/projects/u/tg001/JSON/programs/'\
                  '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nSUBMIT TO MEME:\n\nDictionary (n = 40000)\n',\
                    meme_sub2, '\n\nAs FASTA file.\n', file = l)

                # FASTA file prepared for submission to meme
                print('\nAs FASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub2.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/JSON/programs/5_Results'\
                     '/FASTA/fasta_json#'+str(seqset)+'.fa','a') as s:
                       print('\n','>',meme_samp,'\n',meme_cdrh3,sep='',file=s)
            else:
                print('\nToo few sequences for submission to MEME')
                with open('/d/projects/u/tg001/JSON/programs/'\
                   '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nToo few sequences for submission to MEME\n', \
                       file = l)
                          
    # Other JSON files
    if not f.endswith('.gz'):
       
        # Use of regex to select sample id and cdrh3 sequence
        with open('/d/projects/u/tg001/JSON/programs/'\
         '5_Results/Files/file_#' + str(seqset) +'.txt', 'rt') as file:
           
            for line in file:           
                cdrh3_flag = re.compile(r'((\",[\s]+\"cdr3\":[\s]+\")([A-Z]+)'\
                     '\",[\s]+\"original_name\":[\s]+\"(.+?)\",[\s]+\"errors)')
                cdrh3_found = re.search(cdrh3_flag, line)
                if cdrh3_found != None:
                    samp_name = cdrh3_found.group(4)
                    cdrh3_seq = cdrh3_found.group(3)

                    # cdrh3 sequences are recorded consecutively
                    samp_names.append(samp_name)  
                    cdrh3_seqs.append(cdrh3_seq)
      
                    # short cdrh3 sequences are removed (meme minimum 8 aa)  
                    if len(cdrh3_seq) > 7:
                        meme_samp = samp_name
                        meme_cdrh3 = cdrh3_seq                     
                        meme_samps.append(meme_samp)  
                        meme_cdrh3s.append(meme_cdrh3)

            # Small sequence sets are removed (minimum of five sequences)
            if len(meme_cdrh3s) > 4 and len(meme_cdrh3s) <= 40000:
                length = len(meme_cdrh3s)
                print('\nNumber of sequences = ', length)         

                # Cumulative list of samples for meme
                meme_sub1 = dict(zip(meme_samps, meme_cdrh3s))
                with open('/d/projects/u/tg001/JSON/programs/'\
                 '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                   print('\nSUBMIT TO MEME:\n\nDictionary (n =',length,')\n',\
                    meme_sub1, '\n\nAs FASTA file.\n', file = l)

                # FASTA file prepared for submission to meme
                print('\nAs FASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub1.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/JSON/programs/5_Results'\
                     '/FASTA/fasta_json#'+str(seqset)+'.fa','a') as s:
                      print('\n','>',meme_samp,'\n',meme_cdrh3,sep='',file=s)
         
            elif len(meme_cdrh3s) > 40000:
                length = len(meme_cdrh3s)
                print('\nSince number of sequences was', length, '\nthe '\
                  'dictionary is reduced: 40000 randomly selected sequences')

                # Cumulative list of samples for meme
                meme_dict = dict(zip(meme_samps, meme_cdrh3s))
                meme_sub2 = dict(random.sample(meme_dict.items(), 40000))
                with open('/d/projects/u/tg001/JSON/programs/'\
                  '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nSUBMIT TO MEME:\n\nDictionary (n = 40000)\n',\
                    meme_sub2, '\n\nAs FASTA file.\n', file = l)
 
                # FASTA file prepared for submission to meme
                print('\nAs FASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub2.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/JSON/programs/5_Results'\
                     '/FASTA/fasta_json#'+str(seqset)+'.fa','a') as s:
                      print('\n','>',meme_samp,'\n',meme_cdrh3,sep='',file=s)
            else:
                print('\nToo few sequences for submission to MEME')
                with open('/d/projects/u/tg001/JSON/programs/'\
                   '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nToo few sequences for submission to MEME\n', \
                      file = l)
                    
# Directory is defined: data file is placed here before running the script
if __name__ == '__main__':
    directory = '/d/projects/u/tg001/JSON/data'

    # Files are labelled and functions are called 
    for f in list_filepaths(directory):
        samp_names = []
        cdrh3_seqs = []
        meme_samps = []
        meme_cdrh3s = []
        meme_sub1 = {}
        meme_sub2 = {}
        meme_samp = ''
        meme_cdrh3 = ''
        seqset += 1
        print('\nFile_ref_#', seqset, '\n', f, sep='')
        with open('/d/projects/u/tg001/JSON/programs/'\
          '5_Results/Paths/paths.txt', 'a') as p:
            print('\nFile_#', seqset, '\n', f, sep='', file=p)
        ungzip_file(f)
        parse_file(f)
        cdrh3_regex(f)
        print('\n-------------\n\nend_file_#', seqset, '\n', \
              '\n-------------\n', sep='')

