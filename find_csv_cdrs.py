#!/bin/bash

# ----------------------------------------------
# IMPORT CSV FILES & EXTRACT CDR-H3 FIELDS
# ----------------------------------------------

# This script takes as input a folder containing compressed CSV (.csv.gz)
# files downloaded from the Observed Antibody Space (OAS) repository. First,
# the data for analysis is copied into folder /d/projects/u/tg001/CSV/data
# ready to run the script.
#
# After determining filepaths of files held by the defined directory, this
# script ungzips any compressed CSV (.csv.gz) files. Once decompressed, each
# CSV file represents one sample data unit i.e. one set of pre-processed
# Ig-seq reads, and its associated sample metadata.
# 
# The script can also be used to extract the same data from regular CSV
# files, and could be adapted for other databases using .csv/.csv.gz format.
#
# The aim was to extract all read sequence IDs & associated CDR-H3 sequences
# from each sample. This was done by loading the data from each sample unit
# into a PANDAS dataframe. Sample metadata was also extracted from OAS files.
#
# Once the script has run, data from each CSV file will have been streamed 
# into the directory /d/projects/u/tg001/CSV/programs/5_Results as separate 
# folders:

# (1) Paths: 	A single file listing all filepaths
# (2) Files:	Numbered files with each file's contents
# (3) Fields:	Numbered files containing certain fields from each dataframe
# (4) Meta:	A single file listing metadata parsed from each file
# (5) CDRH3s:	Numbered files each containing a dictionary of the sample 
#               CDR-H3s i.e.(sample name {key}/CDR-H3 sequence {values})
# (6) FASTA:	Numbered FASTA files of CDR-H3 sequences from each file, stored
#               on the server, ready for entry into a MEME Suite application
# -----------------------------------------------------------------------------



import os
import gzip
import re
import random
   
import pandas as pd
pd.set_option('display.max_rows', None)

seqset = 0



# Define function to find filepaths of files in a directory
def list_filepaths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))


# Define function to open and save contents of CSV files
def ungzip_file(src):

    # Compressed CSV files are first ungzipped
    if f.endswith('.gz'):
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as u:
            print('\nFile_ref_#', seqset, ' ungzipped CSV file '\
              'contents:\n', f, sep='', file=u)
        for line in gzip.open(src,'rb'):
            with open('/d/projects/u/tg001/CSV/programs/'\
              '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as u:
                print('\n', line, '\n------------', file=u)
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as u:
            print('\nend_file_#', seqset, \
              '\n\n-------------\n', sep='', file=u)
            
    # Other CSV files are then opened and contents saved
    if not f.endswith('.gz'):
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as n:
            print('\nFile_ref_#', seqset, ' CSV file '\
              'contents\n', f, sep='', file=n)
        for line in open(src,'rb'):
            with open('/d/projects/u/tg001/CSV/programs/'\
              '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as n:
                print('\n', line, '\n------------', file=n)
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Files/file_#' + str(seqset) + '.txt', 'a') as n:
            print('\nend_file_#', seqset, \
              '\n\n-------------', sep='', file=n)

# Define function to save metadata from compressed CSV files
def metadata(src):
    if f.endswith('.gz'):
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Meta/meta_file.txt', 'a') as m:
            print('File_ref_#',seqset,'Patient metadata:\n',f,file=m)

        with gzip.open(f) as g:  
            df = pd.read_csv(g)
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Meta/meta_file.txt', 'a') as m:
            print('\n', df.head(0).to_string(),'\n-------------\n\n'\
              'end_metadata_file_#',seqset,'\n\n-------------', sep='',file=m)           
 
# Define function to create pandas dataframes from file contents
def select_fields(src):

    # Compressed CSV files
    if f.endswith('.gz'):
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Fields/fields_file_#' + str(seqset) + '.txt', 'a') as d:
            print('File_ref_#', seqset, ' Pandas df - selected columns: \n'\
              'Sample name & CDR-H3 sequence\n', f, sep='', file=d)
            
        # Sequence ID and CDR-H3 columns are selected. Metadata is ignored.
        with gzip.open(f) as g:  
            df = pd.read_csv(g, skiprows=1)
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Fields/fields_file_#' + str(seqset) + '.txt', 'a') as d:
            print('\n', df[['sequence_id','cdr3_aa']],'\n------------',\
              '\nend_file_#', seqset, '\n-------------', sep='', file=d)

    # Other CSV files
    if not f.endswith('.gz'):      
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Fields/fields_file_#' + str(seqset) + '.txt', 'a') as d:    
            print('File_ref_#', seqset, ' Pandas df - selected columns: \n'\
               'Sample name & CDR-H3 sequence\n', f, sep='', file=d)
            
        # Name and CDR-H3 columns are selected. No metadata in these files.
        df = pd.read_csv(f)
        print('\n', df.columns.values, '\n')
        print('\n', df[['Name', 'CDRH3']],'\n', \
          '\n-------------', sep='')
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Fields/fields_file_#' + str(seqset) + '.txt', 'a') as d:
            print('\n', df[['Name', 'CDRH3']], '\n-------------',\
              '\nend_file_#', seqset, '\n-------------', sep='', file=d)


# Define function to extract cdrh3 sequences and save them in fasta format
def cdrh3_regex(src):

    # Compressed CSV files
    if f.endswith('.gz'):
        
        # Use of regex to select sample id and cdrh3 sequence
        with open('/d/projects/u/tg001/CSV/programs/'\
         '5_Results/Fields/fields_file_#' + str(seqset) +'.txt', 'rt') as file:
    
            for line in file:
                cdrh3_flag = re.compile(r"^((\d+)(\s+)"\
                  "([^\s]+)(\s+)([A-Z]+)(\s+))", re.M)
                cdrh3_found = re.search(cdrh3_flag, line)
                if cdrh3_found != None:        
                    samp_name = cdrh3_found.group(4)
                    cdrh3_seq = cdrh3_found.group(6)

                    # cdrh3 sequences for each sample are recorded
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
                with open('/d/projects/u/tg001/CSV/programs/'\
                  '5_Results/CDRH3s/cdrh3_#'+ str(seqset) + '.txt','a') as l:
                    print('\nSUBMIT TO MEME:\n\nDictionary (n =', length,')',\
                     '\n', meme_sub1, '\n\nAs FASTA file.\n', file = l)
 
                # FASTA file prepared for submission to meme
                print('\nFASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub1.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/CSV/programs/5_Results'\
                     '/FASTA/fasta_csv#'+str(seqset)+'.fa','a') as s:
                       print('\n','>',meme_samp,'\n',meme_cdrh3,sep='',file=s)
         
            elif len(meme_cdrh3s) > 40000:
                length = len(meme_cdrh3s)
                print('\nSince number of sequences was', length, '\nthe '\
                  'dictionary is reduced: 40000 randomly selected sequences')

                # Cumulative list of samples for meme
                meme_dict = dict(zip(meme_samps, meme_cdrh3s))
                meme_sub2 = dict(random.sample(meme_dict.items(), 40000))
                with open('/d/projects/u/tg001/CSV/programs/'\
                  '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nSUBMIT TO MEME:\n\nDictionary (n = 40000)\n',\
                    meme_sub2, '\n\nAs FASTA file.\n', file = l)
                    
                # FASTA file prepared for submission to meme
                print('\nFASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub2.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/CSV/programs/5_Results'\
                     '/FASTA/fasta_csv#'+str(seqset)+'.fa','a') as s:
                      print('\n','>', meme_samp,'\n',meme_cdrh3,sep='',file=s)

            else:
                print('\nToo few sequences for submission to MEME')
                with open('/d/projects/u/tg001/CSV/programs/'\
                   '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nToo few sequences for submission to MEME\n', \
                       file = l)
                   
    # Other CSV files
    if not f.endswith('.gz'):
        
        # Use of regex to select sample id and cdrh3 sequence
        with open('/d/projects/u/tg001/CSV/programs/'\
          '5_Results/Fields/fields_file_#' + str(seqset)+'.txt','rt') as file:
            for line in file:
                cdrh3_flag = re.compile(r"^((\d+)(\s+)"\
                  "([^\s]+)(\s+)([A-Z]+)(\s+))", re.M)
                cdrh3_found = re.search(cdrh3_flag, line)
                if cdrh3_found != None:        
                    samp_name = cdrh3_found.group(4)
                    cdrh3_seq = cdrh3_found.group(6)

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
                with open('/d/projects/u/tg001/CSV/programs/'\
                  '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nSUBMIT TO MEME:\n\nDictionary (n =', length,')',\
                    '\n', meme_sub1, '\n\nAs FASTA file.\n', file = l)

                # FASTA file prepared for submission to meme
                print('\nFASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub1.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/CSV/programs/5_Results'\
                     '/FASTA/fasta_json#'+str(seqset)+'.fa','a') as s:
                       print('\n','>',meme_samp,'\n',meme_cdrh3,sep='',file=s)
         
            elif len(meme_cdrh3s) > 40000:
                length = len(meme_cdrh3s)
                print('\nSince number of sequences was', length, '\nThe '\
                  'dictionary was reduced: 40000 randomly selected sequences')

                # Cumulative list of samples for meme
                meme_dict = dict(zip(meme_samps, meme_cdrh3s))
                meme_sub2 = dict(random.sample(meme_dict.items(), 40000))
                with open('/d/projects/u/tg001/CSV/programs/'\
                  '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nSUBMIT TO MEME:\n\nDictionary (n = 40000)\n',\
                    meme_sub2, '\n\nAs FASTA file.\n', file = l)

                # FASTA file prepared for submission to meme
                print('\nFASTA file:')             
                for meme_samp, meme_cdrh3 in meme_sub2.items():
                    print('\n','>', meme_samp, '\n', meme_cdrh3, sep='')                
                    with open('/d/projects/u/tg001/CSV/programs/5_Results'\
                     '/FASTA/fasta_csv#'+str(seqset)+'.fa','a') as s:
                        print('\n','>',meme_samp,'\n',meme_cdrh3,sep='',file=s)
            else:
                print('\nToo few sequences for submission to MEME')
                with open('/d/projects/u/tg001/csv/programs/'\
                   '5_Results/CDRH3s/cdrh3_#' + str(seqset) + '.txt','a') as l:
                    print('\nToo few sequences for submission to MEME\n', \
                       file = l)

# Directory is defined: data file is placed here before running the script
if __name__ == '__main__':
    directory = '/d/projects/u/tg001/CSV/data'

    # Files are labelled and functions are called   
    for f in list_filepaths(directory):
        samp_names = []
        cdrh3_seqs = []
        meme_samps = []
        meme_cdrh3s = []
        meme_sub2 = {}
        meme_sub1 = {}
        meme_samp = ''
        meme_cdrh3 = ''
        seqset += 1
        print('\nFile_ref_#', seqset, '\n', f, sep='')
        with open('/d/projects/u/tg001/CSV/programs/'\
         '5_Results/Paths/path_file.txt', 'a') as p:
            print('\nFile_ref_#', seqset, '\n', f, sep='', file=p)   
        ungzip_file(f)
        metadata(f)
        select_fields(f)
        cdrh3_regex(f)
        print('\n-------------\n\nend_file_#', seqset, '\n', \
              '\n-------------\n', sep='')

