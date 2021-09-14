#!/bin/bash

# ---------------------------------------------
# IMPORT TSV FILES & EXTRACT CDR-H3 FIELDS
# ---------------------------------------------

# This script takes as input a folder containing zip archives (.zip) downloaded
# from the iReceptor Gateway repository. First, the data for analysis is copied
# into folder /d/projects/u/tg001/TSV/data ready to run this script from the
# folder /d/projects/u/tg001/TSV/programs
#
# After determining filepaths of files held by the defined directory, this
# script unzips any ZIP archive (.zip) files, Each ZIP archive represents one
# sample data unit, and contains three associated files:
#
# (1) The main TSV (.tsv) file, which contains the sample data unit itself i.e.
# one set of pre-processed Ig-Seq sequences;
# (2) A JSON (.json) file containing associated sample metadata and study
# metadata in an extended form;
# (3) A TEXT (.txt) file containing associated sample metadata and study
# metadata in a shortened form.
#
# (4) In addition to the ZIP archives representing each sample data unit, there
# is also a zip archive containing general study metadata in the form of a TSV
# file, which was also processed and saved.
#
# Many fields of processed data are provided by iReceptor Gateway, but the
# specific aim here was to extract all sequence_ids & associated CDR-H3
# sequences from each sample unit. This was done by loading the data from each
# .tsv file into a PANDAS dataframe. Sample metadata was stored by saving the
# .json and .txt files directly to the server.
# Once the script has run, data from each ZIP file has been streamed into the
# directory /d/projects/u/tg001/TSV/programs as separate folders:
# 
# 1_archpaths:    Having a single file listing all filepaths to each archive
# 2_archfiles:    Numbered files listing files contained in each zip archive
# 3_filesort:
#     meta_files: Numbered files containing metadata for each sample unit
#     summ_file:  Numbered files containing brief metadata on each sample unit
#     seq_files:  Numbered folders containing a sequence file for each sample
#                 unit in 2 different file formats - a TSV file and a CSV file
# 4_regexes:
#     fields_file:   Numbered files containing data fields extracted from
#                    Pandas dataframes (sequence_id and junction_aa)
#     meta_seqs:     A FASTA file of CDR-H3 sequences from each file, stored on
#                    the server ready for entry into a MEME Suite application
# -----------------------------------------------------------------------------

import os
import re
import zipfile
from zipfile import ZipFile
import pandas as pd
pd.set_option('display.max_rows', None)

seqset = 0

# Define function to find filepaths of zip archives
def list_archpaths(directory):
   for dirpath,_,filenames in os.walk(directory):
      for f in filenames:
         yield os.path.abspath(os.path.join(dirpath, f))

# Define function to display names of files in zip archive and sort
def unzip_arch(src):

   # A ZipFile Object is created
   if f.endswith('.zip'):
      with ZipFile(f,'r') as zipObj:
         zip_contents = zipObj.namelist()

         # Files from the unzipped archive are listed
         print('\nzip Archive_#' + str(seqset) + \
            ' contained the following file(s):\n')
         with open('/d/projects/u/tg001/TSV/programs/'\
             '2_archfiles/arch_files_#' + str(seqset) + '.txt', 'w') as z:
            for filename in zip_contents:
               print('    ', filename)
               with open('/d/projects/u/tg001/TSV/programs/'\
                  '2_archfiles/arch_files_#' + str(seqset) + '.txt', 'a') as z:
                    print(filename, file=z)
                      
               # The files are sorted by filetype
               if filename.endswith('.txt'):
                  zipObj.extract(filename, '3_filesort/summ_files/'\
                     'summ_file_#' + str(seqset))
               elif filename.endswith('.json'):
                  zipObj.extract(filename, '3_filesort/meta_files/'\
                     'meta_file_#' + str(seqset))
               elif filename.endswith('.tsv'):          
                  zipObj.extract(filename, '3_filesort/seq_files/'\
                     'seq_file_#' + str(seqset))
               else:           
                  zipObj.extract(filename, '3_filesort/misc_files/'\
                     'unkwn_file_#' + str(seqset))
                  
         print('\nfor transfer to the appropriate folder(s).\n')

# Define function to create pandas dataframes from TSV file contents
def select_fields(src):

   # Saved TSV sequence files are selected
   print('\nPandas dataframe - column headings and selected '\
      'column contents:\n')       
   with open('/d/projects/u/tg001/TSV/programs/'\
     '4_regexes/fields_file_#' + str(seqset) + '.txt', 'a') as d:    
      print('Seq_file_#', seqset, ' Pandas df - selected columns: \n'\
         'Sequence_id & junction_aa\n', f, sep='', file=d)

   # TSV file is converted and saved as CSV file for printing
   df0 = pd.read_csv('/d/projects/u/tg001/TSV/programs/3_filesort/'\
     'seq_files/seq_file_#'+str(seqset)+'/airr-covid-19-1.tsv',sep='\t')
   df1 = df0.to_csv('/d/projects/u/tg001/TSV/programs/3_filesort/'\
     'seq_files/seq_file_#'+str(seqset)+'/airr-covid-19-1.csv')
   df = pd.read_csv('/d/projects/u/tg001/TSV/programs/3_filesort/'\
     'seq_files/seq_file_#'+str(seqset)+'/airr-covid-19-1.csv')
   print(df.columns.values, '\n')
   
   # Sequence_id and junction_aa columns are selected
   print('\n', df[['sequence_id', 'junction_aa']],'\n', \
      '\n\n-------------', sep='')
   with open('/d/projects/u/tg001/TSV/programs/'\
     '4_regexes/fields_file_#' + str(seqset) + '.txt', 'a') as d:
      print('\n', df[['sequence_id','junction_aa']], '\n-------------',\
         '\nend_file_#', seqset, '\n-------------', sep='', file=d)

# Define function to determine cdrh3 sequences and save them in fasta format
def junc_regex(src):

   # Use of regex to select sample id and determine junction sequence
   print('\nMEME CDR-H3 FASTA samples:\n--------------------------')
   with open('/d/projects/u/tg001/TSV/programs/'\
      '4_regexes/fields_file_#' + str(seqset) +'.txt', 'rt') as file:
    
      for line in file:
         junc_flag = re.compile(r"^((\d+)(\s+)"\
                  "([^\s]+)(\s+)([A-Z]+)(\s+))", re.M)
         junc_found = re.search(junc_flag, line)
         if junc_found != None:
            samp_name = junc_found.group(4)
            junc_seq = junc_found.group(6)

            # Slice junction sequence for cdrh3_seq
            junc_len = len(junc_seq)
            cdrh3_len  = junc_len-2
            cdrh3_seq = junc_seq[0:cdrh3_len]
            
            # cdrh3 sequences for each sample are recorded
            samp_names.append(samp_name)  
            cdrh3_seqs.append(cdrh3_seq)

            # short cdrh3 sequences are removed (meme minimum 8 aa)  
            if len(cdrh3_seq) > 7:
               meme_samp = samp_name
               meme_cdrh3 = cdrh3_seq                     
               print('\n','>', meme_samp, sep='')
               with open('/d/projects/u/tg001/TSV/programs/4_regex'\
                  'es/meme_seqs_#'+str(seqset)+'.txt', 'a') as s:
                  print('\n','>', meme_samp, sep='', file=s)
               meme_samps.append(meme_samp)  
               print(meme_cdrh3, '\n')
               with open('/d/projects/u/tg001/TSV/programs/4_regex'\
                  'es/meme_seqs_#' + str(seqset)+'.txt','a') as c:
                  print(meme_cdrh3, file = c)
               meme_cdrh3s.append(meme_cdrh3)

   # Cumulative lists of outputs are also printed
   print('\n--------------------------')
   print('\nALL SAMPLE OUTPUTS:\n\nSample_id list:\n', samp_names)
   print('\nCDR-H3 sequence list:\n', cdrh3_seqs, '\n')
   with open('/d/projects/u/tg001/TSV/programs/'\
      '4_regexes/samp_list_#' + str(seqset) + '.txt','a') as l:
      print('seq_file_#', seqset, '\n', f, sep='', file = l)
      print('\nALL SAMPLE OUTPUTS:\n\nSample_id list:\n', samp_names, \
         '\n\nCDR-H3 sequence list:\n', cdrh3_seqs, '\n', file = l)
                    
   # Cumulative list of meme samples
   print('\nSUBMITTED TO MEME:\n\nMEME_sample list:\n', meme_samps, '\n')
   print('MEME_CDR-H3 list:\n', meme_cdrh3s, '\n') 
   with open('/d/projects/u/tg001/TSV/programs/'\
     '4_regexes/samp_list_#' + str(seqset) + '.txt','a') as l:
      print('SUBMITTED TO MEME:\n\nMEME sample list:\n', meme_samps, \
        '\n\nMEME CDR-H3 list:\n', meme_cdrh3s, '\n', file = l)

                           
# Directory is defined: data file is placed here before running the script
if __name__ == '__main__':
   directory = '/d/projects/u/tg001/TSV/data'
  
   # Archives are labelled and functions are called
   for f in list_archpaths(directory):
      samp_names = []
      cdrh3_seqs = []
      meme_samps = []
      meme_cdrh3s = []
      cdrh3_seq = ''
      meme_samp = ''
      meme_cdrh3 = ''
      seqset += 1
      print('\nArchive_#', seqset, '\n', f, sep='')
      with open('/d/projects/u/tg001/TSV/programs/'\
         '1_archpaths/path_arch.txt', 'a') as p:
         print('\nArchive_#', seqset, '\n', f, sep='', file=p)
      unzip_arch(f)
      select_fields(f)
      junc_regex(f) 
      print('---------------\n\nend_Seq_file_#', seqset, '\n', \
          '\n---------------\n', sep='')


