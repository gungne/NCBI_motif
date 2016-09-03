# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 16:45:18 2016

@author: Gungnir
"""
#==============================================================================
# search and batch download .gb from NCBI
#==============================================================================
import os
import re
from Bio import Entrez as entrez
from Bio import SeqIO

cur_dir = os.getcwd()
entrez.email = 'zhxzcj@hotmail.com'
entrez.tool = 'getGenBank'
#handle = entrz.esearch(db='pub')
print(cur_dir)
keyword = 'pde1c'
storage_dir = '_'.join([keyword,'predicted'])
filter_search = ' and "mRNA"  not "BAC" not "whole genome"  not "partial" not "synthetic" not "TSA" not "STS" not "patent" not "complete sequence"'
handle = entrez.esearch(db="nucleotide",term=' '.join([keyword,filter_search]) ,RetMax = '1000')
match_list = entrez.read(handle) 
print(match_list)


if not os.path.isdir(storage_dir):
    os.mkdir(storage_dir)

os.chdir(storage_dir)

for gene_key in match_list['IdList']:
    handle_gb = entrez.efetch(db="nucleotide", id = gene_key, rettype="gb", retmode="text")
    temp_gb =handle_gb.read()   
    reg_pattern = "(?<=DEFINITION {2}).*(\n?(?!ACCESSION)(?<=\t)*.*)*"
    title = re.search(reg_pattern ,temp_gb,re.M).group(0)
    title = re.sub('\n','',title)
    title = re.sub(' {6}','',title)
    title = re.sub('/','-',title)
    print(title)
    output_name = '.'.join([title,'gb'])
    if not os.path.isfile(output_name):
        out_handle = open(output_name, "w")
        out_handle.write(temp_gb)
        out_handle.close()
        
        
#genbank
os.chdir(cur_dir)