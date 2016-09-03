# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 16:45:18 2016

@author: Gungnir
"""
#==============================================================================
# search and batch download .fasta(Protein) from NCBI
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
keyword = 'plk1'
storage_dir = '_'.join([keyword,'fastaP'])
#not "predicted"

filter_search = '[Title] and "mrna" not "BAC" not "whole genome"  not "partial" not "synthetic" not "TSA" not "STS" not "patent" not "complete sequence"'
handle = entrez.esearch(db="nucleotide",retmax=2500,term=''.join([keyword,filter_search]) )
match_list = entrez.read(handle) 
print(match_list)


if not os.path.isdir(storage_dir):
    os.mkdir(storage_dir)

os.chdir(storage_dir)

for gi in match_list['IdList']:
    temp_feat=[]
    temp_fasta= []
    content= [] 
    handle_fasta = entrez.efetch(db="nucleotide", id = gi, rettype="fasta", retmode="text")
    temp_fasta = handle_fasta.read()
    handle_fasta.close()
#    print(temp_fasta)
    handle_gb = entrez.efetch(db="nucleotide", id = gi, rettype="gb", retmode="text")
    rec = SeqIO.read(handle_gb,'gb')
    for feature in rec.features:
#        print(rec.features)
        if feature.type == "CDS":
            temp_feat = feature.qualifiers['translation']
        #    print(temp_fasta)
    reg_pattern = ">gi\|[0-9]*\|[a-zA-Z]*\|[A-Z0-9._]*\|.*\n?(?![ATCG])"
    title = re.search(reg_pattern ,temp_fasta,re.M).group(0)
    title = re.sub('\n','',title)
    title = re.sub('>gi\|[0-9]*\|[a-zA-Z]*\|[A-Z0-9._]*\| ','',title)
    title = re.sub('/','-',title)
#    print(temp_fasta)
    content = re.sub('(?<=\n)([ATCG]*\n)*[ATCG]*\n',''.join(temp_feat),temp_fasta)
#   print(content)
    output_name = '.'.join([title,'fasta'])
    if not os.path.isfile(output_name):
        out_handle = open(output_name, "w")
        out_handle.write(content)
        out_handle.close()
        
        
#genbank
os.chdir(cur_dir)