# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:31:12 2016

@author: Gungnir
"""

import os
import re
from Bio import SeqIO 
from Bio import AlignIO
import regex 

def listdir_nohidden(path):
    list_t = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            list_t.append(f)
    return list_t   
    

work_dir = os.getcwd()
#print (work_dir)
ref_seq = 'MYRRTSNMVGLSYPPAVIEALKDVDKWSFDVFSLNEASGDHALKFIFYELLTRYDLISRFKIPISALVSFVEALEVGYSKHKNPYHNLMHAADVTQTVHYLLYKTGVANWLTELE'

pdb_dir = os.path.join(work_dir,'pde1c_predicted')
#os.chdir(db_dir)
#db_list = listdir_nohidden(db_dir)
#db_list = db_list.remove('.DS_Store') 
#print(db_list[0])
#Cad_motifset=dict()
#Cai_motifset=dict()

oligo_motif= regex.findall('.{10}',ref_seq, overlapped = False) 

os.chdir(pdb_dir)
pdb_list = listdir_nohidden(pdb_dir)

scoreboard=dict()

    
for motif in oligo_motif:
    scoreboard[motif]=0;
    for file in pdb_list :
#            if file == db_list[0]:
#                continue
#            else:
        handle_read = open(file) 
        rec = SeqIO.read(handle_read, 'gb')
        for feature in rec.features:
            if feature.type == "CDS":
                seq_temp = ''.join(feature.qualifiers['translation'])
#                        print(seq_temp)
            result_temp = regex.findall(motif,seq_temp, overlapped = True) 
#                still got a chance that motif are exactly the same
        if result_temp != []:
            scoreboard[motif]=scoreboard[motif]+1


#Cad_pattern.remove(pattern)
#print(Cad_motifset)
#score_sorted= sorted(scoreboard.iteritems(), key=lambda d:d[1], reverse = True)
print(scoreboard)
#            print(feature.qualifiers['translation'])
#print(alignment)
#AlignIO
os.chdir(work_dir)