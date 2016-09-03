# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 16:12:00 2016

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
    
#==============================================================================
# Patterns
Cad_pattern = (['']*9)
Cai_pattern = (['']*5)
#Ca2+dependent
Cad_pattern[0] = '[FILVW].{8}[FILVW]'
Cad_pattern[1] = '[FILVW].{3}[FILVW].{4}[FILVW]'
Cad_pattern[2] = '[RK][RK][RK][FAILVW].{3}[FILV].{4}[FILVW]'
Cad_pattern[3] = '[FILVW].{10}[FILVW]'
Cad_pattern[4] = '[FILVW].{12}[FILVW]'
Cad_pattern[5] = '[FILVW].{6}[FAILVW].{5}[FILVW]'
Cad_pattern[6] = '[FILVW].{3}[FAILVW].{2}[FAILVW].{5}[FILVW]'
Cad_pattern[7] = '[RK][RK][RK][FILVW].{6}[FAILV].{5}[FILVW]'
Cad_pattern[8] = '[FILVW].{14}[FILVW]'
#Ca2+ independent
Cai_pattern[0] = '[FILV]Q.{3}[RK]G.{3}[RK].{2}[FILVWY]'
Cai_pattern[1] = '[FILV]Q.{3}[RK]'
Cai_pattern[2] = '[FILV]Q.{3}R.{4}[VL][KR].{1}'
Cai_pattern[3] = '[IL]Q.{2}G.{4}K.R.W'
Cai_pattern[4] = '[IVL]Q.{3}R.{4}[RK].{2}[FILVWY]'


#==============================================================================

work_dir = os.getcwd()
#print (work_dir)
db_dir = os.path.join(work_dir,'pde5a_nonpredict')
pdb_dir = os.path.join(work_dir,'pde5a_predict')
os.chdir(db_dir)
db_list = listdir_nohidden(db_dir)
#db_list = db_list.remove('.DS_Store') 
#print(db_list[0])
Cad_motifset=dict()
Cai_motifset=dict()
#==============================================================================
# round 1 filter cross out pattern dont even exist in canonical sequence
#==============================================================================
for pattern in Cad_pattern:
    file = db_list[0]
    handle_read = open(file) 
    rec = SeqIO.read(handle_read, 'gb')
#        print(rec.features)
    for feature in rec.features:
        if feature.type == "CDS":
            seq_temp = ''.join(feature.qualifiers['translation'])
            result = regex.findall(pattern ,seq_temp, overlapped = True)
#            print(result)
            if result==[]:
                Cad_pattern.remove(pattern)
            else:
                Cad_motifset[pattern]=result

#print(Cad_pattern)
#print(Cad_motifset)
for pattern in Cai_pattern:
    file = db_list[0]
    handle_read = open(file) 
    rec = SeqIO.read(handle_read, 'gb')
#        print(rec.features)
    for feature in rec.features:
        if feature.type == "CDS":
            seq_temp = ''.join(feature.qualifiers['translation'])
            result = regex.findall(pattern ,seq_temp, overlapped = True)
            if result==[]:
                Cai_pattern.remove(pattern)
            else:
                Cai_motifset[pattern]=result
                    

#==============================================================================
# round 2 filter for rest of the canonical 
#==============================================================================
handle_read = open(db_list[0]) 
rec = SeqIO.read(handle_read, 'gb')
for feature in rec.features:
    if feature.type == "CDS":
        seq_ref = ''.join(feature.qualifiers['translation'])

        
for pattern in Cad_motifset:
    for motif in Cad_motifset[pattern]:
        for file in db_list :
            if file == db_list[0]:
                continue
            else:
                handle_read = open(file) 
                rec = SeqIO.read(handle_read, 'gb')
                for feature in rec.features:
                    if feature.type == "CDS":
                        seq_temp = ''.join(feature.qualifiers['translation'])
#                        print(seq_temp)
                result_temp = regex.findall(motif,seq_temp, overlapped = True) 
#                still got a chance that motif are exactly the same
                if result_temp == []:
#                    print(motif)
                    try:
                        Cad_motifset[pattern].remove(motif)
                    except:
                        continue

for pattern in Cai_motifset:
    for motif in Cai_motifset[pattern]:
        for file in db_list :
            if file == db_list[0]:
                continue
            else:
                handle_read = open(file) 
                rec = SeqIO.read(handle_read, 'gb')
                for feature in rec.features:
                    if feature.type == "CDS":
                        seq_temp = ''.join(feature.qualifiers['translation'])
#                        print(seq_temp)
                result_temp = regex.findall(motif,seq_temp, overlapped = True) 
#                still got a chance that motif are exactly the same
                if result_temp == []:
#                    print(motif)
                    try:
                        Cai_motifset[pattern].remove(motif)
                    except:
                        continue

#in the case of the e
print(Cai_motifset)
#==============================================================================
# round 3 running through some predicted data
#==============================================================================

os.chdir(pdb_dir)
pdb_list = listdir_nohidden(pdb_dir)

scoreboard=dict()

for pattern in Cad_motifset:
    for motif in Cad_motifset[pattern]:
        scoreboard[motif]=0;
        for file in pdb_list :
            if file == db_list[0]:
                continue
            else:
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