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
target1_pattern = (['']*3)
target2_pattern = (['']*3)
#Ca2+dependent
target1_pattern[0] = '.{1}[S][ST][P].{1}'
target1_pattern[1] = '.{2}[S][ST][P]'
target1_pattern[2] = '[S][ST][P].{2}'
#Ca2+ independent
target2_pattern[0] = '.{1}[S][ST].{1}'
target2_pattern[1] = '.{2}[S][ST]'
target2_pattern[2] = '[S][ST].{2}'




#==============================================================================
os.chdir('Data')
work_dir = os.getcwd()
#print (work_dir)
keyword='HNRNPH1'
db_dir = os.path.join(work_dir,'_'.join([keyword,'nonpredict']))
pdb_dir = os.path.join(work_dir,'_'.join([keyword,'predict']))
os.chdir(db_dir)
db_list = listdir_nohidden(db_dir)
#db_list = db_list.remove('.DS_Store') 
#print(db_list[0])
target1_motifset=dict()
target2_motifset=dict()
#==============================================================================
# round 1 filter cross out pattern dont even exist in canonical sequence
#==============================================================================
for pattern in target1_pattern:
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
                target1_pattern.remove(pattern)
            else:
                target1_motifset[pattern]=result

#print(target1_pattern)
#print(target1_motifset)
for pattern in target2_pattern:
    file = db_list[0]
    handle_read = open(file) 
    rec = SeqIO.read(handle_read, 'gb')
#        print(rec.features)
    for feature in rec.features:
        if feature.type == "CDS":
            seq_temp = ''.join(feature.qualifiers['translation'])
            result = regex.findall(pattern ,seq_temp, overlapped = True)
            if result==[]:
                target2_pattern.remove(pattern)
            else:
                target2_motifset[pattern]=result
                    

#==============================================================================
# round 2 filter for rest of the canonical 
#==============================================================================
handle_read = open(db_list[0]) 
rec = SeqIO.read(handle_read, 'gb')
for feature in rec.features:
    if feature.type == "CDS":
        seq_ref = ''.join(feature.qualifiers['translation'])

        
for pattern in target1_motifset:
    for motif in target1_motifset[pattern]:
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
                        target1_motifset[pattern].remove(motif)
                    except:
                        continue
print(target1_motifset)
for pattern in target2_motifset:
    for motif in target2_motifset[pattern]:
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
                        target2_motifset[pattern].remove(motif)
                    except:
                        continue

#in the case of the e
print(target2_motifset)
#==============================================================================
# round 3 running through some predicted data
#==============================================================================

os.chdir(pdb_dir)
pdb_list = listdir_nohidden(pdb_dir)

scoreboard=dict()

for pattern in target1_motifset:
    for motif in target1_motifset[pattern]:
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

for pattern in target2_motifset:
    for motif in target2_motifset[pattern]:
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
                        
print(len([name for name in os.listdir('.') if os.path.isfile(name)]))
#target1_pattern.remove(pattern)
#print(target1_motifset)
#score_sorted= sorted(scoreboard.iteritems(), key=lambda d:d[1], reverse = True)
print(scoreboard)
#            print(feature.qualifiers['translation'])
#print(alignment)
#AlignIO
os.chdir(work_dir)