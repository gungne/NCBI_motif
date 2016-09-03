# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 16:12:00 2016

@author: Gungnir
"""

import os
import re
from Bio import SeqIO 
from skbio.alignment import local_pairwise_align_ssw as SSWAlign

def listdir_nohidden(path):
    list_t = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            list_t.append(f)
    return list_t     
      
work_dir = os.getcwd()
#print (work_dir)
db_dir = os.path.join(work_dir,'pde5a_nonpredict')
os.chdir(db_dir)
db_list = listdir_nohidden(db_dir)
#db_list = db_list.remove('.DS_Store') 
dataset=[]
for file in db_list :
    handle_read = open(file) 
    rec = SeqIO.read(handle_read, 'gb')
#    print(rec.features)
    for feature in rec.features:
        if feature.type == "CDS":
            dataset.append(feature.qualifiers['translation'])
            print(feature.qualifiers['translation'])
alignment , score = SSWAlign(dataset[0],dataset[1])
print(alignment)
#AlignIO
os.chdir(work_dir)