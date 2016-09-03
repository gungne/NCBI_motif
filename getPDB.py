# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 12:16:35 2016

@author: Gungnir
"""

import os
import re
from Bio import ExPASy
from Bio import SwissProt
import gzip
from Bio.SwissProt import KeyWList
import urllib


keywlist= open("keywlist.txt")
records = KeyWList.parse(keywlist)
for record in records:
    print(record['ID'])
    print(record['DE'])

work_dir = os.getcwd()
#data_dir = os.path.join()

db_pdb = open('uniprot_sprot.dat')
#descriptions = [record.description for record in SwissProt.parse(db_pdb)]
print(dir(db_pdb))
