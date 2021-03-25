#!/usr/bin/python

# (1) How many records are in the file?

f = open("/Users/jkcollection/Downloads/dna.fasta", "r")
lines = f.readlines()

c = 0

for line in lines:
    if line.startswith(">"):
        c += 1

print(c)

# (2) What are the lengths of the sequences in the file? 
# What is the longest sequence and what is the shortest sequence? 
# Is there more than one longest or shortest sequence?
# What are their identifiers?

import Bio
from Bio import SeqIO

d = {}

for record in SeqIO.parse("/Users/jkcollection/Downloads/dna.fasta", "fasta"):
    name = record.id
    l = len(record)
    d[name] = l

import operator 

sorted_d = sorted(d.items(), key = operator.itemgetter(1))
print (sorted_d)

# (3) ORFs

from Bio import SeqIO

d = {}

for record in SeqIO.parse("/Users/jkcollection/Downloads/dna.fasta", "fasta"):
    name = record.id
    s = str(record.seq)
    d[name] = s

def find_ORF(sequence, frame): # frame value is 0, 1 or 2

    start_indexs = []
    stop_indexs = []

    # START codon indexs
    for i in range(frame, len(sequence), 3):
        if sequence[i:i+3] == "ATG":
            start_indexs.append(i)

    # STOP codon indexs
    for j in range(frame, len(sequence), 3):
        if sequence[j:j+3] in ["TAA", "TGA", "TAG"]:
            stop_indexs.append(j)

    ORF = []
    mark = 0
    for i in range(0, len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                k = sequence[start_indexs[i]:stop_indexs[j]+3]
                ORF.append(k)
                mark = stop_indexs[j]+3
                break
    return ORF

def find_start_position_for_ORF(sequence, frame): # frame value is 0, 1 or 2

    start_indexs = []
    stop_indexs = []

    # START codon indexs
    for i in range(frame, len(sequence), 3):
        if sequence[i:i+3] == "ATG":
            start_indexs.append(i)

    # STOP codon indexs
    for j in range(frame, len(sequence), 3):
        if sequence[j:j+3] in ["TAA", "TGA", "TAG"]:
            stop_indexs.append(j)

    mark = 0
    start_position = {}
    for i in range(0, len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                k = sequence[start_indexs[i]:stop_indexs[j]+3]
                start_position[len(k)] = start_indexs[i]
                mark = stop_indexs[j]+3
                break
    return start_position

# What is the length of the longest ORF in the file? 
# What is the identifier of the sequence containing the longest ORF?

from collections import defaultdict

orfs_d = {}
orfs_d_len = defaultdict(list)
orfs_d_max = {}

for k, v in d.items():
    orfs = find_ORF(v, 0)
    orfs_d[k] = orfs
for k, v in orfs_d.items():
    for x in v:
        orfs_d_len[k].append(len(x))
for k, v in orfs_d_len.items():
    h = max(v)
    orfs_d_max[k] = h

longest_ORF = max(orfs_d_max, key = orfs_d_max.get)
print (longest_ORF, orfs_d_max[longest_ORF])

# What is the starting position of the longest ORF in the sequence that contains it?

s_pos = {}
for k, v in d.items():
    orfs = find_start_position_for_ORF(v, 0)
    s_pos[k] = orfs
print (s_pos[longest_ORF][orfs_d_max[longest_ORF]])


# For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? 

identifier = "gi|142022655|gb|EQ086233.1|91"
for k, v in orfs_d_max.items():
    if k == identifier:
        print (v)

#(4) Repeats

# Given a length n, the program should be able to identify all repeats of length n in all sequences in the FASTA file. 
# The program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.

from Bio import SeqIO

d = {}

for record in SeqIO.parse("/Users/jkcollection/Downloads/dna.fasta", "fasta"):
    name = record.id
    s = str(record.seq)
    d[name] = s

def repeats(sequence):
    l = len(sequence)
    repeats = []
    for i in range(l):
        repeats.append(sequence[i:i+5]) # n = 5
    return repeats

r = []
for v in d.values():
    repeats_list = repeats(v)
    for x in repeats_list:
        if len(x) == 5: # n = 5
            r.append(x)

def most_common(lst):
    return max(set(lst), key=lst.count)

print (most_common(r)) # the most frequently occuring repeat in all sequences
print (r.count(most_common(r))) # how many times does it occur in all sequences

from collections import Counter

c = Counter(r)
max_value = max(c.values())
l = [k for k,v in c.items() if v == max_value]
print (str(l) + "\n" + str(len(l))) # the most frequent repeat and how many same length repeats occur max times
