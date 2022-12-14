# Marlow
# This code is designed to take an input file and do the following:
# 1. Take as input a DNA sequence
# 2. Identify simulated interactions - copying & sequencing
# 3. replicate for storage and introduce errors for synthesis, decay, and PCR processes
# 4. Export fasta files for each sample, incorporating errors into sequence and errors to generate quality scores

import numpy as np
import random

BASE = np.array(['A','C','G','T'])
QUANT = {'A': 0, 'C':1, 'G':2, 'T':3}

# Synthesis routine
def synthesize(dna):
    #prob_s = 0.004
    prob_i = 0.00081
    prob_d = 0.0024
    prob_s_matrix = [
        [0.999865, 0.00002, 0.000011, 0.000104],
        [0.00002, 0.99976, 0.00021, 0.00001],
        [0.000021, 0.0001535, 0.9998165, 0.000009],
        [0.002569, 0.000021, 0.000243, 0.997167],
    ]
    print('Synthesizing')
    return errorTool(dna,prob_i,prob_d,prob_s_matrix)


def copy_routine(dnas):      # doubles the amount of oligos, selects half at random to return (should be 1024)
    cpy = pcr(dnas, 1)
    print('Copying')
    return np.random.choice(np.array(cpy), int(len(cpy)/2), replace=False)


def sequence_routine(dnas, depth, naming):     # sequences the oligos at a defined depth, returns 1024 selected at random
    copy = pcr(dnas, depth)
    random.shuffle(copy)      # shuffle the oligos
    rem = copy[:1024]         # remainder keeps 1024 oligos
    readseq = copy[1024:]     # the rest to be sequenced

    # write fasta of sequenced oligos
    print('writing sequencing alignment file')
    align = open('output2/sequences-' + naming + '.fasta', 'w')
    for idx, i in enumerate(readseq):
        if idx <= 5000:          # Limit alignment to first 5000
            align.write(">ih|seq" + naming + "|oligo" + str(idx) + "\n")
            align.write(i + '\n')
    align.close()
    print('output2/sequences-' + naming + '.fasta successfully written')

    #return remainder
    return rem


# PCR routine
def pcr(dnas, cycles):
    prob_i = 0.0024
    prob_d = 0.0024

    prob_s_matrix = [
        [0.9997, 0.00001, 0.00027, 0.00002],
        [0.00001, 0.99974, 0, 0.00025],
        [0.00065, 0, 0.99935, 0],
        [0.00001, 0.00064, 0.00001, 0.99934],
    ]
    new = []
    for x in range(cycles):
        new.clear()
        for dna in dnas:
            new.append(errorTool(dna, prob_i, prob_d, prob_s_matrix))
        dnas = dnas + new
    return dnas

def errorTool(d,prob_i,prob_d,prob_s_matrix):
    res = []
    nts = list(d)
    # Determine substitution
    for nt in nts:
        for i, base in enumerate(['A', 'C', 'G', 'T']):
            if nt == base:
                sub = np.random.choice(['A', 'C', 'G', 'T'], p=prob_s_matrix[i])
                if sub != nt:
                    nt = sub
        res.append(nt)  # result with substitution errors
    # Deletion and Insertion
    delP = np.where(np.random.choice([False, True], size=len(nts), p=[1 - prob_d, prob_d]))[0]
    insP = np.where(np.random.choice([False, True], size=len(nts), p=[1 - prob_i, prob_i]))[0]
    for j, pos in enumerate(delP): (res.pop(pos-j))
    for pos in insP: res.insert(pos, np.random.choice(['A', 'T', 'C', 'G']))

    # return string version of sequence
    return ''.join(str(x) for x in res)

# Intake the DNA encoded data sequence
orig = input("Input Original Sequence: ")

copies = int(input("How many times will the sequence be copied/sequenced? "))
toseq = [int(x) for x in input("Which steps will be sequencing? ").split()]

# Build copy/sequence guide list
guide = ['c'] * (copies + 1)    # add one to account for the original
for i in toseq: guide[i] = 's'
print(guide)

# synthesys
syn = [synthesize(orig)]      # create primary synthesized oligo
print(syn)
synbk = [synthesize(orig)]    # create backup synthesized oligo

print('running initial PCR')
synpcr = pcr(syn, 10)
synbkpcr = pcr(synbk, 10)

currentset = []
currentset = currentset + synpcr
#print(currentset)

# Sequence backup
sequence_routine(synbkpcr, 5, "syn_backup")

# Run through copy/sequencing plan

for idx, i in enumerate(guide):
    if i == 'c':    # copy
        currentset = copy_routine(currentset).tolist()

    if i == 's':    # sequence
        name = 'step' + str(idx)
        currentset = sequence_routine(currentset, 5, name)
    #print(currentset)

print('processing complete')