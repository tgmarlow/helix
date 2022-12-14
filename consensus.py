# Thomas Marlow
# Tool to create a consensus alignment from fasta files provided in a directory
# Padding is provided as needed to ensure sequence lengths are equal

import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from Bio import Seq
import os

file_list = os.listdir('output/')
print(file_list)
input('ok?')

final_align = []
seq_names = []
for x in file_list:
    if x.endswith(".fasta"):
        input_file = 'output/'+ x
        records = SeqIO.parse(input_file, 'fasta')
        records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen
        maxlen = max(len(record.seq) for record in records)
        print('max length:' + str(maxlen))
        # pad sequences so that they all have the same length
        for record in records:
            if len(record.seq) != maxlen:
                sequence = str(record.seq).ljust(maxlen, '.')
                record.seq = Seq.Seq(sequence)
        assert all(len(record.seq) == maxlen for record in records)

        # write to temporary file and do alignment
        output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
        with open(output_file, 'w') as f:
            SeqIO.write(records, f, 'fasta')
        alignment = AlignIO.read(output_file, "fasta")
        print('padded alignment written for ' + input_file)

        alignment = AlignIO.read(os.path.splitext(input_file)[0]+'_padded.fasta', 'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus(float(0.5))
        print(consensus)
        final_align.append(consensus)
        seq_names.append(x)


# make the multiple alignment of the consensus sequences

print('writing consensus alignment file')
align = open('output/consensus_alignment.fasta', 'w')
for idx, i in enumerate(final_align):
    align.write(">ih|file-" + seq_names[idx] + "|consensus" + "\n")
    #i = i[:199]     # truncate to 199 bases to keep length consistent
    align.write(str(i) + '\n')
align.close()
print('consensus sequences fasta successfully written')


# increase initial sequence to 5x
#x reduce sequencing to 5000
#x try emboss cons to see if consensus is better
