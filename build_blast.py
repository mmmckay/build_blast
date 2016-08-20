import os
import glob
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

fastas = glob.glob('*.fasta')

with open('all.fasta', 'w') as all:
    for file in fastas:
        for line in open(file):
            all.write(line)

os.system('makeblastdb -in all.fasta -out allblastdb -parse_seqids -dbtype prot')

for file in fastas:
    if file != 'all.fasta':
        os.system('blastp -db allblastdb -query {} -outfmt "10 qseqid sseqid pident qcovs evalue" -out {}vsall.csv'.format(file,file.split('.')[0]))

gn_dict = defaultdict(list)
for file in fastas:
    if file != 'all.fasta':
        input_fasta = SeqIO.parse(open(file, 'r'), 'fasta')
        gn_dict = defaultdict(list)
        for seq in input_fasta:
            gn_dict['id'].append(seq.id)
            gn_dict['translation'].append(''.join(seq.seq))
            gn_dict['name'].append(seq.description[seq.description.index(' ')+1:])

        df = pd.DataFrame(gn_dict)
        df = df[['name', 'id', 'translation']]
        df.to_csv('{}_genes.csv'.format(file.split('.')[0]),index=False)
