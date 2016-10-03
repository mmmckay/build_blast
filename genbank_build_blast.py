from Bio import SeqIO
import os
import glob
from multiprocessing import Pool as ThreadPool
from datetime import datetime

startTime = datetime.now()
base = os.getcwd()
print(base)
os.chdir(base+'/genbank_input')
os.mkdir(base+'/output')
gbk_files = glob.glob('*.gbk')
seqids = []

for file in gbk_files:
    input_gbk = SeqIO.parse(open('{}'.format(file), 'r'), 'genbank')
    db_output = open('all.fasta', 'a')
    indv_output = open('{}.fasta'.format(file.split('.')[0]),'w')

    count = 0
    for sequence in input_gbk:
        for feature in sequence.features:
            if feature.type == 'CDS':
                if feature.qualifiers['locus_tag'][0] not in seqids:
                    try:
                        db_output.write('>{}\n{}\n'.format(feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0]))
                        indv_output.write('>{}\n{}\n'.format(feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0]))
                        count += 1
                    except KeyError:
                        continue
                seqids.append(feature.qualifiers['locus_tag'][0])
    print(count)

os.system('makeblastdb -in {}/all.fasta -out {}/protblastdb -parse_seqids -dbtype prot'.format(base+'/genbank_input',base+'/output'))

fasta_files = glob.glob('*.fasta')

blastp_list = []
for file in fasta_files:
    if file != 'all.fasta':
        string = 'blastp -db {} -query {} -outfmt "10 qseqid sseqid pident qcovs evalue" -out {}vsall.csv'.format(base+'/output/protblastdb',base+'/genbank_input/'+file,file.split('.')[0])
        blastp_list.append(string)
pool = ThreadPool(len(fasta_files))
pool.map(os.system, blastp_list)
pool.close()
pool.join()

print(datetime.now()-startTime)