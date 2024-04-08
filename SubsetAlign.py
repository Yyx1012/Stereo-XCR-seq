import pandas as pd
import gzip
from glob import glob
import os
import sys 

def read_fastq(cln,fastq):
    rnames = list()
    if fastq.endswith('.gz'):
        f = gzip.open(fastq, 'rt')
    else:
        f = open(fastq, 'r')
    for idx, l in enumerate(f):
        if idx % 4 == 0:
            rnames.append([cln,l.split()[0][1:]])
    return rnames


sampleid = sys.argv[1]
mixcr_align = sys.argv[2]  # "./mc38_0304/mc38_0304.align.tsv"
fastq_need_subset = sys.argv[3]    #"/storage/liuyi/09.ma_tcr/p1/split_fix/MC38-TRBC_mixcr_CDR3.r2.fastq.gz"
clone_fastq_dir = sys.argv[4]   #"./mc38_0304/clone_fq"
clone_info = sys.argv[5]      #"./mc38_0304/mc38_0304.contigs.tsv"
outdir = sys.argv[6]
barcode_file = sys.argv[7]   # "/storage/liuyi/09.ma_tcr/p1/split_fix/barcode.forward.csv.gz"


out_fastq = f"{outdir}/{sampleid}.align.r2.fastq.gz"
out_barcode = f"{outdir}/{sampleid}.align.barcode.csv.gz"
out_read_clone = f"{outdir}/{sampleid}.ReadsCDR3.csv"

# load file
mixcr_align = pd.read_csv(mixcr_align,sep = '\t')
clone_info = pd.read_csv(clone_info,sep = '\t')
barcode = pd.read_csv(barcode_file,compression='gzip',header=None,names = ['rname','barcode'])
cln_fastqs = glob(f'{clone_fastq_dir}/*.gz')

# get clone read in fastq_dir
total_clone_df_list = []
for i in clone_info['cloneId']:
    basename = os.path.basename(cln_fastqs[i])
    cln_num = int(basename.split('.')[0].split(sampleid+'_')[1].split('cln')[1])
    rnames = read_fastq(cln_num,cln_fastqs[i])
    clone_df =  pd.DataFrame(rnames,columns=['clone','rnames'])
    total_clone_df_list.append(clone_df)
total_clone_df = pd.concat(total_clone_df_list,ignore_index=True)
total_clone_df = total_clone_df.sort_values('clone',ignore_index=True)
subset_barcode = barcode[barcode['rname'].isin(total_clone_df['rnames'])].copy()

# write
align_dict = dict(zip(subset_barcode.index,[1] * subset_barcode.shape[0]))
if fastq_need_subset.endswith('.gz'):
    f = gzip.open(fastq_need_subset, 'rt')
else:
    f = open(fastq_need_subset, 'r')

with gzip.open(out_fastq,'wb') as fq_f:
    for idx, l in enumerate(f):
        
        if align_dict.get(idx//4,False):
            if idx % 4 == 0:
                rnames = l.split()[0][1:]
                fq_f.write(f'@{rnames}\n'.encode('utf-8'))
            if idx % 4 == 2:
                fq_f.write('+\n'.encode('utf-8'))
            if idx % 4 == 1:
                seqs = l.rstrip()
                fq_f.write(f'{seqs}\n'.encode('utf-8'))
            if idx % 4 == 3:
                qs = l.rstrip()
                fq_f.write(f'{qs}\n'.encode('utf-8'))
                
subset_barcode.to_csv(out_barcode,compression = 'gzip',header=None,index=None)

#stat
clone_info['func'] = clone_info['aaSeqCDR3'].apply(lambda x : True if '_' in x or '*' in x else False)
total_read = barcode.shape[0]
align_read = mixcr_align.shape[0]
clone_read = subset_barcode.shape[0]
func_dict = clone_info['func'].value_counts().to_dict()
total_clone = func_dict[True]+func_dict[False]
clone_func_dict = dict(zip(clone_info['cloneId'],clone_info['func']))
total_clone_df['func'] = total_clone_df['clone'].map(clone_func_dict)

with open(f'{outdir}/step2.Mixcr.report','w') as f:
    f.write(f'## MIXCR Reads ##\nTotalReads\t{total_read}\t100.00%\nAlignReads\t{align_read}\t{align_read/total_read*100:.2f}%\nCloneReads\t{clone_read}\t{clone_read/total_read*100:.2f}% <<-\n\n')
    f.write(f'## Clone Type ##\nTotalClone\t{total_clone}\t100.00%\nFunc_Clone\t{func_dict[False]}\t{func_dict[False]/total_clone*100:.2f}%\nNonFunc_Clone\t{func_dict[True]}\t{func_dict[True]/total_clone*100:.2f}%\n\n')

total_clone_df['CDR3_aa'] = total_clone_df['clone'].map(dict(zip(clone_info['cloneId'],clone_info['aaSeqCDR3'])))
total_clone_df['CDR3_n'] = total_clone_df['clone'].map(dict(zip(clone_info['cloneId'],clone_info['targetSequences'])))
total_clone_df.to_csv(out_read_clone,index=None)
