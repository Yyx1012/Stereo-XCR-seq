import sys
import pandas as pd
import gzip
import re
import multiprocessing

fastq_file = sys.argv[1]
match_file = sys.argv[2]
cdr3_file = sys.argv[3]
align_file = sys.argv[4]
modify_file = sys.argv[5]
incdr3type = sys.argv[6]

def read_fastq(fastq):
    rnames = list()
    if fastq.endswith('.gz'):
        f = gzip.open(fastq, 'rt')
    else:
        f = open(fastq, 'r')
    for idx, l in enumerate(f):
        if idx % 4 == 0:
            rnames.append(l.split()[0][1:])

    return rnames

def get_Hits(x):
    if isinstance(x,float):
        return 'NA'
    return re.sub(r'\(.*?\)', '', x)

rnames = read_fastq(fastq_file)
raw_match = pd.read_csv(match_file)
cdr3 = pd.read_csv(cdr3_file)
align = pd.read_csv(align_file,sep = '\t')


raw_match['loc'] = raw_match['bam_x'].map(str) + '_' + raw_match['bam_y'].map(str)
match = raw_match[['rname','loc']].copy()
match.columns = ['rnames','loc']
loc_dict = dict(zip(match['rnames'],match['loc']))
cdr3['umi'] = cdr3['rnames'].apply(lambda x : x.split('|||')[1].split('UB:')[1])
cdr3['loc'] = cdr3['rnames'].map(loc_dict)
match_cdr3 = cdr3.dropna()
cdr3_dict = dict(zip(cdr3['rnames'],cdr3['CDR3_aa']))
for typ in ['V','D','J','C']:
    align[f'{typ}Hits'] = align[f'all{typ}HitsWithScore'].apply(get_Hits)
align['rnames'] = align['readId'].map(lambda x : rnames[x])

simple_align = align[['rnames','VHits','DHits','JHits','CHits']].copy()
simple_align = simple_align.merge(match ,on = 'rnames',how = 'left')
simple_align = simple_align.dropna()
simple_align['cdr3'] = simple_align['rnames'].map(cdr3_dict)
simple_align = simple_align.dropna()
simple_align['umi'] = simple_align['rnames'].apply(lambda x : x.split('|||UB:')[1])
dis_dict = dict(zip(raw_match['rname'],raw_match['dis']))
simple_align['dis'] = simple_align['rnames'].map(dis_dict)
simple_align['CDR3Type'] = incdr3type
simple_align.to_csv(modify_file,index=None,sep = '\t')
