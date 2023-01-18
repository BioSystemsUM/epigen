import pandas as pd
import os

home_path = os.getcwd()

METH_FLD = 'data/methylation'
COEF_PTH =

fld = os.listdir(METH_FLD)
fld.remove('DNA__Illumina_450K_methylation_Gene_average.xls')
fld.remove('CCLE_RRBS_TSS1kb_20181022.txt')
f_d = dict()
for f in fld:
    fln = f.split('_')[2]
    df = pd.read_csv(os.path.join(METH_FLD, f), sep='\t', index_col=0, skipinitialspace=True)
    df = df.dropna(how='all').drop('CpG_sites_hg19', axis=1)
    df = df.iloc[:, 1:].apply(lambda x: ((x * df.avg_coverage).sum())/(df.avg_coverage.sum()) ,axis=0) # weights don't add up to one/100, so need to divide by sum of weights
    f_d[fln] = df
a = pd.DataFrame(f_d)
#
fld = os.listdir(METH_FLD)
fld.remove('DNA__Illumina_450K_methylation_Gene_average.xls')
fld.remove('CCLE_RRBS_TSS1kb_20181022.txt')
f_l = list()
for f in fld:
    fln = f.split('_')[2]
    df = pd.read_csv(os.path.join(METH_FLD, f), sep='\t', index_col=0, skipinitialspace=True)
    df = df.dropna(how='all').drop('CpG_sites_hg19', axis=1)
    f_l.append(df)
cf = pd.concat(f_l, axis=0)
cf = cf.iloc[:, 1:].apply(lambda x: ((x * cf.avg_coverage).sum())/(cf.avg_coverage.sum()) ,axis=0) # weights don't add up to one/100, so need to divide by sum of weights
#
fld = METH_FLD
h = pd.read_excel(os.path.join(METH_FLD, 'DNA__Illumina_450K_methylation_Gene_average.xls'), index_col=0, skiprows=10)
h = h.loc[h['Chromosome f']!='Un'].drop(['Entrez gene id e', 'Chromosome f', 'Cytoband f'], axis=1)
h = h.iloc[:, 2:].apply(lambda x: ((x * (h['End f'] - h['Start f'])).sum())/((h['End f'] - h['Start f']).sum()) ,axis=0) # weights don't add up to one/100, so need to divide by sum of weights


