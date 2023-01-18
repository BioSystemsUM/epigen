home_path = '/home/tbarata'

import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])
import pandas as pd
from scipy.io import savemat
import numpy as np
from generic_model import GenericMdOperations

TRANSCRIPTOMICS='data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.txt'
METH_MD_FLD='support/models'
MEDIUM_PATH='data/hams_medium_composition.json'
TRANSCRIPTOMICSMAT='support/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.mat'

raw_exp = pd.read_csv(TRANSCRIPTOMICS, sep='\t')
raw_exp.drop('transcript_ids', axis=1, inplace=True)
raw_exp.gene_id = raw_exp.gene_id.apply(lambda x: x.split('.')[0])
raw_exp.set_index('gene_id', inplace=True)
cellid = np.array(list(raw_exp.columns), dtype=object)
genes = np.array(list(raw_exp.index), dtype=object)
tpm = raw_exp.values
depmap = {'cellID': cellid, 'genes': genes, 'tpm': tpm}
transcdir = '/'.join(TRANSCRIPTOMICSMAT.split('/')[:-1])
if not os.path.exists(transcdir):
    os.makedirs(transcdir)
savemat(TRANSCRIPTOMICSMAT, depmap)

# convert dna methylation model to mat:
for m_nm in os.listdir(METH_MD_FLD):
    # m_nm=os.listdir(METH_MD_FLD)[0]
    if m_nm.endswith('.xml'):
        print(m_nm)
        sbml_path = f'{METH_MD_FLD}/{m_nm}'
        gen = GenericMdOperations(xml_path=sbml_path, fea_tol=1e-9)
        gen.convert_xml_to_mat(xml_path=sbml_path, closed=False, medium_path=MEDIUM_PATH)

