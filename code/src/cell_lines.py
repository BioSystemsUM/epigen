import os.path

import pandas as pd
import re
import json

class CellLineIds:
    '''
    - class to deal with cell line ids
    '''

    @staticmethod
    def convert_cell_line_names(file, smp_info, **kwargs):
        '''
        - convert alias names of cell lines in columns of a provided file to official CCLE cell line names
        :param file: file with alias of cell line names in columns, the first column can have other information
        :param smp_info: path to file with cell line info
        :param kwargs no_fst_clmn: boolean. True if the first column does Not have other information besides cell line names
        :return: the same file with the official CCLE cell line names
        '''
        # repl_cl_nm = {
        #     'A549/ATCC': 'A549_LUNG',
        #     'MDAMB435': 'MDAMB435S_SKIN',
        #     'OVCAR3': 'NIHOVCAR3_OVARY',
        #     'HL60(TB)': 'HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE',
        #     'U251': 'U251MG_CENTRAL_NERVOUS_SYSTEM',
        #     'MDAMB231/ATCC': 'MDAMB231_BREAST',
        #     'SR': 'SR786_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'}
        repl_cl_nm = {
            'A549ATCC': 'A549_LUNG',
            'MDAMB435': 'MDAMB435S_SKIN',
            'OVCAR3': 'NIHOVCAR3_OVARY',
            'HL60TB': 'HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE',
            'U251': 'U251MG_CENTRAL_NERVOUS_SYSTEM',
            'MDAMB231ATCC': 'MDAMB231_BREAST',
            'SR': 'SR786_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE',
            '7860': '786O_KIDNEY'}
        cl_info = pd.read_csv(smp_info)
        alias = cl_info.copy()
        alias = alias[['CCLE_Name', 'alias']].dropna(subset=['CCLE_Name', 'alias'])
        alias.set_index('CCLE_Name', inplace=True)
        alias = alias.iloc[:, 0].to_dict()
        alias = {el: k for k, v in alias.items() if ',' in v for el in v.split(', ')}
        alias_a = {k[1:]: v for k, v in alias.items() if k.startswith(' ')}
        alias_a.update({k: v for k, v in alias.items() if not k.startswith(' ')})
        nms = cl_info[['CCLE_Name', 'cell_line_name']].dropna(subset=['cell_line_name', 'CCLE_Name'])
        nms.set_index('cell_line_name', inplace=True)
        nms = nms.iloc[:, 0].to_dict()
        strp = cl_info[['CCLE_Name', 'stripped_cell_line_name']].dropna(subset=['stripped_cell_line_name', 'CCLE_Name'])
        strp.set_index('stripped_cell_line_name', inplace=True)
        strp = strp.iloc[:, 0].to_dict()
        strp.update(alias_a)
        strp.update(nms)
        strp.update(repl_cl_nm)
        n_col = list()
        for nm in file.columns:
            nm = ''.join(re.split(r'-|\s|_|\/|\(|\)', nm))
            if nm[-2:] == '.1' and nm[:-2] in strp:
                n_col.append(strp[nm[:-2]]+'.1')
            elif nm in strp:
                n_col.append(strp[nm])
            elif nm[-2:] == '.1':
                n_col.append('not_found.1')
            else:
                n_col.append('not_found')
        if 'no_fst_clmn' in kwargs:
            if kwargs['no_fst_clmn']:
                file.columns = n_col
                if 'not_found' in file.columns:
                    file.drop(['not_found'], axis=1, inplace=True)
        else:
            n_col[0] = file.columns[0]
            file.columns = n_col
            if ('not_found' in file.columns) and ('not_found.1' in file.columns):
                file.drop(['not_found', 'not_found.1'], axis=1, inplace=True)


    def intersect_cell_lines(self, envcond, smp_info, methlt, transcriptomics):
        '''
        - finds cell lines common to gene expression, metabolomics and richelle environmental constraints
        :param smp_info: path to file with cell line info
        :param methlt: path to one of the files with DNA methylation data
        :param envcond: path to file with envrionmental constrainsts on exchange fluxes
        :param transcriptomics: path to file with transcriptomics data
        :return: cell lines common to all three datasets
        '''
        # env conditions according to Richelle/Zielenski:
        CellLineIds.convert_cell_line_names(file=envcond, smp_info=smp_info)
        richelle_cl = {cl for cl in envcond.columns[1:] if not cl.endswith('.1')}
        # trasnciptomics data:
        raw_exp = pd.read_csv(transcriptomics, sep='\t')
        raw_exp.drop('transcript_ids', axis=1, inplace=True)
        raw_exp.gene_id = raw_exp.gene_id.apply(lambda x: x.split('.')[0])
        raw_exp.set_index('gene_id', inplace=True)
        # methylation data:
        raw_meth = pd.read_csv(methlt, sep='\t')
        cl_mt = [el for el in raw_meth.columns[3:]]
        # intersection of cells lines at each dataset:
        cl_int = set(raw_exp.columns).intersection(set(cl_mt)).intersection(richelle_cl)
        return cl_int

    def find_cell_lines_to_use(self, smp_info, envcond, cl_path):
        '''
        - find common cell lines to available datasets and saves them
        :param smp_info: path to file with cell line info
        :param envcond: path to file with envrionmental constrainsts on exchange fluxes
        :param cl_path: path to cell lines to use (intersection of cell lines of all datasets)
        '''
        envcond_richelle = pd.read_excel(envcond, sheet_name='Constraint_Used', index_col=2).iloc[:, 2:]
        # find common cell lines to available datasets and update column cell line names by the standard CCLE cell line names:
        cl_lines = self.intersect_cell_lines(envcond=envcond_richelle, smp_info=smp_info)
        # save list of the common cell lines:
        fld = '/'.join(cl_path.split('/')[:-1])
        if not os.path.exists(fld):
            os.makedirs(fld)
        out_file = open(cl_path, mode='w')
        json.dump(list(cl_lines), out_file)
        out_file.close()


