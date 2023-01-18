import os
import cobra
import pandas as pd
from cobra import Metabolite, Reaction
import json
from cobra.io import load_yaml_model, save_matlab_model
import re
from itertools import product

class GenericMdOperations:
    '''
    - class to deal with generic models
    :param yml_pth: path to file with .yml model
    '''
    yml_pth = '/home/tbarata/epigen/data/models/Human-GEM.yml'
    def __init__(self, xml_path, **kwargs):
        '''
        :param file_path: path where .xml model exist or where to save the .xml model
        '''
        self.xml_path = xml_path
        self.kwargs = kwargs
        if not os.path.exists(self.xml_path):
            self.create_xml()
        else:
            md = cobra.io.read_sbml_model(xml_path)
            if 'opt_tol' in kwargs:
                md.solver.configuration.tolerances.optimality = kwargs['opt_tol'] # 1e-09
            if 'fea_tol' in kwargs:
                md.solver.configuration.tolerances.feasibility = kwargs['fea_tol'] # 1e-09
            # md.solver.configuration.tolerances.integrality = 1e-09
            self.model = md

    def create_xml(self):
        '''
        - load human1.xx version develop branch in format .yaml and save in xml, then load in xml format
        :param xml_path: path to save and then load generic xml model
        '''
        #md_link = 'https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/develop/model/Human-GEM.yml'
        #modelPath, _ = urlretrieve(md_link)
        #m = cobra.io.load_yaml_model(modelPath)
        m = load_yaml_model(GenericMdOperations.yml_pth)
        cobra.io.write_sbml_model(m, self.xml_path)
        md = cobra.io.read_sbml_model(self.xml_path)
        if 'opt_tol' in self.kwargs:
            md.solver.configuration.tolerances.optimality = self.kwargs['opt_tol']  # 1e-09
        if 'fea_tol' in self.kwargs:
            md.solver.configuration.tolerances.feasibility = self.kwargs['fea_tol']  # 1e-09
        # md.solver.configuration.tolerances.integrality = 1e-09
        self.model = md

    def processed(self, model, meth_rc):
        '''
        - removes blocked reactions and associated orphans (orphanated genes, metabolites)
        - checks if there is cell growth and model is feasible upon positive flux values of biomass
        - confirms that no (de)methylation reaction is blocked upon positive flux values of biomass
        :param model: model
        :param meth_rc: list of reaction ids associated with DNA (de)methylation
        :return:
        '''
        # set model objective as 'adaptbiomass':
        model.objective = 'adaptbiomass'
        model.reactions.get_by_id('adaptbiomass').bounds = (0.0, 1000)
        # check whether (de)methylation reactions are blocked:
        blck_rc = cobra.flux_analysis.find_blocked_reactions(model)
        if not [id for id in meth_rc if id in blck_rc]:
            print('NO (de)methylation reactions are blocked')
            # remove all blocked reactions:
            model.remove_reactions(blck_rc, remove_orphans=True)
            model.solver = 'cplex'
            # check cell growth/feasibility:
            sol = model.optimize()
            obj_val = sol.objective_value
            if obj_val > 1e-7 and sol.status == 'optimal': # optimality tolreance is 1e-7 (while feasibility tolerance  is 1e-9)
                print(f'model is feasible and grows with obj_val: {obj_val}')
                return model


    @staticmethod
    def set_metab(mtlst, model, rid):
        '''
        - add metabolites described in metabolite string
        :param df: string with metabolites of a reactions with the form 'DNA5CaCn:-1.0, MAM02039n:-1.0, MAM01721n:1.0, MAM01596n:1.0'
        :param model: model with the reaction to which we add metabolites
        :param rid: id of reaction to which we add metabolites
        '''
        mt_d = dict()
        for el in mtlst.split(', '):
            k, v = el.split(':')[0], el.split(':')[1]
            mt_d[k] = float(v)
        r = model.reactions.get_by_id(rid)
        r.add_metabolites(mt_d)

    def create_generic_meth_md(self, extra_path, dna_meth_rc_path, meth_md_fld):
        '''
        - creates and saves xml model of the generic human1 model with extra reactions needed for DNA (de)methylation
        - also adapts gene rules with complexes of isozymes
        :param extra_path: path to file with info on reactions related to DNA methylation that are to be added to generic model
        :param dna_meth_rc_path: path to .json file with/where to save ids of DNA (de)methylation reactions common to all generic methylation models
        :param meth_md_fld: path to folder where to save the generic methylation models
        '''
        # get human1.xx model
        m = self.model

        # correct/adapt gene rules with complexes of isozymes:
        self.replace_gene_rule()

        # remove reaction of DNA methylation in the cytoplasm,
        # there is no DNA in cytoplasm for DNA to be methylated there:
        tst_md = m.copy()
        tst_md.remove_reactions(['MAR04072'], remove_orphans=True)

        # close the generic biomass reaction:
        tst_md.reactions.get_by_id('MAR13082').bounds = (0.0, 0.0)

        # update gene rule of DNA methylation reaction in the nucleus (explanation on 'gene_rules' tab of file /data/md_modify.py):
        rc_modify = pd.read_excel(extra_path, sheet_name='rc_modify')
        tst_md.reactions.get_by_id('MAR08641').gene_reaction_rule = rc_modify.loc[0, 'updated_gene_rule']

        # create new metabolites:
        mtb_include = pd.read_excel(extra_path, sheet_name='mtb_include', usecols=list(range(7)), skipfooter=4)
        mtb_include.apply(lambda row: tst_md.add_metabolites([Metabolite(id=row['metabolite_id'], name=row['name'], formula=row['formula'], charge=row['charge'], compartment='n')]), axis=1)

        # create new reactions:
        rc_include = pd.read_excel(extra_path, sheet_name='rc_include', usecols=[0, 4, 5, 8, 10], skipfooter=2)
        rc_include.fillna(value='', inplace=True)
        rc_include.loc[rc_include['gene_rule'] == 'uncatalyzed', 'gene_rule'] = ''
        rc_include.set_index('reaction_id', inplace=True)
        rc_include.reversible = rc_include.reversible.apply(lambda x: bool(x))
        # rc_common = rc_include.drop(rc_include.index[rc_include.index.str.contains('prodDNAtot')], axis=0)
        rc_common = rc_include
        rc_common.apply(lambda row: tst_md.add_reactions([Reaction(id=row.name, subsystem=row['subsystem'], lower_bound=0.0, upper_bound=1000.0)]), axis=1)
        for rid in rc_common.index[rc_common.reversible]: tst_md.reactions.get_by_id(rid).lower_bound = -1000.0
        for rid in rc_common.index: tst_md.reactions.get_by_id(rid).gene_reaction_rule = rc_common.loc[rid, 'gene_rule']

        # save DNA (de)methylation ids:
        common_ids = list(rc_common.index)
        with open(dna_meth_rc_path, mode='w') as f:
            json.dump(common_ids, f)

        # add metabolites and gene_rules to created reactions:
        for rid in rc_common.index:
            mtlst = rc_common.loc[rid, 'add_metabolites']
            GenericMdOperations.set_metab(mtlst=mtlst, model=tst_md, rid=rid)
            tst_md.gene_reaction_rule = rc_common.loc[rid, 'gene_rule']

        if not os.path.exists(meth_md_fld):
            os.mkdir(meth_md_fld)
        # next remove blocked reactions and test models ability to produce biomass:
        tst_md.objective = 'adaptbiomass'
        common_ids.append('MAR08641')
        prc_md = self.processed(model=tst_md, meth_rc=common_ids)
        cobra.io.write_sbml_model(prc_md, os.path.join(meth_md_fld, 'prodDNAtot')+'.xml')

        # create models, each model containing a different reaction for the formation of DNAtot - with different metabolites:
        # for rid in rc_include.index[0:6]:
        #     print(rid)
        #     mtlst = rc_include.loc[rid, 'add_metabolites']
        #     with tst_md as ct_md:
        #         ct_md.add_reactions([Reaction(id='prodDNAtot', subsystem=rc_include.loc[rid, 'subsystem'], lower_bound=0.0, upper_bound=1000.0)])
        #         GenericMdOperations.set_metab(mtlst=mtlst, model=ct_md, rid='prodDNAtot')
        #         ct_md.objective = 'adaptbiomass'
        #         # next remove blocked reactions and test models ability to produce biomass:
        #         prc_md = self.processed(ct_md)
        #         cobra.io.write_sbml_model(prc_md, os.path.join(meth_md_fld, rid)+'.xml')

    def convert_xml_to_mat(self, xml_path, closed, **kwargs):
        '''
        - save xml model as .mat model
        :param xml_path: path to xml model
        :param closed: bool whether to close (close uptake of non-medium metabolites) the model or not
        :param kwargs - medium_path: path to .json file with medium composition. necessary when 'closed' is True
        '''
        print(xml_path)
        md = self.model
        if closed:
            with open(kwargs['medium_path'], 'r') as f:
                med = list(json.load(f).values())
            for ex in md.exchanges:
                if ex.id in med:
                    ex.bounds = (-1000.0, 1000.0)
                else:
                    ex.bounds = (0.0, 1000.0)
            # blck_rc = cobra.flux_analysis.find_blocked_reactions(md)
            # md.remove_reactions(blck_rc, remove_orphans=True)
        mat_path = xml_path.split('.')[0]
        save_matlab_model(md, mat_path + '.mat')

    def replace_gene_rule(self):
        '''
        - changes gene rules with complexes of isozymes
         *  "((G1 and G2) or (G2 and G3)) and ((G5 and G6) or (G7 and G8))" --> "(G1 and G2 and G5 and G6) or (G1 and G2 and G7 and G8) or (G2 and G3 and G5 and G6) or (G2 and G3 and G7 and G8)"
         *  "(G1 or G2) and (G3 or G4)" --> "(G1 and G3) or (G1 and G4) or (G2 and G3) or (G2 and G4)"
        :param self.model: model to adapt gene rules
        :return: changes gene rules in place - in the model
        '''
        print('replace_gene_rule is adapted for gene rules of human 1.12 only . \n For other models change gene rules manually.')
        for r in self.model.reactions:
            gr = r.gene_reaction_rule
            if (')) and ((' in gr) and [p for p in gr.split(')) and ((') if ' or ' in p]:
                print(f"{r.id} contains '(( or )) and (( or ))' - complexes of isozymes")
                lst = [st.split(') or (') for st in [re.sub(r'\(\(|\)\)', '', el) for el in re.split('\)\) and \(\(|\)\) and | and \(\(', gr)]]
                ngr = '(' + ') or ('.join([' and '.join(i) for i in product(*lst)]) + ')'
                r.gene_reaction_rule = ngr
                print(r.id, r.gene_reaction_rule)
            # if 'MAR07161' == r.id or 'MAR07162' == r.id:
            elif ((') and ' in gr) or (') and (' in gr)) and (
                    [p for p in gr.split(') and ') if ' or ' in p] or [p for p in gr.split(') and (') if ' or ' in p]):
                print(f"{r.id} contains '( or ) and ( or )' - complexes of isozymes")
                lst = [st.split(' or ') for st in [re.sub(r'\(|\)', '', el) for el in re.split('\) and \(|\) and | and \(', gr)]]
                ngr = '(' + ') or ('.join([' and '.join(i) for i in product(*lst)]) + ')'
                r.gene_reaction_rule = ngr
                print(r.id, r.gene_reaction_rule)
