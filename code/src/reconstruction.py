import sys, os
home_path = '/home/tbarata'
sys.path.extend([os.path.join(home_path, pkg, 'src') for pkg in ('cobamp', 'troppo')])

import numpy as np
import json
import pandas as pd
import copy
import os
import cobra
import generic_model
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.core import OmicsContainer
from cobamp.utilities.parallel import batch_run
from multiprocessing import cpu_count
from cell_lines import CellLineIds
from pathos.multiprocessing import ProcessingPool as Pool
from itertools import repeat
from mewpy.simulation import get_simulator
from generic_model import GenericMdOperations
from cobra.io import save_matlab_model



class Reconstruction:
    '''
    - class to reconstruct models
    '''
    def __init__(self, model, cl_path):
        '''
        :param model: the generic model to work with
        :param cl_path: path json to file with cell lines to use
        '''
        self.model = model
        with open(cl_path, mode='r') as file:
            cl_lines = json.load(file)
        self.cl_lines = cl_lines
        self.nthreads = int(cpu_count() - (cpu_count() / 4))

    def preprocess_transcriptomics(self, data):
        '''
        - process transcriptomics data
        :param data: path to file with transcriptomics data to preprocess
        :return:
        '''
        data = pd.read_csv(data, sep = '\t')
        self.cl_lines = list(set(data.columns).intersection(set(self.cl_lines)))
        data['gene_id'] = data.gene_id.apply(lambda x: x.split('.')[0])
        data.set_index('gene_id', inplace=True)
        data = data.loc[set(data.index).intersection(set([g.id for g in self.model.genes])), self.cl_lines]
        # note: there are 30 genes in model that do not exist in the transcriptomics dataset
        # reactions of those genes may be protected:
        no_data_genes = set([g.id for g in self.model.genes]) - set(data.index)
        # save genes with 0 expression on all cell lines for later and remove them from the dataset - to avoid division by zero when calc. gene scores:
        all_zero_genes = data[np.sum(data == 0, axis=1) == data.shape[1]].index
        data.drop(data[np.sum(data == 0, axis=1) == data.shape[1]].index, axis=0, inplace=True)
        self.gene_data, self.no_data_genes, self.all_zero_genes = data, no_data_genes, all_zero_genes
        #return data, no_data_genes, all_zero_genes

    def get_gene_scores(self, ug_p, lg_p, lc_p, protect, gene_path):
        '''
        - calculates gene scores for certain parameters
        :param self.gene_data: dataframe with processed transcriptomics data (the output of 'preprocess_transcriptomics')
        :param self.no_data_genes: set with ids of genes that exist in the generic model but for which there isn't transcriptomics information
        :param self.all_zero_genes: set with ids of genes that have zero expression across all conditions
        :param ug_p: percentile used to obtain the upper global threshold
        :param lg_p: percentile used to obtain the lower global threshold
        :param lc_p: percentile used to obtain the local threshold. if False, the mean is used instead
        :param protect: whether to protect genes for which there isn't gene score information. either True or False
        :param gene_path: path to file where to save gene scores
        '''
        # determine gene thresholds:
        if not lc_p:
            lc_thr = self.gene_data.apply(lambda x: np.mean(x), axis=1)
        elif isinstance(lc_p, int) or isinstance(lc_p, float):
            lc_thr = self.gene_data.apply(lambda x: np.percentile(x, lc_p), axis=1)
        ug_thr = np.percentile(self.gene_data, ug_p)
        lg_thr = np.percentile(self.gene_data, lg_p)
        thr = copy.deepcopy(lc_thr)
        thr[thr >= ug_thr] = ug_thr
        thr[thr <= lg_thr] = lg_thr
        # get gene scores:
        gene_sc = (self.gene_data.T / thr).T  # same as data.div(thr, axis=0)
        gene_sc = 5 * np.log(1 + gene_sc)
        # whether to protect genes for which there aren't gene score information or to give them score corresponding to 0 expression:
        no_data_dct = dict()
        if protect:
            for k in self.gene_data.columns: no_data_dct[k] = {g: np.max(gene_sc[k]) for g in self.no_data_genes}
        else:
            for k in self.gene_data.columns: no_data_dct[k] = {g: 5 * np.log(1) for g in self.no_data_genes}
        # add the genes that have no expression in any of the conditions, by giving them a score corresponding to zero expression:
        no_exp = pd.DataFrame.from_dict({k: {g: 5 * np.log(1) for g in self.all_zero_genes} for k in self.gene_data.columns})
        gene_sc = pd.concat([gene_sc, pd.DataFrame.from_dict(no_data_dct), no_exp])
        self.gene_scores = gene_sc
        gene_sc.to_csv(gene_path, sep='\t')

    def rc_scr_one(self, score_tuple, params):
        '''
        - gets reaction scores for gene scores of one model
        - funct used in 'get_reaction_scores'
        :param score_tuple: (('and' function, 'or' function), dict with gene ids as keys and gene score as values)
        :param params: {'rw': rw}
        :return: reaction scores for gene scores of one model
        '''
        aofx, data_dict = score_tuple
        # aofx is tupple with and/or rule functions
        # data_dict is dict with gene id vs gene score of a combination to test
        # create omics container:
        oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x', data=data_dict, nomenclature='custom')
        rw = params['rw']  # load parameters
        mr = rw.model_reader
        rcsc = oc_sample.get_integrated_data_map(model_reader=mr, and_func=aofx[0], or_func=aofx[1])
        return rcsc

    def get_reaction_scores(self, int_str, protect, path, medium_path, **kwargs):
        '''
        - gets reaction scores for gene scores of different models / cell lines / group of parameters
        :param self.gene_scores: gene scores dataframe obtained from 'get_gene_scores'
        :param self.model: generic model
        :param int_str: integration strategy to convert gene scores to reaction scores using the gene-reaction (AND/OR) rules
                        (function_for_AND, function_for_OR)
        :param protect: boolean that indicates whether to protect (give highest score) to reactions that do not have gene reaction rule (like some transport reactions)
        :param path: path to where to save the result
        :param medium_path: path to file with medium composition reactions
        :param smp_info: path to file with cell line info
        :param kwargs: numbthreads - number of threads default value is object attribute
        '''
        if 'numbthreads' in kwargs:
            nthreads = kwargs['numbthreads']
        else:
            nthreads = self.nthreads
        fld = '/'.join(path.split('/')[:-1])
        if not os.path.exists(fld):
            os.makedirs(fld)
        runs = {cb: (int_str, sc_d) for cb, sc_d in self.gene_scores.to_dict().items()}
        labs, iters = tuple(runs.keys()), tuple(runs.values())
        rw = ReconstructionWrapper(self.model, ttg_ratio=9999)
        # apply reaction scores to all models:
        output = batch_run(self.rc_scr_one, iters, {'rw': rw}, threads=min(len(runs), nthreads))
        lstout = list()
        for out in output:  # out corresponds to one models / cell lines / group of parameters
            sc_d = dict()
            for rc, sc in out.get_scores().items():  # for each reaction id and reaction score
                sc_d[rc] = sc
            lstout.append(sc_d)
        batch_rc_res = dict(zip(labs, lstout))
        batch_rc_res = pd.DataFrame.from_dict(batch_rc_res)
        # always protect exchange reactions of medium composition and give score=threshold to the remaining exchange reactions:
        with open(medium_path, 'r') as f:
            envcond = list(json.load(f).values())
        constr_ex_ids = [ex.id for ex in self.model.exchanges if ex.id in envcond]
        no_constr_ex_ids = [ex.id for ex in self.model.exchanges if ex.id not in envcond]
        spont_rc = [r.id for r in self.model.reactions if not r.gene_reaction_rule]
        batch_rc_res.drop(spont_rc, axis=0, inplace=True)  # drop spontaneous reactions cause they are NAs and they're added below
        constr_ex_dc = pd.DataFrame.from_dict({k: {rid: np.max(batch_rc_res[k]) for rid in constr_ex_ids} for k in batch_rc_res.columns})
        no_constr_ex_dc = pd.DataFrame.from_dict({k: {rid: 5 * np.log(2) for rid in no_constr_ex_ids} for k in batch_rc_res.columns})
        # whether to protect the spontaneous reactions that are Not exchange reactions:
        spont_no_ex_ids = set(spont_rc) - {ex.id for ex in self.model.exchanges}
        if protect:
            spont_frm = pd.DataFrame.from_dict({k: {rid: np.max(batch_rc_res[k]) for rid in spont_no_ex_ids} for k in batch_rc_res.columns})
        else:
            spont_frm = pd.DataFrame.from_dict({k: {rid: 5 * np.log(2) for rid in spont_no_ex_ids} for k in batch_rc_res.columns})
        final = pd.concat([batch_rc_res, constr_ex_dc, no_constr_ex_dc, spont_frm])
        final.to_csv(path, sep='\t')
        self.reaction_scores = final

    def weight_fc(self, rc_sc_path, alg, protect_rc, path_rc_wgts, essen_tsk_lst, **kwargs):
        '''
        - calculates reaction weights
        :param rc_sc_path: path to file with reaction scores
        :param alg: either 'init' or 'fastcore'
        :param protect_rc: list of reactions to protect - give maximum score for that condition/model. If False, nothing to protect.
        :param path_rc_wgts: path to where to save the result
        :param essen_tsk_lst: list of tasks that should always be protected - tasks essential for cell survival
                              with mandatory activity (the reactions needed for that task)
        :param kwargs: tsk_protect_path - path to file with tasks that should/should not be protected in each condition/cell line
                                          only use if user wants tasks to be protected
                       tsk_lst - list of tasks (all that can be executed in the generic model) with mandatory activity (the reactions needed for that task)
                                 only use if user wants tasks to be protected
        '''
        fld = '/'.join(path_rc_wgts.split('/')[:-1])
        if not os.path.exists(fld):
            os.makedirs(fld)
        data = pd.read_csv(rc_sc_path, sep='\t', index_col=0)
        # user specified reactions to protect - e.g. biomass:
        if protect_rc:
            for pr in protect_rc: data.loc[pr, :] = data.max()
        # reactions of essential tasks are always protected in all cell lines:
        essen_tsk_rc = {r for t in essen_tsk_lst for r in t.mandatory_activity}
        for r in essen_tsk_rc: data.loc[r, :] = data.max()
        # reactions of task consensus may be protected:
        if 'tsk_protect_path' in kwargs:
            # tsk_keep_df = pd.read_csv(tsk_protect_path, sep='\t', index_col=0)
            tsk_keep_df = pd.read_csv(kwargs['tsk_protect_path'], sep='\t', index_col=0)
            '''
            tsk_keep_df.drop(list(range(171,174)), inplace=True) # demethylation tasks are not anymore tissue specific, they are generic
            '''
            # tsk_dct = {t.name: t.mandatory_activity for t in tsk_lst}
            tsk_dct = {t.name: t.mandatory_activity for t in kwargs['tsk_lst']}
            for cond in tsk_keep_df:
                tsks_cnd = list(tsk_keep_df[tsk_keep_df[cond] == 1].index)
                rc_keep = {r for t in tsks_cnd for r in tsk_dct[str(t)]}
                data.loc[rc_keep, cond] = data[cond].max()
        if alg == 'init':
            data[data <= 5 * np.log(2)] = -8.0
        elif alg == 'fastcore':
            pass
        else:
            print('only fastcore and init are available')
        data.to_csv(path_rc_wgts, sep='\t')
        self.reaction_weights = data
        return data

    def apply_alg_one(self, rc_wt, params):
        '''
        - applies reconstruction algorithm to one model
        :param rc_wt: dictionary with reaction ids as keys and reaction weights as values
        :param params: 'rw' is reconstruction wrapper of generic model;
                       'alg' is algorithm to apply. either 'init' or 'fastcore'
                       'thr' is threshold to use with 'fastcore'. reactions above this thrshold will be considered core reactions.
        :return:
        '''
        rw = params['rw']
        alg = params['alg']
        # model = params['model']
        thr = params['thr']
        try:
            if alg == 'init':
                # print(rc_wt)
                # return rw.run(tINITProperties(rc_wt, solver='CPLEX'))
                # scores_indx = {model.reactions.index(k): v for k, v in scores.items()}
                scores = {k: v if (v is not None) else -8 for k, v in rc_wt.items()}
                return rw.run_from_omics(omics_data=scores, integration_strategy='', algorithm='tinit', solver='CPLEX')
            elif alg == 'fastcore':
                # core_indx = [model.reactions.index(rid) for rid in core]
                # return rw.run(FastcoreProperties(core_indx, flux_threshold=1e-7, solver='CPLEX'))
                core = [[k for k, v in rc_wt.items() if v > thr]]
                return rw.run_from_omics(omics_data=core, integration_strategy='', algorithm='fastcore', solver='CPLEX')
        except Exception as e:
            print(e)
            return {r: False for r in rw.model_reader.r_ids}

    def apply_alg_all(self, model, alg, path, **kwargs):
        '''
        - apply reconstruction algorithm to reaction weights in different conditions/models
        :param self.reaction_weights: dataframe with weights for reactions. rows are reactions and columns are conditions/models
        :param model: generic model to use
        :param alg: algorithm to use. either 'init' or 'fastcore'
        :param path: path to where to save the result
        :param kwargs: 'thr' is the threshold to use with 'fastcore'. reactions above this thrshold will be considered core reactions.
                       'nthreads' is the number of threads, default is self.nthreads
        '''
        fld = '/'.join(path.split('/')[:-1])
        if not os.path.exists(fld):
            os.makedirs(fld)
        if 'numbthreads' in kwargs:
            nthreads = kwargs['numbthreads']
        else:
            nthreads = self.nthreads
        rw = ReconstructionWrapper(model, ttg_ratio=9999)
        labs, iters = list(), list()
        for cond in self.reaction_weights:
            labs.append(cond)
            iters.append(self.reaction_weights[cond].to_dict())
        thr = kwargs['thr']
        if alg == 'init':
            output = batch_run(self.apply_alg_one, iters, {'rw': rw, 'alg': 'init', 'thr': thr}, threads=min(len(iters), 1))
        elif alg == 'fastcore':
            output = batch_run(self.apply_alg_one, iters, {'rw': rw, 'alg': 'fastcore', 'thr': thr}, threads=min(len(iters), nthreads))
        batch_res = dict(zip(labs, output))  # create dictionary with conditions' labels and dictionaries mapping reactionId to bolean (True - means reaction to keep after reconstruction)
        final = pd.DataFrame.from_dict(batch_res, orient='index')
        final.to_csv(path, sep='\t')
        self.algo_reactions = final

    @staticmethod
    def apply_exch_const(model, cond, med_ex_ids, **kwargs):
        '''
        - applies to model the environmental constraints of exchange reactions in constr_ex,
          opens with bounds -1000, 1000 exchange reactions in med_ex_ids
          closes the remaining exchange reactions
        :param model: model to use
        :param cond: condition or cell line. string equal to a column name of envcond
        :param med_ex_ids: exchange reactions of media to open from -1000 to 1000
        :param kwargs: - envcond: dataframe with bounds for some exchange reactions in different cell lines
                       - constr_ex: exchange reactions to which apply flux bound constraints
        :return: model with the mentioned changes
        '''
        if 'constr_ex' in kwargs:
            for ex in model.exchanges:
                if ex.id in kwargs['constr_ex']:
                    print(ex.id)
                    model.reactions.get_by_id(ex.id).bounds = (kwargs['envcond'].loc[ex.id, cond], kwargs['envcond'].loc[ex.id, cond + '.1'])
                elif ex.id in med_ex_ids:
                    model.reactions.get_by_id(ex.id).bounds = (-1000.0, 1000.0)
                else:
                    model.reactions.get_by_id(ex.id).bounds = (0.0, 1000.0)
        else:
            for ex in model.exchanges:
                if ex.id in med_ex_ids:
                    model.reactions.get_by_id(ex.id).bounds = (-1000.0, 1000.0)
                else:
                    model.reactions.get_by_id(ex.id).bounds = (0.0, 1000.0)
        return model

    def find_gapfill_ids(self, univ_md, incomp_md):
        '''
        - find ids of reactions that need to be added to gapfill using cobrapy gapfill
        :param univ_md: generic model
        :param incomp_md: reconstructed
        :return gapfill_sol_ids: ids of reactions that are needed to gapfill
        '''
        gapfiller = cobra.flux_analysis.gapfilling.GapFiller(incomp_md, universal=univ_md, exchange_reactions=False,
                                                             demand_reactions=False, integer_threshold=1e-9,
                                                             lower_bound=1e-6)  # lower_bound = 1e-6:  need for objective reaction flux to be positive
        gapfiller.model.solver.problem.parameters.timelimit.set(1200)  # gapfill can take a max of 20 min (1200 sec) otherwise it will continue without solution
        gapfiller.model.solver.configuration.tolerances.feasibility = 1e-9
        gapfiller.model.solver.configuration.tolerances.integrality = 1e-9
        gapfiller.model.solver.configuration.tolerances.optimality = 1e-9
        gapfill_sol = gapfiller.fill()[0]
        gapfill_sol_ids = [rc.id for rc in gapfill_sol]
        return gapfill_sol_ids

    def gapfill_md(self, generic_path, context_md, cond, med_ex_ids, rc_keep, **kwargs):
        '''
        - gapfill infeasible models
        :param generic_path: path to generic xml model
        :param context_md: path to context specific model
        :param cond: condition or cell line. string equal to a column name of envcond
        :param med_ex_ids: exchange reactions of media to open from -1000 to 1000
        :param rc_keep: pandas series with boolean result of reconstruction algorithm. True means the reaction should be kept in reconstruction. False means it shouldn't.
        :param kwargs - constr_ex: exchange reactions to which apply flux bound constraints
                      - envcond: dataframe with bounds to apply for some exchange reactions in different cell lines
        :return: gapfilled model
        '''
        # set generic model:
        md_generic = generic_model.GenericMdOperations(xml_path=generic_path).model  # creates another independent instance of generic and obtains the model
        md_generic.objective = context_md.objective
        obj_rid_lst = [el.to_json()['name'] for el in md_generic.objective.variables if 'reverse' not in el.to_json()['name']]
        for rid in obj_rid_lst:
                md_generic.reactions.get_by_id(rid).bounds = (1e-6, 1000.0)
        if 'constr_ex' in kwargs:
            Reconstruction.apply_exch_const(model=md_generic, cond=cond, med_ex_ids=med_ex_ids, envcond=kwargs['envcond'], constr_ex=kwargs['constr_ex'])
        else:
            Reconstruction.apply_exch_const(model=md_generic, cond=cond, med_ex_ids=med_ex_ids)
        # remove reactions in context model - the False reactions except if they are medium or constrained exchanges:
        rc_to_remove = list({r.id for r in context_md.reactions if not rc_keep[r.id]} - set(med_ex_ids) - constr_ex)
        context_md.remove_reactions(rc_to_remove)
        # gapfill:
        gap_ids = self.find_gapfill_ids(univ_md=md_generic, incomp_md=context_md)
        print(f'{cond} needs to gapfill with:', gap_ids)
        gap_rc = [md_generic.reactions.get_by_id(rid).copy() for rid in gap_ids]
        print('before gapfill', len(context_md.reactions))
        context_md.add_reactions(gap_rc)
        print('after gapfill', len(context_md.reactions))
        return context_md

    def build_rc_md(self, reconst_df, cond, m, med_ex_ids, generic_path, fba, **kwargs):
        '''
        - returns a flux distribution of context reconstructed model.
          the objective used in simulation is the generic model pre-defined objective
        - informs whether the model is able to grow or not
        - gapfills when model is infeasible,
          after assuming a positive value for a objective, regardless if that obj. is growth or other.
        - if gapfill does not work, then returns the name of the condition/cell line
        :param reconst_df: dataframe with True/Falses when reaction is to be kept/or not, in accordance to reconstruction algorithm.
                           columns are conditions or cell lines
        :param cond: condition or cell line.
                     string equal to a column name of envcond and of reconst_df.
        :param m: generic model
        :param med_ex_ids: exchange reactions of media to open from -1000 to 1000
        :param generic_path: path to generic xml model
        :param fba: boolean. if True means FBA is done. If False, pFBA is done instead.
        :param kwargs - envcond: dataframe with bounds to apply for some exchange reactions in different cell lines
                      - constr_ex: exchange reactions to which apply flux bound constraints
        :return: context model and flux distribution (dictionary reaction:flux )
        '''
        r_c = reconst_df[cond]
        with m as md:
            # build context specific model:
            for r in md.reactions:
               if not r_c[r.id] and r.id not in [ex.id for ex in md.exchanges]:
                   md.reactions.get_by_id(r.id).bounds = (0.0, 0.0)
            # constrain with medium + 3 metab of human1:
            if 'constr_ex' in kwargs:
                Reconstruction.apply_exch_const(model=md, cond=cond, med_ex_ids=med_ex_ids, envcond=kwargs['envcond'], constr_ex=kwargs['constr_ex'])
            else:
                Reconstruction.apply_exch_const(model=md, cond=cond, med_ex_ids=med_ex_ids)
            # simulate:
            # sets the objective reactions with a minimum positive flux:
            obj_rid_lst = [el.to_json()['name'] for el in md.objective.variables if 'reverse' not in el.to_json()['name']]
            for rid in obj_rid_lst:
                md.reactions.get_by_id(rid).bounds = (1e-6, 1000.0)
            try:
                if fba:
                    sol = md.optimize()
                else:
                    sol = cobra.flux_analysis.pfba(md)
                flx = sol.fluxes.to_dict()
                bf = flx['adaptbiomass']
                print(f'{cond} is growing: {bf}')
                return flx
            except:
                print(f' except {cond} not working . trying gapfill')
                if 'constr_ex' in kwargs:
                    md = self.gapfill_md(generic_path=generic_path, context_md=md, cond=cond, med_ex_ids=med_ex_ids, rc_keep=r_c, envcond=kwargs['envcond'], constr_ex=kwargs['constr_ex'])
                else:
                    md = self.gapfill_md(generic_path=generic_path, context_md=md, cond=cond, med_ex_ids=med_ex_ids, rc_keep=r_c)
                try:
                    if fba:
                        sol = md.optimize()
                    else:
                        sol = cobra.flux_analysis.pfba(md)
                    flx = sol.fluxes.to_dict()
                    bf = flx['adaptbiomass']
                    print(f'{cond} is growing with gapfill: {bf}')
                    return flx
                except:
                    print(f'{cond} not working with gapfill')
                    return cond

    def reconstruct_mds(self, algo_reactions_pth, medium_path, m, generic_path, obj_id, smp_info, flx_path, fba, **kwargs):
        '''
        - gets pFBA flux distributions for each condition/model using a user-defined objective
        - informs whether each model is able to grow or not, gapfills when models are infeasible
        :param self.algo_reactions_pth: dataframe with result from reconstruction algorithm.
                                    True are reactions to keep. columns are reactions and rows are conditions
        :param medium_path: path to .json file with medium composition, which is a dict with metabolite: reaction
        :param m: generic model
        :param generic_path: path to generic xml model
        :param obj_id: dict representing an objective (with one or more reactions), of type reaction_id:weight
        :param smp_info: path to file with cell line info
        :param flx_path: path to where to save the results
        :param fba: boolean. if True means FBA is done. If False, pFBA is done instead.
        :param kwargs - envcond: path to dataframe with bounds to apply for some exchange reactions in different cell lines
                      - constr_ex: exchange reactions to which apply flux bound constraints
        '''
        out_fld = '/'.join(flx_path.split('/')[:-1])
        if not os.path.exists(out_fld):
            os.makedirs(out_fld)
        # apply objective to generic model:
        m.objective = {m.reactions.get_by_id(rid): w for rid, w in obj_id.items()}
        # build reconstructed models:
        algo_reactions = pd.read_csv(algo_reactions_pth, sep='\t', index_col=0)
        reconst_df = algo_reactions.T
        with open(medium_path, mode='r') as f:
            hams_med_cmp = json.load(f)
        med_ex_ids = [ex.id for ex in m.exchanges if ex.id in hams_med_cmp.values()]  # exchange reactions that are part of media
        all_flx_d = dict()
        infeasbl = list()
        for cond in reconst_df:
            print(cond)
            # reconstruct context model:
            if 'constr_ex' in kwargs:
                envcond_df = pd.read_excel(kwargs['envcond'], sheet_name='ExpFlux_lb_ub', index_col=2).iloc[:, 2:]
                CellLineIds.convert_cell_line_names(file=envcond_df, smp_info=smp_info)
                envcond_df.dropna(inplace=True)
                flx_res = self.build_rc_md(reconst_df=reconst_df, cond=cond, m=m, med_ex_ids=med_ex_ids, generic_path=generic_path, fba=fba, envcond=envcond_df, constr_ex=kwargs['constr_ex'])
            else:
                flx_res = self.build_rc_md(reconst_df=reconst_df, cond=cond, m=m, med_ex_ids=med_ex_ids, generic_path=generic_path, fba=fba)
            if type(flx_res) == dict: # if a flux distribution with flux in the objective was generated (feasible)
                all_flx_d[cond] = flx_res
            else:
                infeasbl.append(cond)
        final = pd.DataFrame.from_dict(all_flx_d)
        final.to_csv(flx_path, sep='\t')
        if infeasbl:
            inf_flx_pth = flx_path.split('.tsv')[0] + '_infeasible_conditions.tsv'
            pd.DataFrame(infeasbl).to_csv(inf_flx_pth, sep='\t')

    @staticmethod
    def build_rc_md_for_batch(cond, obj_id, reconst_df, md_pth, med_ex_ids, fba, include_algo_res, md_info_pth, md_info_sheet, tss_spc_md_fld_path, algor, with_tsk, **kwargs):
        '''
        - same as 'build_rc_md' but to be used with processes
        - also allows to build context model based only on constraining some exchange fluxes (generic model is used),
          if include_algo_res is False
        :param obj_id: dict representing an objective (with one or more reactions), of type reaction_id:weight
        :param reconst_df: dataframe with True/Falses when reaction is to be kept/or not, in accordance to reconstruction algorithm.
                           columns are conditions or cell lines
        :param cond: condition or cell line.
                     string equal to a column name of envcond and of reconst_df.
        :param md_path: path to generic model
        :param med_ex_ids: exchange reactions of media to open from -1000 to 1000
        :param fba: boolean. if True means FBA is done. If False, pFBA is done instead.
        :param include_algo_res: if True means that a context specific model will be reconstructed based on reconst algorithm results and flux constraints of some exchange reactions
                                 if False means that only flux constraints of some exchange reactions will be used to build context specific models (the complete/generic model is used)
        :param md_info_pth: path to excel file with tissue specific DNA methylation information
        :param md_info_sheet: name of excel sheet with tissue specific DNA methylation information
        :param tss_spc_md_fld_path: path to folder where to save cell specific reconstructed models
        :param algor: algorithm name used
        :param with_tsk: bool, whether highest reaction score was given to necessary reaction of tissue specific tasks or not
        :param kwargs - constr_ex: exchange reactions to which apply flux bound constraints
                      - envcond_df: dataframe with bounds to apply for some exchange reactions in different cell lines
        :return: context model and flux distribution (dictionary reaction:flux )
        '''
        print(cond)
        r_c = reconst_df[cond]
        md = GenericMdOperations(xml_path=md_pth, fea_tol=1e-9).model
        # apply objective to generic model:
        md.objective = {md.reactions.get_by_id(rid): w for rid, w in obj_id.items()}
        if 'constr_ex' in kwargs:
            # constrain with medium + some metab exchange values:
            exch_const = {ex.id: (kwargs['envcond_df'].loc[ex.id, cond], kwargs['envcond_df'].loc[ex.id, cond + '.1']) if ex.id in kwargs['constr_ex'] else (-1000.0, 1000.0) if ex.id in med_ex_ids else (0.0, 1000.0) for ex in md.exchanges}
        else:
            exch_const = {ex.id: (-1000.0, 1000.0) if ex.id in med_ex_ids else (0.0, 1000.0) for ex in md.exchanges}
        # if we want to include algo reconstruction results (we could just apply flux const on exchanges of generic model instead):
        if include_algo_res:
            ex_ids = [ex.id for ex in md.exchanges]
            cup = {r.id: (0.0, 0.0) for r in md.reactions if (not r_c[r.id]) and (r.id not in ex_ids)}
        # remove reactions considered 'False' by the algorithm,
        # except those from medium and non-catalyzed reactions of DNA demethylation that
        # (because of being non-catalyzed and non-essential for the task) don't have a gene score and therefore are always removed during reconstruction:
        to_remove = list(set(cup.keys()) - set(['consdirectDNA5CaC', 'consdirectDNA5fC']))
        md.remove_reactions(to_remove, remove_orphans=True)
        # exclude models for which no DNA methylation reaction information is known:
        tss_spc_df = pd.read_excel(md_info_pth, sheet_name=md_info_sheet, skipfooter=11, index_col=0, skiprows=6)
        if np.isnan(tss_spc_df.loc['MAM01722n', cond]): # if no information on DNA methylation is given
            print(f' No DNA methylation info available for {cond}')
        else:
            # simulate:
            simul = get_simulator(md, envcond=exch_const)
            simul.objective = obj_id
            if fba:
                sol = simul.simulate(method='FBA')
            else:
                sol = simul.simulate(method='pFBA')
            flx = sol.fluxes
            if sol.status.value == 'Optimal' and flx['adaptbiomass'] > 1e-9 and flx['MAR08641'] > 1e-9:
                print(f'{cond} is growing: {flx["adaptbiomass"]}, and methylates: {flx["MAR08641"]} ')
                # save reconstructed model as .mat:
                if with_tsk:
                    fld = os.path.join(tss_spc_md_fld_path, algor, 'including_tsks')
                else:
                    fld = os.path.join(tss_spc_md_fld_path, algor)
                if not os.path.exists(fld):
                    os.makedirs(fld)
                if cond == '786O_KIDNEY':
                    cond = 'KIDNEY_786O' # because gecko model will be build with matlab script and no matalab variable can start with numbers
                md.id = cond
                save_matlab_model(md, os.path.join(fld, cond + '.mat'))
                return cond, flx
            else:
                print(f'{cond} not working')
        return cond

    @staticmethod
    def reconstruct_mds_batch(algo_reactions_pth, medium_path, generic_path, obj_id, flx_path, fba, include_algo_res, md_info_pth, md_info_sheet, tss_spc_md_fld_path, envcond, with_tsk, **kwargs):
        '''
        - does the same as 'reconstruct_mds' but using many processes
        - also allows to build context models based only on constraining some exchange fluxes (generic model is used),
          if include_algo_res is False
        :param algo_reactions_pth: dataframe with result from reconstruction algorithm.
                                    True are reactions to keep. columns are reactions and rows are conditions
        :param medium_path: path to .json file with medium composition, which is a dict with metabolite: reaction
        :param m: generic model
        :param obj_id: dict representing an objective (with one or more reactions), of type reaction_id:weight
        :param smp_info: path to file with cell line info
        :param flx_path: path to where to save the results
        :param fba: boolean. if True means FBA is done. If False, pFBA is done instead.
        :param include_algo_res: if True means that a context specific model will be reconstructed based on reconst algorithm results and flux constraints of some exchange reactions
                                 if False means that only flux constraints of some exchange reactions will be used to build context specific models (the complete/generic model is used)
        :param md_info_pth: path to excel file with tissue specific DNA methylation information
        :param md_info_sheet: name of excel sheet with tissue specific DNA methylation information
        :param tss_spc_md_fld_path: path to folder where to save cell specific reconstructed models
        :param envcond: path to dataframe with experimental bounds of growth rate to apply
        :param with_tsk: bool, whether highest reaction score was given to necessary reaction of tissue specific tasks or not
        :param kwargs - nb_processes: number of processes to use. if not provided 1/4 of the OS processes will be used
        '''
        out_fld = '/'.join(flx_path.split('/')[:-1])
        if not os.path.exists(out_fld):
            os.makedirs(out_fld)
        algor = flx_path.split('/')[-1].split('_')[1]
        # build reconstructed models:
        algo_reactions = pd.read_csv(algo_reactions_pth, sep='\t', index_col=0)
        reconst_df = algo_reactions.T
        with open(medium_path, mode='r') as f:
            hams_med_cmp = json.load(f)
        med_ex_ids = list(hams_med_cmp.values())  # exchange reactions that are part of media
        all_flx_d = dict()
        infeasbl = list()
        cond_lst = list(reconst_df.columns)
        if 'nb_processes' in kwargs:
            total_processes = kwargs['nb_processes']
        else:
            total_processes = int(cpu_count() - (cpu_count() / 4))
        nb_nodes = min(total_processes, len(cond_lst))
        print('number_of_nodes:', nb_nodes)
        pool = Pool(nodes=nb_nodes)
        flx_lst = pool.map(Reconstruction.build_rc_md_for_batch, cond_lst, repeat(obj_id), repeat(reconst_df), repeat(generic_path),repeat(med_ex_ids), repeat(fba), repeat(include_algo_res), repeat(md_info_pth), repeat(md_info_sheet), repeat(tss_spc_md_fld_path), repeat(algor), repeat(with_tsk))
        pool.close() # close pool
        pool.join() # master process waits for the worker processes to finish
        pool.clear()
        print('WE ARE HERE!!!')
        # reconstruct context model:
        for flx_tup in flx_lst:
            if type(flx_tup) == str: # if infeasible
                infeasbl.append(flx_tup)
            elif flx_tup is None:
                continue
            else: # if a flux distribution with flux in the objective was generated (feasible)
                all_flx_d[flx_tup[0]] = flx_tup[1]
        final = pd.DataFrame.from_dict(all_flx_d)
        if not final.empty:
            final.to_csv(flx_path, sep='\t')
        if infeasbl:
            inf_flx_pth = flx_path.split('.tsv')[0] + '_infeasible_conditions.tsv'
            pd.DataFrame(infeasbl).to_csv(inf_flx_pth, sep='\t')
        print('REACHED THE END OF RECONST MDS BATCH')


