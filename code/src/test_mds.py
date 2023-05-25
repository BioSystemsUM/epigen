import sys, os
home_path = '/home/tbarata'
sys.path.extend([os.path.join(home_path, pkg, 'src') for pkg in ('cobamp', 'troppo')])

import pandas as pd
import numpy as np
from troppo.omics.core import IdentifierMapping, TypedOmicsMeasurementSet
from cobamp.wrappers.external_wrappers import get_model_reader
from cobamp.core.optimization import BatchOptimizer
from cobamp.utilities.parallel import MP_THREADS
import seaborn as sns
import matplotlib.pyplot as plt
import re
import json
import os
from multiprocessing import cpu_count
from reconstruction import Reconstruction
from sklearn.metrics import matthews_corrcoef
from scipy.stats import pearsonr
from cell_lines import CellLineIds
from scipy.stats.mstats import spearmanr

class TestMds:
    '''
    - class to make tests on models
    '''

    def __init__(self, generic_model):
        '''
        :param generic_model: generic model
        '''
        self.generic_model = generic_model
        self.generic_cobamp = get_model_reader(generic_model, ttg_ratio=100000).to_cobamp_cbm('CPLEX')
        self.generic_cobamp_genes = self.generic_cobamp.gpr.get_genes()
        self.nthreads = int(cpu_count() - (cpu_count() / 4))

    def prepare_achilles_set(self, achilles_path):
        '''
        - preprocess Achilles dataset with experimental lethality gene scores of different cell lines:
        :param achilles_path: path to the achilles dataset with CERES normalized score
        :return: - same dataset with ensembl gene ids instead of gene names and entrez ids
                 - genes not in metabolic model are dropped
        '''
        ach_df = pd.read_csv(achilles_path, index_col=0)  # achilles 20Q1
        hgnc_df = pd.read_csv('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt', index_col=0, sep='\t')  # gene id conversion dataframe
        hgnc_df['entrez_id'] = [str(int(x)) if not np.isnan(x) else np.nan for x in hgnc_df['entrez_id']]
        mapping = IdentifierMapping('human_transcriptomics', hgnc_df)
        # Create an omics measurement set object with cell lines' dataset components:
        ach_mset = TypedOmicsMeasurementSet(ach_df.index, ach_df.columns, ach_df.values, mapping)
        # keep entrez gene id only in column names:
        ensembl_patt = re.compile('\([0-9]*\)')
        ach_mset.column_names = [str(ensembl_patt.findall(k)[0].replace('(', '').replace(')', '')) for k in ach_mset.column_names]
        # convert entrez to ensembl:
        conv_dict = mapping.get_id_table(ach_mset.column_names, 'entrez_id').set_index('entrez_id')['ensembl_gene_id'].to_dict()  # dict of achiles_entrez_gene_ids: ensembl gene ids
        ach_mset.column_names = [conv_dict[k] if k in conv_dict else np.nan for k in ach_mset.column_names]  # keep only ensembl gene ids, when convertion from entrez didn't work it was put as nan - conversion
        # drops columns corresponding to genes not in metabolic model:
        ach_mset.drop(columns=ach_mset.data.columns[~ach_mset.data.columns.isin(self.generic_cobamp_genes)])  # '~' symbol turns Trues to Falses
        self.ach_mset = ach_mset
        return ach_mset

    def inact_react_gn_ko(self):
        '''
        - get reactions that become inactive when a gene is knockout
        :param self.generic_cobamp: generic model in cobamp format
        :param self.generic_cobamp_genes: all genes of cobamp_model
        :return kos_to_test: list of combination of reactions that are inactive, depending on which gene is knockout
        :return reaction_states_inactive: dict where key is gene to knockout and value is id of reaction(s) that are inactive in universal model when that knockout happens
        '''
        # func bellow creates dict where all model genes are True except gene we want to knockout which is False:
        def get_state_dict_from_ko(ko):
            d = {k: True for k in self.generic_cobamp.gpr.get_genes()}
            for k in ko:
                d[k] = False
            return d
        # create dict {gene_to_ko:{gene:bool}}. bool is True for all genes except for the knockout:
        state_dicts = {gko: get_state_dict_from_ko([gko]) for gko in self.generic_cobamp_genes}
        # create dict {gene_to_ko:{reaction:value}. 'reaction' is the reaction that is inactive or without GPR associated and value is respectively 'False' or 'None':
        reaction_states = {gko: {i: j for i, j in
                                 {k: self.generic_cobamp.gpr.eval_gpr(ind, state_dicts[gko]) for ind, k in
                                  enumerate(self.generic_cobamp.reaction_names)}.items() if not j} for gko in
                           state_dicts.keys()}
            # above, {k: cobamp_model.gpr.eval_gpr(ind, state_dicts[gko]) for ind, k in enumerate(cobamp_model.reaction_names)}
            # is a dict for a certain gene knockout of type {reaction: True/False/None}.  True/False/None is active/inactive/no GPR rule associated

            # {i: j for i, j in {k: cobamp_model.gpr.eval_gpr(ind, state_dicts[gko]) for ind, k in enumerate(cobamp_model.reaction_names)}.items() if not j}
            # is a dict for a certain gene knockout of type {reaction:False/None}. False/None is inactive/no GPR rule associated

        # create dict {gene_ko: reaction_id} where reaction_id is the id of the reaction(s) that are inactive in generic model for that gene_ko:
        reaction_states_inactive = {k: tuple(frozenset({i for i, j in v.items() if j is not None})) for k, v in reaction_states.items()}
        # create dict {inactive_reaction(s): gene(s)}
        # can be read as, reaction(s) is/are inactive when any of the genes in corresponding set is knockout:
        rko_sets = {}
        for k, v in reaction_states_inactive.items():
            if v not in rko_sets.keys():
                rko_sets[v] = {k}
            else:
                rko_sets[v] |= {k}  # |= in dict means update: if inactive reaction is already in rko_sets, joins the new gene to the set of other genes
        # exclude empty reactions, for when a gene knockout doesn't have inactive associated reaction:
        kos_to_test = [l for l in list(rko_sets.keys()) if len(l) > 0]
        # kos_to_test is a list of combination of reactions that are inactive, depending on which gene is knockout
        return kos_to_test, reaction_states_inactive

    def build_ctx_dict(self, rc_rc_keep, cond, envcond, constr_ex, med_ex_ids):
        '''
        - builds dictionary {reaction_id: (lb, ub)} for a certain condition/model based on reconstruction algorithm results
        and medium constraints to apply
        :param rc_rc_keep: dataframe with boolean indicating whether reaction is to be kept or not in accordance with reconstruction algorithm
                           columns are reactions and rows are conditions
        :param cond: condition/cell to take into account. has to be one of the rownames of rc_rc_keep
        :param envcond: dataframe with bounds to apply for some exchange reactions in different cell lines
        :param constr_ex: exchange reactions to which apply flux bound constraints
        :param med_ex_ids: exchange reactions of media to open from -1000 to 1000
        '''
        cond_keep = rc_rc_keep[cond]
        context_dict = dict()
        for r in self.generic_model.reactions:
            if not cond_keep[r.id] and r.id not in [ex.id for ex in self.generic_model.exchanges]:
                context_dict[r.id] = (0.0, 0.0)
            else:
                context_dict[r.id] = self.generic_model.reactions.get_by_id(r.id).bounds
        #####context_dict = {rid: (0.0, 0.0) if not v else (-1000.0, 1000.0) for rid, v in cond_keep.items()}
        # override the reaction bounds with exchange bounds adapted for medium composition and constraints:
        for ex in self.generic_model.exchanges:
            if ex.id in constr_ex:
                context_dict[ex.id] = (envcond.loc[ex.id, cond], envcond.loc[ex.id, cond + '.1'])
            elif ex.id in med_ex_ids:
                context_dict[ex.id] = (-1000.0, 1000.0)
            else:
                context_dict[ex.id] = (0.0, 1000.0)
        return context_dict

    def batch_ko_optimize(self, kos, context, objective, objective_sense, nthreads):
        '''
        - for a context model, create a list with scores,
          each score representing if a group of reactions inactivated by a gene ko is essential (0) or affects biomass (<1) or is non-essential (>=1)
          position in the list corresponds to position in list 'kos'
        :param kos: list of groups of reactions that are inactive, depending on which gene is knockout.
        :param context: dictionary {reaction_id: (lb, ub)} for a certain condition/model based on reconstruction algorithm results
                        and medium constraints to apply
        :param objective: objective to apply to cobamp reconstructed models
        :param objective_sense: when False means to maximize the objective
        :param nthreads: number of threads
        '''
        cb = self.generic_cobamp
        with cb as cb_gm:
            # reconstruct a context cobamp model and do FBA:
            cb_gm.set_reaction_bounds('adaptbiomass', lb=1e-6, ub=1000.0)
            for rid, bds in context.items():
                cb_gm.set_reaction_bounds(rid, lb=bds[0], ub=bds[1])  # change reaction bounds
            sol = cb_gm.optimize(objective, objective_sense)
            ofv = sol.objective_value()
            stat = sol.status()
            print(ofv, stat)
            if sol.status() == 'optimal' and ofv > 0:
                # get list of dicts where each dict corresponds to set of reactions ko when a gene is ko
                # [{{index_reactionA, index_reactionB, ...}: (0,0)}, ...]:
                gko_to_rko = [{cb_gm.decode_index(rko, 'reaction'): (0, 0) for rko in ko_set} for ko_set in kos]  # ko_set is comb of reactions that are ko when a gene is ko
                # simulate several models for the same reconstructed model (1 condition) where one gene is ko at each time:
                bopt = BatchOptimizer(cb_gm.model, threads=int(min(len(gko_to_rko), nthreads, MP_THREADS)))
                opt_result = bopt.batch_optimize(gko_to_rko, [{self.generic_cobamp.decode_index(rko, 'reaction'): v for rko, v in objective.items()}] * len(gko_to_rko), [objective_sense] * len(gko_to_rko))
                    # objective = {'biomass_human': 1}, '1' is weight of objective
                    # objective_sense is False, which means to optimize
                    # opt_result - list of 'optimize' solutions - for each of reaction/reaction combinations
                return [k.objective_value() / ofv if (k.status() == 'optimal') and ofv > 0 else 0 for k in opt_result]  # if infeasible, model can't get a steady-state distribution, so ko is lethal
            else:
                print('Reconstructed model objective is 0... skipping')
                return []

    def sim_exp_leth(self, kos_to_test, reaction_states_inactive, ach_mset, rc_rc_keep, envcond, constr_ex,  med_ex_ids, ach_id_dct, **kwargs):
        '''
        - get dict {condition: dataframe} where dataframe has simulated and experimental lethality scores on different columns and genes are rows
        :param kos_to_test: list of groups of reactions that are inactive, depending on which gene is knockout.
        :param reaction_states_inactive: dict where key is gene to knockout and value is id of reaction(s) that are inactive in universal model when that knockout happens
        :param ach_mset: preprocessed Achilles dataset with experimental lethality gene scores of different cell lines
        :param rc_rc_keep: dataframe with boolean indicating whether reaction is to be kept or not in accordance with reconstruction algorithm
                           columns are reactions and rows are conditions
        :param envcond: dataframe with bounds to apply for some exchange reactions in different cell lines
        :param constr_ex: exchange reactions to which apply flux bound constraints
        :param med_ex_ids: exchange reactions of media to open from -1000 to 1000
        :param ach_id_dct: dictionary {CCLE_cell_line_name: DepMap_id}
        :param kwargs nthreads: number of threads
        :return:
        '''
        if 'nthreads' in kwargs:
            nthreads = kwargs['nthreads']
        else:
            nthreads = self.nthreads
        corr_coef_dict = dict()
        rc_rc_keep = rc_rc_keep.T
        for cond in rc_rc_keep:
            context_dict = self.build_ctx_dict(rc_rc_keep=rc_rc_keep, cond=cond, envcond=envcond, constr_ex=constr_ex, med_ex_ids=med_ex_ids)
            # get dict {{reactionA, reactionB, ...}: score} where
            # with scores representing if the group of inactive reactions (corresponding to a gene ko)
            # was essential (0) or affected biomass (<1) or not (>=1).
            bkopt = self.batch_ko_optimize(kos=kos_to_test, context=context_dict, objective={'adaptbiomass': 1}, objective_sense=False, nthreads=nthreads)
            opt_res = dict(zip(kos_to_test, bkopt))
            # create dict where key is gene ko and value is score representing if was essential (0) or affected biomass (<1) or not (>=1):
            of_res = {k: opt_res[v] if v in opt_res.keys() else 1 for k, v in reaction_states_inactive.items()}
                # if else condition  guarantees that genes
                # which don't have associated inactive reaction when ko get a score of 1 - meaning they don't affect biomass
            # create dataframe where:
            # 1st col is experimental score centered on zero, where negative values are lethal genes
            # 2nd col is score calculated from reconstructed model
            # where score 0 / ]0,1[ / >=1 are essential genes/ genes that affected biomass / genes not decreasing biomass:
            ddf = ach_mset.data.reindex(index=[ach_id_dct[cond]], columns=of_res).append(pd.Series(of_res, name='simulated')).T
            # create dict {condition: value} where value is a dataframe with experimental and simulation scores:
            corr_coef_dict[cond] = ddf.dropna()  # removes rows with NAN in Achilles dataset
            # exclude and save to a list empty conditions (those for which there isn't Achilles lethality experimental scores):
            corr_coef_final = dict()
            empty = list()
            for k, v in corr_coef_dict.items():
                if not len(v):
                    empty.append(k)
                else:
                    corr_coef_final[k] = v
        return corr_coef_final, empty

    def mcc_dif_lth_thrs(self, corr_coef_final, ach_id_dct, final_path, test_lth_thrs, **kwargs):
        '''
        - get series with Mathews Correlation Coefficient (MCC) score for each condition, and the median of those values,
          using the lethality thresholds (experimental and simulated) that give the highest median of MCC (across all conditions)
        :param corr_coef_final: dict of type {cond: dataframe}.
               dataframe is lethality scores where columns are experimental or simulated, and rows are genes' ensembl ids.
        :param ach_id_dct: dict {CCLE_cell_line: Achiles_cell_line_id}
        :param final_path: path where to save the result
        :param test_lth_thrs: whether to test (or not) several lethality thresholds to see which is better. boolean.
        :param kwargs - el: user provided threshold for Achilles experimental lethality score
                      - sl: user provided threshold for simulated lethality score
        :return md_lth: lethality experimental and simulated thresholds with the highest median of MCC across different conditions
        '''
        corr_coef_params = {}
        if test_lth_thrs:
            # determine correlation values for dif. lethality thresholds:
            # get dict  {cond:{k:v}} where k is number representing a score threshold to consider gene as lethal and v is correlation coef:
            for cond, ddf in corr_coef_final.items():
                corr_coef_params[cond] = {}
                for el in list(np.arange(-0.1, -1, -0.1)):  # thresholds for Achilles experimental lethality score
                    for sl in [0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.999]:  # thresholds for simulated lethality score
                        corr_coef_params[cond][(round(el, 1), sl)] = matthews_corrcoef((ddf['simulated'] < sl), (ddf[ach_id_dct[cond]] < el))
            # get dataframe of correlation coef where column is condition and row is comb of score cutoffs to consider gene lethal:
            scores_df = pd.DataFrame.from_dict(corr_coef_params)
            # get the lethality experimental and simulated thresholds with the highest median of MCC across different conditions.
            # cannot use idxmax, in case there is a draw. upon a draw, just arbitrarily select the one with lowest thresholds:
            mdan = list(scores_df.median(axis=1)[scores_df.median(axis=1) == scores_df.median(axis=1).max()].index)
            sl = np.min([e[1] for e in [e for e in mdan if e[0] == np.min([e[0] for e in mdan])]])
            md_lth = [e for e in [e for e in mdan if e[0] == np.min([e[0] for e in mdan])] if e[1] == sl][0]
            # get the MCC for each condition - using experimental and simulated lethality threshold with the highest MCC median:
            fsc = scores_df.loc[md_lth, :]
            fsc_md = pd.Series([fsc.median()], index=['median'])
            # also get the median MCC score value across all conditions (using as lethality thresholds those with highest MCC median):
            rs = pd.concat([fsc, fsc_md])
            rs.to_csv(final_path)
            return md_lth
        elif not test_lth_thrs:
            # use provided lethality thresholds:
            # get dict  {cond:v} where cond is the condition/cell line, and v is correlation coef:
            for cond, ddf in corr_coef_final.items():
                corr_coef_params[cond] = [matthews_corrcoef((ddf['simulated'] < kwargs['sl']), (ddf[ach_id_dct[cond]] < kwargs['el']))]
            # get series with MCC score for different conditions/cell lines - lethality thresholds are the ones provided by user:
            fsc = pd.DataFrame.from_dict(corr_coef_params)
            # also get the median MCC score value across all conditions:
            fsc_md = fsc.median(axis=1)
            fsc_md.index = ['median']
            rs = pd.concat([fsc.T, fsc_md], axis=0)
            rs.to_csv(final_path)

    def gene_lethal(self, envcond_path, smp_info, medium_path, achilles_path, rc_rc_path,
                               final_path, constr_ex, test_lth_thrs, **kwargs):
        '''
        - get the Mathews Correlation Coeficient (MCC) between experimental and simulated lethal genes for each condition
        - gives the median MCC score value across all conditions
        - returns lethality experimental and simulated thresholds with the highest median of MCC across different conditions
        :param envcond_path: path to dataframe with bounds to apply for some exchange reactions in different cell lines
        :param smp_info: path to file with cell line info
        :param medium_path: path to .json file with medium composition. contains a dict of the format {metabolite_id: reaction_id}
        :param self.generic_model: generic model
        :param achilles_path: path to file with gene essentiality data - achilles dataset
        :param rc_rc_path: path to file with reactions to keep/ or not during model reconstruction
        :param final_path: path to file where to save the results
        :param constr_ex: exchange reactions to which apply flux bound constraints
        :param test_lth_thrs: whether to test (or not) several lethality thresholds to see which gives better median MCC. boolean.
        :param kwargs - el: user provided threshold for Achilles experimental lethality score
                      - sl: user provided threshold for simulated lethality score
        :return: saves results to final_path and returns lethality experimental and simulated thresholds with the highest median of MCC across different conditions
        '''
        # load datasets needed:
        rc_rc_keep = pd.read_csv(rc_rc_path, sep='\t', index_col=0)
        envcond_df = pd.read_excel(envcond_path, sheet_name='Constraint_Used', index_col=2).iloc[:, 2:]
        CellLineIds.convert_cell_line_names(file=envcond_df, smp_info=smp_info)
        with open(medium_path, mode='r') as f:
            hams_med_cmp = json.load(f)
        med_ex_ids = [ex.id for ex in self.generic_model.exchanges if ex.id in hams_med_cmp.values()]  # exchange reactions that are part of media
        ach_id_df = pd.read_csv(smp_info)
        ach_id_df.set_index('CCLE_Name', inplace=True)
        ach_id_dct = ach_id_df['DepMap_ID'].to_dict()
        # preprocess dataframe with experimental lethality scores:
        ach_mset = self.prepare_achilles_set(achilles_path=achilles_path)
        # get inactive reactions when a gene is knockout:
        kos_to_test, reaction_states_inactive = self.inact_react_gn_ko()
        # get dict {condition: dataframe} where dataframe has simulated and experimental lethality scores on different columns and genes are rows:
        corr_coef_final, empty = self.sim_exp_leth(kos_to_test=kos_to_test,
                                                  reaction_states_inactive=reaction_states_inactive, ach_mset=ach_mset,
                                                  rc_rc_keep=rc_rc_keep, envcond=envcond_df, constr_ex=constr_ex,
                                                  med_ex_ids=med_ex_ids, ach_id_dct=ach_id_dct)
        if test_lth_thrs:
            mcc_dif = self.mcc_dif_lth_thrs(corr_coef_final, ach_id_dct, final_path, test_lth_thrs)
            return mcc_dif
        elif not test_lth_thrs:
            self.mcc_dif_lth_thrs(corr_coef_final, ach_id_dct, final_path, test_lth_thrs, sl=kwargs['sl'], el=kwargs['el'])

    @classmethod
    def create_file_list(cls, path, incl_wrd):
        '''
            - create list with paths to files which names include a provided word
            :param path: base folder
            :param incl_wrd: word that the names of files we want to list must include
        '''
        file_lst = list()
        for root, fldrs, files in os.walk(path):
            for f in files:
                if (incl_wrd in f) and (('.png' not in f) and ('.pdf' not in f) and ('.svg' not in f)):
                    file_lst.append(os.path.join(root, f))
        return file_lst

    @classmethod
    def compare_leth_scores(cls, base_fld):
        '''
        - gives the best data processing choices (e.g. protect tasks, constraint certain exchange reac., etc)
          that give the best gene lethality MCC scores
        :param base_fld: root folder where to look recursively for files with gene lethality MCC scores
        :return max_median: list of processing choices that gave the maximum median MCC score (across all conditions/cell lines)
        :return cond_more_mx_score: list of processing choices that ar most frequent, among those that gave the best MCC score for each cell line/condition
        '''
        # create list of paths to files storing MCC score for gene lethality results:
        file_lst = cls.create_file_list(path=base_fld, incl_wrd='lethality')
        df_lst = list()
        # build a dataframe with gene lethality MCC scores for all strategies tested:
        for file_pth in file_lst:
            df = pd.read_csv(file_pth, index_col=0)
            fn = '__'.join(file_pth.split('/')[2:-1] + [file_pth.split('/')[-1].split('.')[0]])
            df.rename(columns={'0': fn}, inplace=True)
            df_lst.append(df)
        fdf = pd.concat(df_lst, axis=1)
        # get the strategy(ies) that gave the maximum median MCC score (across all conditions/cell lines):
        max_median = list(fdf.median()[fdf.median() == fdf.median().max()].index)
            # when more than one value is max idxmax() selects the first..so not using this
        # check which strategy is most frequent, among those that gave the best MCC score for each cell line/condition:
        fdf_max = [el for l in fdf.apply(lambda r: list(r[r == r.max()].index), axis=1) for el in l]
        cd = dict()
        vlst = list()
        for el in set(fdf):
            cd[el] = fdf_max.count(el)
            vlst.append(fdf_max.count(el))
        cond_more_mx_score = [k for k, v in cd.items() if v == max(vlst)]
        return max_median, cond_more_mx_score

    @classmethod
    def compare_growth(cls, base_fld, gr_path, smp_info):
        '''
        - gives the best data processing choices (e.g. protect tasks, constraint certain exchange reac., etc)
          that give the lowest relative error between experimentally determined and simulated growth rates
        - and those that give the highest spearman correlation with p-value below 0.05
        :param base_fld: root folder where to look recursively for files with fluxes
        :param gr_path: path to file with experimental growth rates
        :param smp_info: path to file with cell line info
        :return min_median: list of processing choices that gave the lowest relative error median between exp and siml growth rates (across all conditions/cell lines)
        :return cond_more_min_error: list of processing choices that are most frequent, among those that gave the lowest relative errors for each cell line/condition
        :return max_corr: list of processing choices that gave the highest spearman correlations between exp and siml growth rates (across all conditions/cell lines)
        '''
        # the 11 "good" cell lines from human1:
        eleven_cl = ['HT29_LARGE_INTESTINE', 'HOP62_LUNG', '786O_KIDNEY', 'MALME3M_SKIN', 'UO31_KIDNEY',
                     'HS578T_BREAST', 'HOP92_LUNG', 'NCIH226_LUNG', 'MDAMB231_BREAST',
                     'RPMI8226_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'SR786_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE']

        # create list of paths to files storing pFBA fluxes:
        file_lst = cls.create_file_list(base_fld, 'fluxes')
        # load file with experimental growth rates:
        grf = pd.read_excel(gr_path, sheet_name='GrowthRate', index_col=0).T
        CellLineIds.convert_cell_line_names(file=grf, smp_info=smp_info, no_fst_clmn=True)

        # - build a dataframe with relative errors between simulated and experimentally measured growth rates
        #   for each cell line/condition (rows) and each data processing decision (columns)
        # - build a dataframe with spearman correlation and p-value between simulated and experimentally measured growth rates
        #   across all cell lines/conditions, for each data processing decision:
        grf = grf.mean()
        df_lst = list()
        corr_lst = list()
        for file_pth in file_lst:
            fn = '__'.join(file_pth.split('/')[2:-1] + [file_pth.split('/')[-1].split('.')[0]])
            save_id = file_pth.split('.')[0]
            df = pd.read_csv(file_pth, index_col=0, sep='\t')
            df_bm = df.loc['adaptbiomass', :]
            grff = grf.loc[grf.index.isin(list(df_bm.index))]
            grff.name = 'experimental_biomass'
            df_bm = df_bm.reindex(index=grff.index)
            rs = pd.concat([df_bm, grff], axis=1)
            relative_error = np.abs(rs['adaptbiomass'] - rs['experimental_biomass']) / rs['experimental_biomass']
            relative_error.name = fn
            df_lst.append(relative_error)
            # also test correlation - use spearman method instead of pearson, cause the distribution is not normal:
                # sns.displot(df_bm)
                # plt.show()
            corrv_spear, p = spearmanr(df_bm, grff)
            corr_s = pd.Series([corrv_spear, p])
            corr_s.name = fn
            corr_lst.append(corr_s)
            # show dispersion graphs - all cell lines:
            plt_f = pd.DataFrame({'simulated': np.log10(df_bm + 1), 'experimental': np.log10(grff.astype(float) + 1)})
            sns.scatterplot(x='experimental', y='simulated', data=plt_f, s=2)
            sns.regplot(x=plt_f['experimental'], y=plt_f['simulated'])
            # plt.show()
            plt.savefig(save_id + '_growthcorr_all_cl.png')
            plt.close()
            # just the "good" seven:
            plt_f_seven = plt_f.loc[plt_f.index[plt_f.index.isin(eleven_cl)], :]  # in fact we use 7 instead of 11 cell lines, cause for some of those there are no methylation data or exp growth rates or other essential data
            sns.scatterplot(x='experimental', y='simulated', data=plt_f_seven, s=2)
            sns.regplot(x=plt_f_seven['experimental'], y=plt_f_seven['simulated'])
            # plt.show()
            plt.savefig(save_id + '_growthcorr_seven_cl.png')
            plt.close()
        fdf = pd.concat(df_lst, axis=1)
        corr_f = pd.concat(corr_lst, axis=1)
        corr_f.index = ['correlation', 'p-value']
        # get the strategy(ies) that gave the minimum median relative error (across all conditions/cell lines):
        min_median = list(fdf.median()[fdf.median()==fdf.median().min()].index)
        # get the strategy(ies) that gave the highest correlation with p-value below 0.05:
        corr_sb = corr_f.loc[:, corr_f.loc['p-value', :] < 0.05]
        max_corr = list(corr_sb.loc['correlation', :][corr_sb.loc['correlation', :] == corr_sb.loc['correlation', :].max()].index)
        # check which strategy(ies) is/are most frequent, among those that gave the lowest relative errors, for each cell line/condition:
        fdf_min = [el for l in fdf.apply(lambda r: list(r[r==r.min()].index), axis=1) for el in l]
        cd = dict()
        vlst = list()
        for el in set(fdf_min):
            cd[el] = fdf_min.count(el)
            vlst.append(fdf_min.count(el))
        cond_more_min_error = [k for k, v in cd.items() if v == max(vlst)]
        # verify which gives lower errors, using all cell lines or just the 11 "good" ones from human1:
        eleven_cl = ['HT29_LARGE_INTESTINE', 'HOP62_LUNG', '786O_KIDNEY', 'MALME3M_SKIN', 'UO31_KIDNEY', 'HS578T_BREAST', 'HOP92_LUNG', 'NCIH226_LUNG', 'MDAMB231_BREAST', 'RPMI8226_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'SR786_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', ]
        seven_cl = fdf.index[fdf.index.isin(eleven_cl)]  # in fact we use 7 instead of 11 cell lines, cause for some of those there are no methylation data or exp growth rates or other essential data
        min_median_val_seven_cl = fdf.loc[seven_cl, :].median().min()
        min_median_val_all_cl = fdf.median().min()
        if min_median_val_seven_cl < min_median_val_all_cl:
            print(f'Error is lower when using only the "good" cell : {min_median_val_seven_cl} < {min_median_val_all_cl}')
        elif min_median_val_all_cl < min_median_val_seven_cl :
            print(f'Error is lower when using all cell lines: {min_median_val_all_cl} < {min_median_val_seven_cl}')
        return min_median, cond_more_min_error, max_corr

    @classmethod
    def compare_exofluxes(cls, base_fld, env_path, smp_info, medium_path, constr_ex):
        '''

        :param base_fld: root folder where to look recursively for files with fluxes
        :param env_path: path to file with experimental fluxes for exometabolites
        :param smp_info: path to file with cell line info
        :param medium path: path to .json file with medium composition
        :param constr_ex: set with ids of exchange reactions which fluxes are constrained
        :return:
        '''
        # load medium composition:
        with open(medium_path, 'r') as f:
            med_comp = json.load(f)
        # create list of paths to files storing pFBA fluxes:
        file_lst = cls.create_file_list(base_fld, 'fluxes')
        # load file with fluxes calculated from experimental data:
        envcond_df = pd.read_excel(env_path, sheet_name='Exometabolomics_ALL', index_col=2).iloc[:, 2:]
        CellLineIds.convert_cell_line_names(file=envcond_df, smp_info=smp_info)
        # the 11 "good" cell lines from human1:
        eleven_cl = ['HT29_LARGE_INTESTINE', 'HOP62_LUNG', '786O_KIDNEY', 'MALME3M_SKIN', 'UO31_KIDNEY',
                     'HS578T_BREAST', 'HOP92_LUNG', 'NCIH226_LUNG', 'MDAMB231_BREAST',
                     'RPMI8226_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'SR786_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE']
        # exclude fluxes that were constrained for the simulation and 'human1_metabolite' column:
        if 'no exchange' in envcond_df.index:
            t_remove = set(constr_ex).union({'no exchange'})
        else:
            t_remove = set(constr_ex)
        envcond_df = envcond_df.drop(list(t_remove), axis=0).drop('human1_metabolite', axis=1)
        envcond_df = envcond_df.drop('adaptbiomass', axis=0)
        envf = envcond_df
        # build dataframe with correlations between simulated and experimental fluxes
        # for each reaction (rows) and each data processing decision (columns):
        corr_d = dict()
        corr_d_nine = dict()
        for file_pth in file_lst:
            # file_pth = file_lst[0]
            fn = '__'.join(file_pth.split('/')[2:-1] + [file_pth.split('/')[-1].split('.')[0]])
            save_id = file_pth.split('.')[0]
            df = pd.read_csv(file_pth, index_col=0, sep='\t')
            new_df = df.loc[envf.index, :]
            new_envf = envf.loc[:, new_df.columns]
            nine_df = new_df[new_df.columns[new_df.columns.isin(eleven_cl)]]  # in fact we use 9 instead of 11 cell lines, cause for some of those there are no methylation data or other essential data
            nine_envf = new_envf[new_envf.columns[new_envf.columns.isin(eleven_cl)]]
            def row_corr(simdf, expdf):
                lst = list()
                for (_, row_df), (_, row_envf) in zip(simdf.iterrows(), expdf.iterrows()):
                    corrv_spear, p = spearmanr(row_df, row_envf)
                    lst.append((corrv_spear, p))
                s = pd.Series(lst, index=simdf.index).apply(lambda x: 0 if x[1] > 0.05 else x[0])
                return s
            corr_d[fn] = row_corr(simdf=new_df, expdf=new_envf)
            corr_d_nine[fn] = row_corr(simdf=nine_df, expdf=nine_envf)
            # also get graph to visualize correlations - all cell lines:
            plt_f = pd.DataFrame({'simulated': np.log10(np.abs(new_df.to_numpy().flatten() + 1)), 'experimental': np.log10(np.abs(new_envf.to_numpy().flatten() + 1))})
            plt_f = pd.DataFrame({'simulated': new_df.to_numpy().flatten(), 'experimental': new_envf.to_numpy().flatten()})
                # np.abs() cause log of negative values gives NaN
            sns.scatterplot(x='experimental', y='simulated', data=plt_f, s=2)
            sns.regplot(x=plt_f['experimental'], y=plt_f['simulated'])
            #plt.show()
            plt.savefig(save_id + '_flxcorr_all_cl.png')
            plt.close()
            #  get graph to visualize correlations - 7 good cell lines:
            plt_f_nine = pd.DataFrame({'simulated': np.log10(np.abs(nine_df.to_numpy().flatten() + 1)), 'experimental': np.log10(np.abs(nine_envf.to_numpy().flatten() + 1))}) # in fact we use 9 instead of 11 cell lines, cause for some of those there are no methylation data or other essential data
            sns.scatterplot(x='experimental', y='simulated', data=plt_f_nine, s=2)
            sns.regplot(x=plt_f_nine['experimental'], y=plt_f_nine['simulated'])
            # plt.show()
            plt.savefig(save_id + '_flxcorr_nine_cl.png')
            plt.close()
        # correlations for all cell lines:
        corr_df = pd.DataFrame.from_dict(corr_d)
        # correlations for the severn "good" cell lines:
        corr_df_nine = pd.DataFrame.from_dict(corr_d_nine)
        return corr_df.median(), corr_df_nine.median()

    @classmethod
    def load_meth_dataset(cls, methlt_pth):
        '''
        - pre-process experimental DNA methylation dataset in txt format (NAs were replaced by the median value of that genome position across all cell lines)
          and then get a total value of DNA methylation for each cell line
        :param methlt_pth: path to DNA methylation dataset
        '''
        raw_df = pd.read_csv(methlt_pth, sep='\t', index_col=0, na_values=['    NaN', '     NA']).iloc[:-1, :]
        meth_df = raw_df.iloc[:, 2:]
        median_mth = meth_df.apply(lambda row: row.dropna().median(), axis=1).to_dict()
        meth_df.apply(lambda row: row.fillna(value=median_mth[row.name], inplace=True), axis=1)
        meth_res = meth_df.sum()
        return meth_res

    @classmethod
    def load_meth_dataset_gene_centric(cls, methlt_pth, smp_info):
        '''
        - pre-process experimental DNA methylation dataset in xls format for genes (NAs were replaced by the median value of that gene all cell lines)
          and then get a total value of DNA methylation for each cell line
        :param methlt_pth: path to DNA methylation dataset
        :param smp_info: path to file with cell line info
        '''
        meth_df = pd.read_excel(methlt_pth, sheet_name='Results', index_col=1, skiprows=10, skipfooter=7, na_values=['-']).iloc[:, 5:]
        median_mth = meth_df.apply(lambda row: row.dropna().median(), axis=1).to_dict()
        meth_df.apply(lambda row: row.fillna(value=median_mth[row.name], inplace=True), axis=1)
        meth_df.columns = [c.split(':')[1] for c in meth_df.columns]
        CellLineIds.convert_cell_line_names(file=meth_df, smp_info=smp_info, no_fst_clmn=True)
        meth_res = meth_df.sum()
        return meth_res

    @classmethod
    def create_corr_files(cls, file_pth, mth_res_lst, list_meth_rc, divide_by_biomass):
        '''
        - create scatter plots with experimental vs simulated fluxes for DNA methylation related reactions
        :param file_pth: path to file with simulated fluxes
        :param mth_res_lst: list of pandas series with experimental methylation data,
                            where each element corresponds to a methylation dataset
        :param list_meth_rc: list with ids of reactions related with DNA methylation
        :param divide_by_biomass: whether divide simulated flux values by biomass simulated values
        '''
        def corrfunc(x, y, **kwargs):
            tkeep = ~np.logical_or(np.isnan(x), np.isnan(y))
            x = x[tkeep]
            y = y[tkeep]
            nmod = tkeep.sum()
            pearson, pp = pearsonr(x, y)
            spearman, ps = spearmanr(x, y)
            ax = plt.gca()
            ax.annotate(f'pearson = {round(pearson, 2)}',
                        xy=(1.05, .9), xycoords=ax.transAxes, fontsize=12)
            ax.annotate(f'p-val = {round(pp, 2)}',
                        xy=(1.05, .8), xycoords=ax.transAxes, fontsize=12)
            ax.annotate(f'spearman = {round(spearman, 2)}',
                        xy=(1.05, .7), xycoords=ax.transAxes, fontsize=12)
            ax.annotate(f'p-val = {round(ps, 2)}',
                        xy=(1.05, .6), xycoords=ax.transAxes, fontsize=12)
            # ax.annotate(f'# cell lines = {nmod}',
            #             xy=(1.05, .5), xycoords=ax.transAxes, fontsize=12)

        sv_fld = '/'.join(file_pth.split('/')[:-1])
        nm = 'fluxcomp'
        if divide_by_biomass:
            list_meth_rc.append('adaptbiomass')
        simul_df = pd.read_csv(file_pth, index_col=0, sep='\t')
        simul_df = simul_df.loc[list_meth_rc, :]
        col_lst = [[col for col in simul_df.columns if col in exp_df.index] for exp_df in mth_res_lst]
        exp_df_lst = [exp_df[lst] for exp_df, lst in zip(mth_res_lst, col_lst)]
        df = pd.concat([simul_df.T] + exp_df_lst, axis=1)
        cols = [col.replace('arm_', '') if 'arm_' in col else re.sub('No[0-9]+', '', col) if 'No' in col else col for col in df]
        df.columns = cols
        df.dropna(inplace=True) # removes cell lines that do not have experimental values for all experimental dataframes
        d = {'Illumina_450K_methylation_Gene_average': 'Genes',
             'DNA_methylation_ratio': 'Global DNA methylation',
             'TSS1kb': 'Upstream of TSS',
             'cgi_CpG_clusters': 'CpG islands (clusters)',
             'enh_CpG_clusters': 'Enhancers (clusters)',
             'tss_CpG_clusters': 'TSS (clusters)'}
        df.rename(columns=d, inplace=True)
        df[df==0] = 1e-15 # 0 is replaced by small value because of log calculation
        df_log = np.log10(df[~df.isna()])
        exp_vars = [d[exp.name] for exp in exp_df_lst]
        # sim_vars = list(set(df.columns) - set(exp_vars) - {'adaptbiomass'})
        sim_vars = [cl for cl in df.columns if (cl not in exp_vars) and (cl != 'adaptbiomass')]
        if divide_by_biomass:
            df[sim_vars] = df[sim_vars].div(df['adaptbiomass'], axis=0)
        # correlation plots:
        fn_fld = os.path.join(sv_fld, nm)
        if not os.path.exists(fn_fld):
           os.makedirs(fn_fld)
        xl = [exp_vars[:int(len(exp_vars)/2)], exp_vars[int(len(exp_vars)/2):]]
        yl = [sim_vars[:int(len(sim_vars)/2)], sim_vars[int(len(sim_vars)/2):]]
        for x in xl:
            for y in yl:
                nmi = x[0].split(' ')[0] + '*' + y[0].split(' ')[0]
                fn_pth = os.path.join(fn_fld, nmi + '_corr.svg')
                g = sns.PairGrid(df_log, x_vars=x, y_vars=y, aspect=1.4)
                g.map(plt.scatter)
                g.map(corrfunc)
                for ax in g.axes.flat:
                    ax.xaxis.label.set_size(14)
                    ax.yaxis.label.set_size(13.5)
                    if np.min(df_log).min() == np.log10(1e-15):
                        ax.set_ylim([-15.1, ax.get_ylim()[1]])
                # for ax in g.axes.flat:
                #     ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
                plt.tight_layout()
                plt.savefig(fn_pth, dpi=300)
                plt.close()
        # histograms log10:
        variables = list(df_log.columns)
        plt.figure(figsize=(17, 10))
        plt.subplots_adjust(hspace=0.5)
        for i, var in enumerate(variables):
            ax = plt.subplot(len(variables)/4, 4, i + 1)
            n, bins, edges = ax.hist(df_log[var].dropna())
            ax.set_title(var, fontsize=26,pad=20)
            # ax.set_xticks(np.round(bins[:len(bins):2], 2))
            ax.set_xticks(np.round(bins[:len(bins):2], 2))
            ax.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True)
            ax.xaxis.set_tick_params(labelsize=11.5)
            ax.yaxis.set_tick_params(labelsize=11.5)
        plt.tight_layout()
        pth_hist = os.path.join(fn_fld, sim_vars[0] + '_hlog.svg')
        plt.savefig(pth_hist, dpi=300)
        plt.close()
        # histogram NOT log10:
        variables = list(df.columns)
        plt.figure(figsize=(17, 10))
        plt.subplots_adjust(hspace=0.5)
        for i, var in enumerate(variables):
            ax = plt.subplot(len(variables)/4, 4, i + 1)
            n, bins, edges = ax.hist(df[var].dropna())
            ax.set_title(var, fontsize=26, pad=20)
            ax.set_xticks(bins[:len(bins):2])
            ax.ticklabel_format(axis='x', scilimits=(0,0), useMathText=True)
            ax.xaxis.set_tick_params(labelsize=10.5)
            ax.yaxis.set_tick_params(labelsize=11.5)
        plt.tight_layout()
        pth_hist = os.path.join(fn_fld, sim_vars[0] + '_h.svg')
        plt.savefig(pth_hist, dpi=300)
        plt.close()

    @classmethod
    def obtain_fluxes(cls, algo, rc_keep_pth, md_pth, mth_corr_fld_path,
                            flux_str, constr_ex, medium_path, envcond,
                            smp_info, obj_lst, simul, include_algo_res):
        '''
        - obtain simulated fluxes for a list of objectives and experimental fluxes of methylation related reactions
        :param algo: name of reconstruction algorithm used (to get the recation to keep)
        :param rc_keep_pth: path to file with reactions to keep from reconstruction algorithm
        :param md_pth: path to file with generic model to use
        :param cl_pth: path json to file with cell lines to use
        :param mth_corr_fld_path: path to folder where to store fluxes and methylation correlation results
        :param flux_str: string representing name of file where fluxes are stored.
                         {algo}, {obj} and {const} are to be replaced by algorithm, objective and constraints
        :param constr_ex: exchange reactions to which apply flux bound constraints
        :param medium_path: path to .json file with medium composition, which is a dict with metabolite: reaction
        :param envcond: dataframe with bounds to apply for some exchange reactions in different cell lines
        :param smp_info: path to file with cell line info
        :param obj_lst: list with objectives to test
        :param simul: string with either pFBA or FBA, which is used in simulation
        :param include_algo_res: if True means that a context specific model will be reconstructed based on reconst algorithm results and flux constraints of some exchange reactions
                                 if False means that only flux constraints of some exchange reactions will be used to build context specific models (the complete/generic model is used)
        '''
        # load a model:
        model_fld_nm = md_pth.split('/')[-1]
        # for each onjective in a list of objectives
        # create a flux distribution - pfba:
        if simul == 'pFBA':
            fba = False
        elif simul == 'FBA':
            fba = True
        for obj_id in obj_lst:
            print('obj_id', obj_id)
            obj_id_name = {k: '{:.0e}'.format(v) for k, v in obj_id.items()}
            if include_algo_res:
                flx_path = os.path.join(mth_corr_fld_path, model_fld_nm, flux_str.format(algo=algo, obj=obj_id_name, const=constr_ex, simul=simul))
            else:
                flx_path = os.path.join(mth_corr_fld_path, model_fld_nm, flux_str.format(algo=algo, obj=obj_id_name, const=constr_ex, simul=simul).split('.')[0]+'_generic_only.tsv')
            Reconstruction.reconstruct_mds_batch(algo_reactions_pth=rc_keep_pth, medium_path=medium_path, envcond=envcond,
                                generic_path=md_pth, obj_id=obj_id, smp_info=smp_info, flx_path=flx_path,
                                constr_ex=constr_ex, fba=fba, include_algo_res=include_algo_res)

    @classmethod
    def obtain_methylation(cls, methlt_dt_fld, smp_info):
        '''
        - obtain list of experimental fluxes of methylation related reactions
        :param methlt_dt_fld: folder where DNA methylation data is
        :param smp_info: path to file with cell line info
        :return: list of experimental fluxes of methylation related reactions
        '''
        # create a list where each item is an experimental dataset
        # with the total values of DNA methylation of each cell line:
        # in this study there are 4 datasets of this type
        # - RRBS_ehn_CpG_clusters: DNA methylation clusters
        # - RRBS_cgi_CpG_clusters: DNA methylation CpG islands
        # - RRBS_ts_CpG_clusters DNA methylation Transcription Site CpG clusters
        # - RRBS_TSS1kb: DNA methylation promoter 1kb upstream TSS
        mth_res_lst = list()
        for file_nm in os.listdir(methlt_dt_fld):
            # file_nm = os.listdir(methlt_dt_fld)[6]
            # obtain total value of DNA methylation for each cell line:
            if file_nm.endswith('.txt'):
                mth_df = cls.load_meth_dataset(methlt_pth=os.path.join(methlt_dt_fld, file_nm))
                mth_df.name = '_'.join(file_nm.split('.')[0].split('_')[2:-1])
                mth_res_lst.append(mth_df)
            # elif file_nm.endswith('.xls'):
            elif file_nm.startswith('DNA__Illumina'):
                mth_df = cls.load_meth_dataset_gene_centric(methlt_pth=os.path.join(methlt_dt_fld, file_nm), smp_info=smp_info)
                mth_df.name = '_'.join(file_nm.split('.')[0].split('_')[2:])
                mth_res_lst.append(mth_df)
            elif file_nm.startswith('total_DNA'):
                methlt_pth = os.path.join(methlt_dt_fld, file_nm)
                mth_df = pd.read_excel(methlt_pth, header=0).iloc[0, :]
                mth_df.name = '_'.join(file_nm.split('.')[0].split('_')[1:])
                mth_res_lst.append(mth_df)
        return mth_res_lst