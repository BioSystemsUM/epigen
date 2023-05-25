home_path = '/home/tbarata'

import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])

import copy
from cobra.io import load_matlab_model, save_matlab_model
from mewpy.simulation import get_simulator
import json
import pandas as pd
import numpy as np
from functools import reduce

def add_necessary_rc(fld_md_pth, final_fld_md_pth, gen_md_base_pth, consen_tsk_pth, base_fld, md_info_pth, md_info_sheet, medium_path, with_tsk, demeth_tsk=False):
    '''
    - if required, adds essential reactions for cell-type specific tasks when those tasks are suppose to be done by the cell line
    - if required, adds essential reactions for the DNA demethylation tasks only when those tasks are suppose to be done by the cell line
    - adds reactions involved in DNA demethylation that do not require enzymes and that wouldn't be otherwise included
      (because they do not have an associated gene score and aren't essential for demethylation tasks)
    - replaces 'prodDNAtot' reaction by an equivalent one reflecting the cell line-specific proportion of DNA methylation
    - test if models are feasible and able to produce positive flux of biomass and DNA methylation
    :param fld_md_pth: path to folder with models to use
    :param final_fld_md_pth: path to folder where to save resulting .m models
    :param gen_md_base_pth: path to generic model
    :param consen_tsk_pth: path to .json file with non-essential tasks done by each cell line (and reactions needed for those tasks)
    :param base_fld: folder with file 'consen_tasks_to_keep.tsv' which in turn contains the tissue specific tasks that should be done by each cell line
    :param md_info_pth: path to excel file with tissue specific DNA methylation information
    :param md_info_sheet: name of excel sheet with tissue specific DNA methylation information
    :param medium_path: path to .json file with medium composition. contains a dict of the format {metabolite_id: reaction_id}
    :param with_tsk: bool, whether necessary reactions of tissue specific tasks were included or not
    :param demeth_tsk: bool, whether to include necessary reactions of demethylation tasks
                       (included only if those are done in specific cell line)
                       when the reactions needed for other cell specific tasks are NOT included.
                       default is False.
    '''
    generic_m = load_matlab_model(gen_md_base_pth)
    # in each reconstructed model:
    for fl in os.listdir(fld_md_pth):
        # fl = os.listdir(fld_md_pth)[0]
        pth = os.path.join(fld_md_pth, fl)
        if fl == 'KIDNEY_786O.mat':  # matlab doesn't accept variables starting with numbers, so the name is inverted
            fl = '786O_KIDNEY.mat'
        md = load_matlab_model(pth)
        md.solver.configuration.tolerances.feasibility = 1e-9
        cond = fl.split('.mat')[0]
        # change reactions like 'MAM01965e --> MAM01965b' to be 'MAM01965e --> ' (automatically sets them as exchanges):
        for rc in md.reactions:
            for mt, cf in rc.metabolites.items():
                if mt.compartment == 'b':
                    rc.subtract_metabolites({mt: cf})
        # remove metabolites with 'b' - boundary compartment (met. of reactions changed above):
        md.remove_metabolites([mt for mt in md.metabolites if mt.id.endswith('b')])
        # get reactions needed for cell line specific tasks:
        with open(consen_tsk_pth, mode='r') as f2:
            consen_tsk_lst = json.load(f2)
        present = {r.id for r in md.reactions}
        # get all reactions needed for non-essential tissue (consensus) tasks:
        consen_tsk_keep_pth = os.path.join(base_fld, 'consen_tasks_to_keep.tsv')
        consen_tsks_keep_df = pd.read_csv(consen_tsk_keep_pth, sep='\t', index_col=0)
        consen_ts = consen_tsks_keep_df[cond].astype(bool)
        consen_tsks_keep = consen_ts.index[consen_ts].astype(str)
        consen_tsks_tss = [t for t in consen_tsk_lst if t['name'] in consen_tsks_keep]
        # if with_tsk=True, all tissue specific tasks will be added (including demethylation tasks) + non-catalyzed demethylation reactions will be added
        # otherwise, non-catalyzed demethylation reactions will be added and two options exist:
        if with_tsk: # with cell-specific tasks
            grp_consen = set(reduce(lambda x, y: x + y, [t['mandatory_activity'] for t in consen_tsks_tss]))
        else: # without cell specific tasks
            if demeth_tsk: # but with demethylation tasks if task is done in the cell line
                consen_tsks_tss = [t for t in consen_tsks_tss if t['annotations']['system'] == 'DNA (DE)METHYLATION']
            else: # without all cell specific tasks (also the dmethylation tasks)
                consen_tsks_tss = list()
            if consen_tsks_tss:
                grp_consen = set(reduce(lambda x, y: x + y, [t['mandatory_activity'] for t in consen_tsks_tss]))
            else:
                grp_consen = set()
        # include un-catalyzed reactions of DNA demethylation (because of being un-catalyzed don't have a gene score),
        # that are non-essential for the task and therefore are almost always excluded during reconstruction:
        tokeep = grp_consen.union({'consdirectDNA5CaC', 'consdirectDNA5fC'}) - present
        tokeep_lst = [copy.deepcopy(generic_m.reactions.get_by_id(rid)) for rid in tokeep]
        md.add_reactions(tokeep_lst)
        # exclude models for which no DNA methylation reaction information is known:
        tss_spc_df = pd.read_excel(md_info_pth, sheet_name=md_info_sheet, skipfooter=11, index_col=0, skiprows=6)
        if np.isnan(tss_spc_df.loc['MAM01722n', cond]):  # if no information on DNA methylation is given
            print(f' No DNA methylation info available for {cond}')
        else:
            print(cond, md.reactions.get_by_id('prodDNAtot').reaction)
            # test objectives:
            md.objective = {md.reactions.get_by_id('adaptbiomass'): 1.0}
            with open(medium_path, mode='r') as f:
                hams_med_cmp = json.load(f)
            med_ex_ids = list(hams_med_cmp.values())  # exchange reactions that are part of media
            exch_const = {ex.id: (-1000.0, 1000.0) if ex.id in med_ex_ids else (0.0, 1000.0) for ex in md.exchanges}
            simul = get_simulator(md, envcond=exch_const)
            flx = {}  # to guarantee that if 'flx' variable is not obtained in one iteration (infeasible) the value of previous iteration is used
            try:
                sol = simul.simulate(method='pFBA', maximize=True)
                flx = sol.fluxes
                if sol.status.value == 'Optimal' and flx['adaptbiomass'] > 1e-9 and flx['MAR08641'] > 1e-9:
                    print(f'{cond} is growing: {flx["adaptbiomass"]}, and methylates: {flx["MAR08641"]} ')
                    if not os.path.exists(final_fld_md_pth):
                        os.makedirs(final_fld_md_pth)
                    if fl == '786O_KIDNEY.mat':
                        fl = 'KIDNEY_786O.mat'  # because gecko model will be build with matlab script and no matalab variable can start with numbers
                    save_matlab_model(md, os.path.join(final_fld_md_pth, fl))
                else:
                    print(f'{cond} not growing: {flx["adaptbiomass"]}')
            except:
                print(f'{cond} infeasible')

if __name__ == "__main__":
    GEN_MD_BASE_PTH = 'support/models/prodDNAtot.mat'
    FLD_MD_PTH = 'data/models_tINIT_human_pipe/init'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    FINAL_FLD_MD_PTH = 'support/models_tINIT_human_pipe/init'
    CONSEN_TSK_JSON_PATH = 'support/tasks/metabolicTasks_processed_consensus.json'
    BS_PTH = 'support/reconstruction/prodDNAtot_25_75_mean/task_protected/no_protect_inner_spont/no_wo_data_protected/fastcore'
    EXTRA_PATH = 'data/md_modify.xlsx'
    # without cell line specific tasks
    # + NO demethylation cell specific tasks:
    with_tsk = False
    final_fld_md_pth = FINAL_FLD_MD_PTH
    add_necessary_rc(fld_md_pth=FLD_MD_PTH, final_fld_md_pth=final_fld_md_pth,
                     gen_md_base_pth=GEN_MD_BASE_PTH,
                     consen_tsk_pth=CONSEN_TSK_JSON_PATH,
                     base_fld=BS_PTH,
                     md_info_pth=EXTRA_PATH,
                     md_info_sheet='DNAtot_coef',
                     medium_path=MEDIUM_PTH,
                     with_tsk=with_tsk)

    # with consensus (cell line specific) tasks:
    with_tsk = True
    if with_tsk: final_fld_md_pth = os.path.join(FINAL_FLD_MD_PTH, 'including_tsks')
    add_necessary_rc(fld_md_pth=FLD_MD_PTH, final_fld_md_pth=final_fld_md_pth,
                     gen_md_base_pth=GEN_MD_BASE_PTH,
                     consen_tsk_pth=CONSEN_TSK_JSON_PATH,
                     base_fld=BS_PTH,
                     md_info_pth=EXTRA_PATH,
                     md_info_sheet='DNAtot_coef',
                     medium_path=MEDIUM_PTH,
                     with_tsk=with_tsk)

    # without cell line specific tasks
    # + WITH demethylation cell specific tasks:
    with_tsk = False
    demeth_tsk = True
    final_fld_md_pth = os.path.join(FINAL_FLD_MD_PTH, 'notsk_wdemethtsk')
    add_necessary_rc(fld_md_pth=FLD_MD_PTH, final_fld_md_pth=final_fld_md_pth,
                     gen_md_base_pth=GEN_MD_BASE_PTH,
                     consen_tsk_pth=CONSEN_TSK_JSON_PATH,
                     base_fld=BS_PTH,
                     md_info_pth=EXTRA_PATH,
                     md_info_sheet='DNAtot_coef',
                     medium_path=MEDIUM_PTH,
                     with_tsk=with_tsk,
                     demeth_tsk=demeth_tsk)

