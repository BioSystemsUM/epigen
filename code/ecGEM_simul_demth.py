home_path = '/home/tbarata'

import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])

import pandas as pd
from cobra.io import load_matlab_model
import json
from test_mds import TestMds
from cell_lines import CellLineIds
from graphs import Graphs
from simulation import  Simulation
from cobra import Reaction
from generic_model import GenericMdOperations
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats.mstats import spearmanr

def simul_ecGEMS(ecgem_md_fld, prot_md, ec_flx_fld, constr_ex, envcond, smp_info, medium_pth,
                 obj_id, gr_const, with_tsk, maxm, algo, flx_pror_rules, prot_limit, cl_spec,
                 md_info_pth, md_info_sheet, **kwargs):
    '''
    - simulations with ecGEMS
    :param ecgem_md_fld: - directory where folders with the reconstructed ecModels are.
                         - each folder name is a model/cell line name.
    :param prot_md: name of file containing the ecModel (either the protein pool or without pool or a mix)
    :param ec_flx_fld: folder where to save simulation results
    :param constr_ex: exchange reactions to which apply flux bound constraints OR False when no env constraints besides opening medium are applied
    :param envcond: path to dataframe with bounds to apply for some exchange reactions in different cell lines
    :param smp_info: path to file with cell line info
    :param medium_pth: path to file with reactions of medium composition
    :param obj_id: metabolic objective to apply
    :param gr_const: whether to constraint cell growth with experimental growth rates
    :param with_tsk: bool, whether necessary reactions of tissue specific tasks were included or not
    :param generic_md_name: name of folder where to save results
    :param maxm: whether to maximize or minimize the objective
    :param algo: either 'fastcore' or 'init'. algorithm with which models were reconstructed
    :param flx_pror_rules: whether to set flux proportion rules or not
    :param prot_limit: whether to set ub to reaction prot_pool_exchange or not.
                       if 'human_default', a default value of total protein for human is used. the value corresponds to the protein mass fraction used in the Human1 biomass reaction
    :param cl_spec: bool, whether to use cell line specific reaction of 'prodDNAtot' (True) or use the generic instead (False).
    :param md_info_pth: path to excel file with tissue specific DNA methylation information
    :param md_info_sheet: name of excel sheet with tissue specific DNA methylation information
    :param kwargs: * rc_remove_lst: list of ids of reactions to remove
                   * rc_add_lst: list of dictionaries with keys 'rc_id', 'subsystem', 'lb', 'ub' and 'mtlst'
                                  each dictionary has the information of one reaction.
                                  'mtlst' is a string, e.g. 'subst_id:-coef, prod_id:coef'
                   * extra_const: dict {rcid: (lb, ub)} for internal reactions where we want to add constrains
                   * prot_wf_raw: pandas series with total protein of each cell line in gprot/gDW
                   * sigma: parameter for saturation of enzymes in vivo. needed if prot_wf_raw is given
                   * frc: fraction of protein with enzyme function in the model. needed if prot_wf_raw is given
                   * enz_pfba: True/False, True when wanting to do an enzymatic pFBA (maximize objective, set the max objective value and minimize total protein)
    '''
    # load list with ids of medium exchange reactions:
    with open(medium_pth, mode='r') as f:
        hams_med_cmp = json.load(f)
    med_ex_ids = list(hams_med_cmp.values())
    # load file with experimental fluxes of some exchange reactions:
    exp_df = pd.read_excel(envcond, sheet_name='ExpFlux', index_col=2)
    exp_df = exp_df.dropna()
    CellLineIds.convert_cell_line_names(file=exp_df, smp_info=smp_info)
    exp_df = exp_df.iloc[:, 3:]
    if constr_ex: exp_df = exp_df.drop(constr_ex, axis=0)
    # load experimental cells' growth rates:
    grf = pd.read_excel(envcond, sheet_name='GrowthRate', index_col=0)
    grf['mean'] = grf.mean(axis=1)
    # load each model:
    if with_tsk: ecgem_md_fld = os.path.join(ecgem_md_fld, algo, 'including_tsks')
    else: ecgem_md_fld = os.path.join(ecgem_md_fld, algo, 'no_tsks')
    md_nms = os.listdir(ecgem_md_fld)
    all_flx_d = dict()
    exc_flx = dict()
    mass_final_flx = dict()
    feas=0
    unfeas=0
    for md_nm in md_nms:
        # md_nm = md_nms[0]
        print(md_nm)
        md_pth = os.path.join(ecgem_md_fld, md_nm, prot_md)
        if md_nm == 'KIDNEY_786O':  # matlab doesn't accept variables starting with numbers, so the name is inverted
            md_nm = '786O_KIDNEY'
        if md_nm == 'MDAMB468_BREAST': # there's no info on exchange fluxes of 'MDAMB468_BREAST' (as cell mass info is not available),
            continue
        # if required, apply constraint to prot_pool_exchange ' --> prot_pool':
        md = load_matlab_model(md_pth)
        md.solver.configuration.tolerances.feasibility = 1e-9
        if prot_limit=='human_default':
            md.reactions.get_by_id('prot_pool_exchange').upper_bound = 0.593 * kwargs['sigma'] * kwargs['frc']
            md.reactions.get_by_id('prot_pool_exchange').lower_bound = 0.0
        elif prot_limit:
            PTotal = kwargs['prot_wf_raw'].loc[md_nm]
            md.reactions.get_by_id('prot_pool_exchange').upper_bound = PTotal * kwargs['sigma'] * kwargs['frc']
            md.reactions.get_by_id('prot_pool_exchange').lower_bound = 0.0
        else:
            md.reactions.get_by_id('prot_pool_exchange').bounds = (0.0, float('inf'))
        # remove reactions needed:
        if 'rc_remove_lst' in kwargs:
            md.remove_reactions(kwargs['rc_remove_lst'])
        # add reactions needed:
        if 'rc_add_lst' in kwargs:
            for r in kwargs['rc_add_lst']:
                # for r in rc_add_lst:
                if r['rc_id'] not in [reac.id for reac in md.reactions]:
                    rc = Reaction(id=r['rc_id'], subsystem=r['subsystem'], lower_bound=r['lb'], upper_bound=r['ub'])
                    md.add_reactions([rc])
                    GenericMdOperations.set_metab(mtlst=r['mtlst'], model=md, rid=r['rc_id'])
                    print(md.reactions.get_by_id(r['rc_id']).reaction)
        if cl_spec:
            # cell lines with both information for DNA5mC and DNA5hmC:
            cl_both_info = ['CAKI1_KIDNEY', 'NCIH226_LUNG', 'NCIH322M_LUNG', 'SKMEL28_SKIN',
            'OVCAR4_OVARY', 'HCT15_LARGE_INTESTINE', 'SF295_CENTRAL_NERVOUS_SYSTEM', 'UACC62_SKIN',
            'NCIH460_LUNG', 'SKOV3_OVARY', 'A498_KIDNEY', 'NIHOVCAR3_OVARY',
            'OVCAR8_OVARY', 'NCIH23_LUNG', '786O_KIDNEY', 'IGROV1_OVARY',
            'MALME3M_SKIN', 'SW620_LARGE_INTESTINE', 'HCT116_LARGE_INTESTINE',
            'SKMEL2_SKIN', 'A549_LUNG', 'NCIH522_LUNG', 'U251MG_CENTRAL_NERVOUS_SYSTEM',
            'ACHN_KIDNEY', 'SKMEL5_SKIN', 'LOXIMVI_SKIN', 'UACC257_SKIN',
            'HOP62_LUNG', 'EKVX_LUNG', 'OVCAR5_OVARY', 'UO31_KIDNEY',
            'SF268_CENTRAL_NERVOUS_SYSTEM', 'SF539_CENTRAL_NERVOUS_SYSTEM', 'SNB75_CENTRAL_NERVOUS_SYSTEM', 'HOP92_LUNG']
            if md_nm not in cl_both_info:
                continue
            # replace DNAtotal formation reaction by one that takes into account tissue specific methylated vs unmethylated DNA ratios:
            md.remove_reactions(['prodDNAtot'])
            tss_spc_df = pd.read_excel(md_info_pth, sheet_name=md_info_sheet, skipfooter=11, index_col=0, skiprows=6)
            tss_spc_df.rename(index={'DNA5hmCn ***': 'DNA5hmCn'}, inplace=True)
            rc = Reaction(id='prodDNAtot', subsystem='', lower_bound=0.0, upper_bound=1000.0)
            md.add_reactions([rc])
            mtblst = ", ".join([mtb + ":" + str(-tss_spc_df.loc[mtb, md_nm]) for mtb in ["MAM01722n", "DNA5hmCn", "MAM01721n", "DNA5fCn"]]) + ", DNAtotn:1.0"
            GenericMdOperations.set_metab(mtlst=mtblst, model=md, rid='prodDNAtot')
        print(md_nm, md.reactions.get_by_id('prodDNAtot').reaction)
        # get environmental constraints:
        if constr_ex: exch_const, constr_ex_name = Simulation.get_env_const_ecGEMs(m=md, m_nm=md_nm, constr_ex=constr_ex, med_ex_ids=med_ex_ids, envcond=envcond, sheet_name='ExpFlux_lb_ub', smp_info=smp_info)
        else: exch_const, constr_ex_name = Simulation.get_env_const_ecGEMs(m=md, m_nm=md_nm, constr_ex=constr_ex, med_ex_ids=med_ex_ids)
        # add measured experimental values of growth rates if requested:
        internal_const = {}
        if gr_const: internal_const.update({'adaptbiomass': (grf.loc[md_nm, 'Growth MIN (1/h)'], grf.loc[md_nm, 'Growth MAX (1/h)'])})
        # do FBA while adding flux proportion constraints for DNA methylation reactions if required:
        try:
            if 'enz_pfba' in kwargs:
                res = Simulation.fba_ecGEMs(md=md, exch_const=exch_const, internal_const=internal_const, obj_id=obj_id, maxm=maxm, feas=feas, unfeas=unfeas, all_flx_d=all_flx_d, md_nm=md_nm, flx_pror_rules=flx_pror_rules, enz_pfba=True)
            else:
                res = Simulation.fba_ecGEMs(md=md, exch_const=exch_const, internal_const=internal_const, obj_id=obj_id, maxm=maxm, feas=feas, unfeas=unfeas, all_flx_d=all_flx_d, md_nm=md_nm, flx_pror_rules=flx_pror_rules)
            if isinstance(res, int):
                unfeas = res
                continue
            else: feas, all_flx_d = res
        except:
            print('weird unfeasible')
            unfeas += 1
            continue
        # update dictionaries with simulated exchange reactions' flux values and growth rates:
        exc_flx, mass_final_flx = Simulation.get_exc_grw_flx(m=md, m_nm=md_nm, exp_df=exp_df, all_flx_d=all_flx_d, exc_flx=exc_flx, mass_final_flx=mass_final_flx, md_type='ecGEM')
    # save dataframe with all fluxes of all models:
    fld_fn = Simulation.save_all_flx_mds(all_flx_d=all_flx_d, obj_id=obj_id, flx_fld=ec_flx_fld, algo=algo, constr_ex_name=constr_ex_name, with_tsk=with_tsk, constr_ex=constr_ex, gr_const=gr_const, cl_spec=cl_spec)
    # scatter plot of exchange reactions fluxes + histograms:
    Graphs.plot_exc_flx(exp_df, exc_flx, fld_fn)
    # scatter plot of biomass reaction fluxes + boxplot with relative errors:
    Graphs.prep_plot_biomass(mass_final_flx, fld_fn, grf)
    print('feasible:', feas)
    print('unfeasible:', unfeas)



#########################################################################################
    # - Cell specific DNA methylation                                                        #
    # - with biomass constraint
    # - Exometabolite UNconstrained with
    # minimize total protein AFTER maximizing different weights of demethylation
    ##########################################################################################
    ###########################
    # - WITH tissue specific DNA methylation proportions
    # - ecGEM
    # - INIT without tasks - the best strategy to ecGEMs
    # - medium open
    # - CONstraints on biomass
    # - minimize total protein AFTER maximize dif. weights of demethylation as objective
    # - No flux proportion rules
    ###########################
    obj_id={'prot_pool_exchange': -1.0}
    obj_id_1st = {'arm_ligateDNA':1.0, 'consdirectDNA5fC':1.0, 'consdirectDNA5CaC':1.0}
    obj_id_2nd = {'prot_pool_exchange': -1.0}
    maxm = True
    algo = 'init'
    with_tsk = False
    prot_limit = False  # no limitation on total enzyme usage
    constr_ex = False
    gr_const = True
    flx_pror_rules = True # No flux rules
    ecGEM_MD_FLD = 'support/ecGEMs_human1'
    EC_FLX_FLD = 'results/ecGEM_simul_human1'
    PROT_POOL_MD = 'ecModel_batch.mat'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    divide_by_biomass = False
    exp_biomass = False # don't compare experimental biomass flux with experimental methylation
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = True  # with cell line specific total DNA reaction
    simul_ecGEMS(ecgem_md_fld=ecGEM_MD_FLD, prot_md=PROT_POOL_MD, ec_flx_fld=EC_FLX_FLD, constr_ex=constr_ex,
                 envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO, medium_pth=MEDIUM_PTH, obj_id=obj_id,
                 gr_const=gr_const, with_tsk=with_tsk, maxm=maxm, algo=algo, flx_pror_rules=flx_pror_rules,
                 prot_limit=prot_limit, cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)

    flx['arm_ligateDNA']


    list_meth_rc = ['arm_MAR08641', 'arm_MAR03875']
    compare_methyl_flx(mth_corr_fld_path=EC_FLX_FLD,
                       methlt_dt_fld=METHLT_FLD, constr_ex=constr_ex,
                       smp_info=ACHILLES_SMP_INFO, list_meth_rc=list_meth_rc, divide_by_biomass=divide_by_biomass,
                       with_tsk=with_tsk, obj_id=obj_id, algo=algo, exp_biomass=exp_biomass, envcond=ENVCOND)
    list_meth_rc = ['consdirectDNA5fC', 'consdirectDNA5CaC', 'arm_prodDNA5CaC', 'arm_prodDNA5mU', 'prodAPsite3No1', 'prodAPsite4No1']
    compare_methyl_flx(mth_corr_fld_path=EC_FLX_FLD,
                       methlt_dt_fld=METHLT_FLD, constr_ex=constr_ex,
                       smp_info=ACHILLES_SMP_INFO, list_meth_rc=list_meth_rc, divide_by_biomass=divide_by_biomass,
                       with_tsk=with_tsk, obj_id=obj_id, algo=algo, exp_biomass=exp_biomass, envcond=ENVCOND)
