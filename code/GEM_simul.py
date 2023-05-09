home_path = '/home/tbarata'

import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])

import pandas as pd
from cobra.io import load_matlab_model
from cobra import Reaction
import json
import numpy as np
import re
from cell_lines import CellLineIds
from metbtasks import MetTask
from graphs import Graphs
from simulation import Simulation
from mewpy.simulation import get_simulator
from generic_model import GenericMdOperations

def simul_GEMS(gem_md_fld, flx_fld, constr_ex, envcond, smp_info, medium_pth, gr_const, with_tsk, algo, sim, cl_spec, md_info_pth, md_info_sheet, **kwargs):
    '''
    - simulations with GEMS
    :param gem_md_fld: - directory where folders with the reconstructed GEM models are.
                       - each folder name is a model/cell line name.
    :param flx_fld: folder where to save simulation results
    :param constr_ex: exchange reactions to which apply flux bound constraints
    :param envcond: path to dataframe with bounds to apply for some exchange reactions in different cell lines
    :param smp_info: path to file with cell line info
    :param obj_lst: list with objectives to analyze
    :param medium_pth: path to file with reactions of medium composition
    :param gr_const: whether to constraint growth with experimental growth rates
    :param with_tsk: bool, whether highest reaction score was given to necessary reaction of tissue specific tasks or not
    :param algo: either 'fastcore' or 'init'. algorithm with which models were reconstructed
    :param sim: simulation type either minimize sum of all fluxes ('min_sum_all_flx') or maximize biomass minimizing sum of all fluxes ('biomass_pFBA')
    :param cl_spec: bool, whether to use cell line specific reaction of 'prodDNAtot' (True) or use the generic instead (False).
    :param md_info_pth: path to excel file with tissue specific DNA methylation information
    :param md_info_sheet: name of excel sheet with tissue specific DNA methylation information
    :param kwargs: * rc_remove_lst: list of ids of reactions to remove
                   * rc_add_lst: list of dictionaries with keys 'rc_id', 'subsystem', 'lb', 'ub' and 'mtlst'
                                  each dictionary has the information of one reaction.
                                  'mtlst' is a string, e.g. 'subst_id:-coef, prod_id:coef'
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
    # load file with lower and upper bounds for experimental fluxes of some exchange reactions:
    if constr_ex: exp_df = exp_df.drop(constr_ex, axis=0)
    # load experimental cells' growth rates:
    grf = pd.read_excel(envcond, sheet_name='GrowthRate', index_col=0)
    grf['mean'] = grf.mean(axis=1)
    # load each model:
    if with_tsk: gem_md_fld = os.path.join(gem_md_fld, algo, 'including_tsks')
    else: gem_md_fld = os.path.join(gem_md_fld, algo)
    md_nms = [el for el in os.listdir(gem_md_fld) if '.mat' in el]
    all_flx_d = dict()
    exc_flx = dict()
    mass_final_flx = dict()
    feas=0
    unfeas=0
    for md_nm in md_nms:
        print(md_nm)
        # md_nm = md_nms[0]
        md_pth = os.path.join(gem_md_fld, md_nm)
        md_nm = md_nm.split('.mat')[0]
        print(md_nm)
        if md_nm == 'KIDNEY_786O':  # matlab doesn't accept variables starting with numbers, so the name is inverted
            md_nm = '786O_KIDNEY'
        if md_nm == 'MDAMB468_BREAST': # there's no info on exchange fluxes of 'MDAMB468_BREAST' (as cell mass info is not available),
            continue
        md = load_matlab_model(md_pth)
        md.solver.configuration.tolerances.feasibility = 1e-9
        # remove reactions if needed:
        if 'rc_remove_lst' in kwargs:
            md.remove_reactions(kwargs['rc_remove_lst'])
        # add reactions if needed:
        if 'rc_add_lst' in kwargs:
            for r in kwargs['rc_add_lst']:
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
                            'SF268_CENTRAL_NERVOUS_SYSTEM', 'SF539_CENTRAL_NERVOUS_SYSTEM',
                            'SNB75_CENTRAL_NERVOUS_SYSTEM', 'HOP92_LUNG']
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
        if constr_ex: env_const, constr_ex_name = Simulation.get_env_const_GEMs(m=md, m_nm=md_nm, constr_ex=constr_ex, med_ex_ids=med_ex_ids, envcond=envcond, sheet_name='ExpFlux_lb_ub', smp_info=smp_info)
        else: env_const, constr_ex_name = Simulation.get_env_const_GEMs(m=md, m_nm=md_nm, constr_ex=constr_ex, med_ex_ids=med_ex_ids)
        # add measured experimental values of growth rates if requested:
        internal_const = {}
        if gr_const: internal_const.update({'adaptbiomass': (grf.loc[md_nm, 'Growth MIN (1/h)'], grf.loc[md_nm, 'Growth MAX (1/h)'])})
        # run simulation:
        try:
            if sim == 'min_sum_all_flx':
                sol= MetTask.minSumFluxes(model=md, envcond=env_const, constraints=internal_const)
            elif sim == 'biomass_pFBA':
                simul = get_simulator(md, envcond=env_const, constraints=internal_const)
                simul.objective = {'adaptbiomass': 1.0}
                sol = simul.simulate(method='pFBA')
        except:
            print('weird unfeasible')
            unfeas += 1
            continue
        if sol.status.name == 'OPTIMAL':
            flx = {}  # to guarantee that if 'flx' variable is not obtained in one iteration (infeasible) the value of previous iteration is used
            if sim == 'min_sum_all_flx':
                flx = sol.values
            elif sim == 'biomass_pFBA':
                flx = sol.fluxes
            flx = {k: v for k, v in flx.items() if not re.search(r'\+$|\-$', k)} # the output produces split reactions of reversible reactions, besides the reversible reactions.
            flx = {k: 0.0 if np.abs(v) < 1e-9 else v for k, v in flx.items()}
            print(sol.status.name)
            feas += 1
            all_flx_d[md_nm] = flx
        else:
            print(sol.status.name)
            unfeas += 1
            continue
        # update dictionaries with simulated exchange reactions' flux values and growth rates:
        exc_flx, mass_final_flx = Simulation.get_exc_grw_flx(m=md, m_nm=md_nm, exp_df=exp_df, all_flx_d=all_flx_d, exc_flx=exc_flx, mass_final_flx=mass_final_flx, md_type='GEM')
    # save dataframe with all fluxes of all models:
    fld_fn = Simulation.save_all_flx_mds(all_flx_d=all_flx_d, obj_id=None, flx_fld=flx_fld, algo=algo, constr_ex_name=constr_ex_name, with_tsk=with_tsk, constr_ex=constr_ex, gr_const=gr_const, cl_spec=cl_spec)
    # scatter plot of exchange reactions fluxes + histograms:
    Graphs.plot_exc_flx(exp_df, exc_flx, fld_fn)
    # scatter plot of biomass reaction fluxes + boxplot with relative errors:
    Graphs.prep_plot_biomass(mass_final_flx, fld_fn, grf)
    print('feasible:', feas)
    print('unfeasible:', unfeas)

if __name__ == "__main__":
    ##############################################
    #  Loosely constrained                       #
    ##############################################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - FASTCORE without tasks
    # - medium open
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo = 'fastcore'
    with_tsk = False
    GEM_MD_FLD = 'support/models_richelle_pipe'
    FLX_FLD = 'results/GEM_simul_richelle/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex = False
    divide_by_biomass = False
    gr_const = False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet='DNAtot_coef'
    cl_spec=False # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
    ###########################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - FASTCORE WITH tasks
    # - medium open
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo = 'fastcore'
    with_tsk = True
    GEM_MD_FLD = 'support/models_richelle_pipe'
    FLX_FLD = 'results/GEM_simul_richelle/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex = False
    divide_by_biomass = False
    gr_const = False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
    ###########################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - INIT without tasks
    # - medium open
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo='init'
    with_tsk = False
    GEM_MD_FLD = 'support/models_tINIT_human_pipe'
    FLX_FLD = 'results/GEM_simul_human1/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv' # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex=False
    divide_by_biomass=False
    gr_const=False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
    ###########################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - INIT WITH tasks
    # - medium open
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo='init'
    with_tsk = True
    GEM_MD_FLD = 'support/models_tINIT_human_pipe'
    FLX_FLD = 'results/GEM_simul_human1/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv' # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex=False
    divide_by_biomass=False
    gr_const=False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
    ###########################

    ##############################################
    #  CONSTRAINED with some medium components   #
    ##############################################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - FASTCORE without tasks
    # - medium DEFINED
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo = 'fastcore'
    with_tsk = False
    GEM_MD_FLD = 'support/models_richelle_pipe'
    FLX_FLD = 'results/GEM_simul_richelle/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex=['MAR09034', 'MAR09135', 'MAR09044'] # glucose, lactate, threonine
    divide_by_biomass = False
    gr_const = False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
    ###########################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - FASTCORE WITH tasks
    # - medium DEFINED
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo = 'fastcore'
    with_tsk = True
    GEM_MD_FLD = 'support/models_richelle_pipe'
    FLX_FLD = 'results/GEM_simul_richelle/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex=['MAR09034', 'MAR09135', 'MAR09044'] # glucose, lactate, threonine
    divide_by_biomass = False
    gr_const = False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
    ###########################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - INIT without tasks
    # - medium Defined
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo='init'
    with_tsk = False
    GEM_MD_FLD = 'support/models_tINIT_human_pipe'
    FLX_FLD = 'results/GEM_simul_human1/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv' # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex=['MAR09034', 'MAR09135', 'MAR09044'] # glucose, lactate, threonine
    divide_by_biomass=False
    gr_const=False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
    ###########################

    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - INIT WITH tasks
    # - medium Defined
    # - no constraints on biomass
    # - pFBA with biomass as objective
    ###########################
    sim = 'biomass_pFBA'
    algo='init'
    with_tsk = True
    GEM_MD_FLD = 'support/models_tINIT_human_pipe'
    FLX_FLD = 'results/GEM_simul_human1/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv' # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex=['MAR09034', 'MAR09135', 'MAR09044'] # glucose, lactate, threonine
    divide_by_biomass=False
    gr_const=False
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)

    ###############################
    #  CONSTRAINED with biomass   #
    ###############################
    ###########################
    # - model without tissue specific DNA methylation proportions
    # - GEM
    # - INIT without tasks
    # - medium open
    # - WITH constraints on biomass
    # - FBA with MIN SUM OF ALL FLUXES
    ###########################
    sim = 'min_sum_all_flx'
    algo = 'init'
    with_tsk = False
    GEM_MD_FLD = 'support/models_tINIT_human_pipe'
    FLX_FLD = 'results/GEM_simul_human1/'
    MEDIUM_PTH = 'data/hams_medium_composition.json'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    constr_ex = False
    divide_by_biomass = False
    gr_const = True
    EXTRA_PATH = 'data/md_modify.xlsx'
    md_info_sheet = 'DNAtot_coef'
    cl_spec = False  # without cell line specific total DNA reaction
    simul_GEMS(gem_md_fld=GEM_MD_FLD, flx_fld=FLX_FLD, constr_ex=constr_ex, envcond=ENVCOND, smp_info=ACHILLES_SMP_INFO,
               medium_pth=MEDIUM_PTH, gr_const=gr_const, with_tsk=with_tsk, algo=algo, sim=sim,
               cl_spec=cl_spec, md_info_pth=EXTRA_PATH, md_info_sheet=md_info_sheet)
