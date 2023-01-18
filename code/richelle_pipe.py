home_path = '/home/tbarata'

import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])
# following only need to include if using init from troppo - use dev versions of troppo and cobamp when using init:
# sys.path.extend([os.path.join(home_path, pkg) for pkg in ['troppo/src', 'cobamp/src']])
import numpy as np
from reconstruction import Reconstruction
from generic_model import GenericMdOperations
from metbtasks import MetTask

def get_rc_models(generic_path, m_name, obj_id, cl_path, transcriptomics, ug_p, lg_p, lc_p, wo_data_protected, gn_sc_path,
                       int_str, protect_inner_spont, rc_sc_path, envcond, medium_path, original_task_path, essen_tsk_json_path,
                       tsk_json_path, path_rc_wgts, reconst_fld, reconst_file_path, flx_file_path, algo, constr_ex,
                       simul, md_info_pth, md_info_sheet, tss_spc_md_fld_path, **kwargs):
    '''
    - obtains reconstructed models for certain combinations of choices - protect/ not protect inner spontaneuous reactions, reactions without gene expression, tasks, etc
    :param generic_path: path to generic xml model
    :param m_name: model name. usefull when dealing with different methylation models
    :param obj_id: objective. it's a dictionary of type reaction_id: reaction weight
    :param cl_path: path json to file with cell lines to use
    :param transcriptomics: path to file with transcriptomics data to preprocess
    :param ug_p: percentile used to obtain the upper global threshold
    :param lg_p: percentile used to obtain the lower global threshold
    :param lc_p: percentile used to obtain the local threshold. if False, the mean is used instead
    :param wo_data_protected: whether to protect genes for which there aren't gene score information. either True or False
    :param gn_sc_path: path to file where to save gene scores
    :param int_str: integration strategy to convert gene scores to reaction scores using the gene-reaction (AND/OR) rules
                        (function_for_AND, function_for_OR)
    :param protect_inner_spont: boolean that indicates whether to protect (give highest score) to reactions that do not have gene reaction rule
    :param rc_sc_path: path to where to save the reaction scores
    :param envcond: path to file with exchange flux constraints and biomass constraints
    :param medium_path: path to .json file with medium composition. contains a dict of the format {metabolite_id: reaction_id}
    :param original_task_path: path to file with original tasks
    :param tsk_json_path: path to list with tasks. each task in troppo format and with mandatory activity
    :param essen_tsk_json_path: path to list with essential tasks. each task in troppo format and with mandatory activity
    :param path_rc_wgts: path to where to save the reaction weights
    :param reconst_fld: folder where to save results
    :param reconst_file_path: path to file where to save the reactions to keep/ or not for the algorithm
    :param flx_file_path: path to file where to save the fluxes after reconstruction and pFBA for provided objective
    :param algo: algorithm to use in the reconstruction
    :param constr_ex: exchange reactions to which apply flux bound constraints. False if there are none.
    :param simul: either pFBA or FBA is used in simulation
    :param md_info_pth: path to excel file with tissue specific DNA methylation information
    :param md_info_sheet: name of excel sheet with tissue specific DNA methylation information
    :param tss_spc_md_fld_path: path to folder where to save cell specific reconstructed models
    :param kwargs - path_tsk_keep: path to file with tasks that should be protected during model reconstruction for each of the cell types/conditions based on the task metabolic score
                  - el: user provided threshold for Achilles experimental lethality score
                  - sl: user provided threshold for simulated lethality score
    :return: base folder where all results are saved
    '''
    lc_p_id = 'mean' if not lc_p else lc_p
    base_fld = os.path.join(reconst_fld, f'{m_name}_{lg_p}_{ug_p}_{lc_p_id}')
    base_fld = os.path.join(base_fld, 'task_protected') if 'path_tsk_keep' in kwargs else os.path.join(base_fld, 'no_task_protected')
    # base_fld = os.path.join(base_fld, 'task_protected')
    for k, v in {'protect_inner_spont': protect_inner_spont, 'wo_data_protected': wo_data_protected}.items():
        base_fld = os.path.join(base_fld, k) if v else os.path.join(base_fld, 'no_' + k)
    # for k, v in {'protect_inner_spont': protect_inner_spont, 'wo_data_protected': wo_data_protected}.items():
    #   base_fld = os.path.join(base_fld, 'no_' + k)
    base_fld = os.path.join(base_fld, algo)
    if not os.path.exists(base_fld):
        os.makedirs(base_fld)
    paths = {k: os.path.join(base_fld, v) for k, v in {'rc_sc_path': rc_sc_path, 'gn_sc_path': gn_sc_path, 'path_rc_wgts': path_rc_wgts}.items()}
    if 'path_tsk_keep' in kwargs:
        path_tsk_keep = os.path.join(base_fld, kwargs['path_tsk_keep'])
    # path_tsk_keep = os.path.join(base_fld, path_tsk_keep)
    nthreads = 32 # optional
    m = GenericMdOperations(xml_path=generic_path, fea_tol=1e-9).model

    ## get the reaction scores with the best threshold combination:
    # obtain cell lines to use and preprocess transcriptomics data.
    # - deepmap/CCLE transcriptomics sample description doesn't have NCIADRRES cell line (so, only 59 cell lines are used at this point):
    rec = Reconstruction(model=m, cl_path=cl_path)
    # - only 47 of the 59 cell lines were on deepmap/CCLE transcriptomics data. So, only those were tested from this point onwards:
    rec.preprocess_transcriptomics(data=transcriptomics)
    rec.get_gene_scores(ug_p=ug_p, lg_p=lg_p, lc_p=lc_p, protect=wo_data_protected, gene_path=paths['gn_sc_path'])
     # protect=True protects (gives highest score) genes for which there aren't gene score information
    rec.get_reaction_scores(int_str=int_str, protect=protect_inner_spont, path=paths['rc_sc_path'], medium_path=medium_path, nthreads=nthreads)
     # protect=True protects (gives highest score) to reactions that do not have gene reaction rule - uncatalyzed reactions

    # check which tasks should be done by each of the 43 models depending on the task metabolic score:
    thr = 5 * np.log(2)
    mtask = MetTask(task_file_path=original_task_path, model=m)
    if 'path_tsk_keep' in kwargs:
        mtask.tasks_to_keep_reconst(path_rc_sc=paths['rc_sc_path'], thr=thr, path_out=path_tsk_keep, tsk_path=tsk_json_path)

    # get reaction weights:
    protected = ['adaptbiomass']
    tsk_lst = mtask.load_tasks(path=tsk_json_path)
    essen_tsk_lst = mtask.load_tasks(path=essen_tsk_json_path)
    # when path_tsk_keep exists, there is task protection, otherwise, no:
    if 'path_tsk_keep' in kwargs:
        rec.weight_fc(rc_sc_path=paths['rc_sc_path'], alg=algo, protect_rc=protected, path_rc_wgts=paths['path_rc_wgts'], essen_tsk_lst=essen_tsk_lst, tsk_protect_path=path_tsk_keep, tsk_lst=tsk_lst)
    else:
        rec.weight_fc(rc_sc_path=paths['rc_sc_path'], alg=algo, protect_rc=protected, path_rc_wgts=paths['path_rc_wgts'], essen_tsk_lst=essen_tsk_lst)

    # get reactions to keep in the reconstructed model:
    reconst_path = reconst_file_path.format(algo=algo)
    rc_rc_path = os.path.join(base_fld, reconst_path)
    rec.apply_alg_all(model=m, alg=algo, path=rc_rc_path, thr=thr, nthreads=nthreads)
    # rec.apply_alg_all(model=m, alg='init', path=rc_rc_path, thr=thr, nthreads=nthreads)

    # build models with tailored objective and do pFBA:
    if simul == 'pFBA':
        fba = False
    elif simul == 'FBA':
        fba = True
    flx_path = os.path.join(base_fld, flx_file_path.format(algo=algo, obj=obj_id, const=constr_ex, simul=simul))
    if 'path_tsk_keep' in kwargs:
        Reconstruction.reconstruct_mds_batch(algo_reactions_pth=rc_rc_path, medium_path=medium_path, generic_path=generic_path, obj_id=obj_id, flx_path=flx_path,
                                             fba=fba, include_algo_res=True, md_info_pth=md_info_pth,
                                             md_info_sheet=md_info_sheet, tss_spc_md_fld_path=tss_spc_md_fld_path, envcond=envcond, with_tsk=True)
    else:
        Reconstruction.reconstruct_mds_batch(algo_reactions_pth=rc_rc_path, medium_path=medium_path, generic_path=generic_path, obj_id=obj_id, flx_path=flx_path,
                                         fba=fba, include_algo_res=True, md_info_pth=md_info_pth, md_info_sheet=md_info_sheet, tss_spc_md_fld_path=tss_spc_md_fld_path, envcond=envcond, with_tsk=False)


if __name__ == "__main__":
    METH_MD_FLD = 'support/models'
    GN_MD_PTH = os.path.join(METH_MD_FLD, 'prodDNAtot.xml')
    TRANSCRIPTOMICS = 'data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.txt'
    NCI60_PTH = 'data/cell_lines.json'
    MD_PATH_SBML = 'data/models/Human-GEM.xml'
    EXTRA_PATH = 'data/md_modify.xlsx'
    DNA_METH_RC_PATH = 'data/DNAmeth_reactions.json'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    GN_SC_FILE_PATH = 'gene_scores.tsv'  # file where to save gene scores
    RC_SC_FILE_PATH = 'reaction_scores.tsv'
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    MEDIUM_PATH = 'data/hams_medium_composition.json'
    CONSEN_TASK_PATH = 'data/tasks/metabolicTasks_CellfieConsensus.txt'
    CONSEN_TSK_JSON_PATH = 'support/tasks/metabolicTasks_processed_consensus.json'
    SUBSYS_PATH = 'data/tasks/pcbi.1006867.s005.xlsx'  # file just for system info (do not use for tasks - as those are from recon2.2.)
    ESENT_TASK_PATH = 'data/tasks/metabolicTasks_Essential.txt'
    ESSEN_TSK_JSON_PATH = 'support/tasks/metabolicTasks_processed_essential.json'
    WEIGHTS_FILE_PATH = 'reaction_weights.tsv'
    FLX_FILE_PATH = 'fluxes_{algo}_{obj}_constraints_{const}_{simul}.tsv'
    RECONST_FLD = 'support/reconstruction'
    RECONST_FILE_PATH = '{algo}_reactions_to_keep.tsv'
    MTH_CORR_FLD_PATH = 'support/methyl_corr'
    TSK_KEEP_FILE_PATH = 'consen_tasks_to_keep.tsv'
    TSS_SPC_MD_FLD_PATH = 'support/models_richelle_pipe'

    # - read human1.12 master branch in yml format and convert to SBML.
    #   if already loaded, it opens the SBML model
    #   note: we didn't use xml directly cause was giving segmentation error when creating the methylation models bellow.
    # - create generic methylation model, remove blocked reactions and check whether it grows and (de)methylation reactions aren't blocked
    # - adapt rules with complexes of isozymes, e.g.:
    #      *  "((G1 and G2) or (G2 and G3)) and ((G5 and G6) or (G7 and G8))" --> "(G1 and G2 and G5 and G6) or (G1 and G2 and G7 and G8) or (G2 and G3 and G5 and G6) or (G2 and G3 and G7 and G8)"
    #      *  "(G1 or G2) and (G3 or G4)" --> "(G1 and G3) or (G1 and G4) or (G2 and G3) or (G2 and G4)"
    #   it's a requirement to use the script that creates gecko models.
    gen = GenericMdOperations(xml_path=MD_PATH_SBML, fea_tol=1e-9)
    gen.create_generic_meth_md(extra_path=EXTRA_PATH, dna_meth_rc_path=DNA_METH_RC_PATH, meth_md_fld=METH_MD_FLD)

    # load (de)methylation model:
    m = GenericMdOperations(xml_path=GN_MD_PTH, fea_tol=1e-9).model

    ## prepare tasks to be used:
    # the consensus tasks - from richelle et al. article
    # task 37 had an NA (cause no corresponding metabolite id in human1 was found for that metabolite)
    # we changed NA in that task to glycogenin (see here https://github.com/SysBioChalmers/Human-GEM/discussions/376)
    # Task name/number in json format is not the same as in data .txt files.
    # So, one should always look to task description to make the correspondence to data .txt files!
    consen_mtask = MetTask(task_file_path=CONSEN_TASK_PATH, model=m)
    consen_mtask.prepare_tasks(path_save=CONSEN_TSK_JSON_PATH, sys_path=SUBSYS_PATH)
    # essential tasks from human1 article:
    essen_mtask = MetTask(task_file_path=ESENT_TASK_PATH, model=m)
    essen_mtask.prepare_tasks(path_save=ESSEN_TSK_JSON_PATH, sys_path='inside')  # task (sub)system information is inside the file

    # Notes: - it's been discovered that NCIADRRES is an ovarian and Not breast cancer cell line as originally thought.
    #        - some people call MDAN to MDAMB435.
    #          However, MDAN (not included in NCI60 cell line list) is a cell line derived from MDAMB435, but it's not an original NCI60 cell line.
    #        - both MDAMB435 and MDAN are melanoma cell lines, although initially were thought to be breast cell lines.

    algo='fastcore'
    obj_id = {'adaptbiomass': 1}
    m_name = 'prodDNAtot'
    ##### protect cell type specific tasks ####
    get_rc_models(generic_path=GN_MD_PTH, m_name=m_name, obj_id=obj_id, cl_path=NCI60_PTH, transcriptomics=TRANSCRIPTOMICS,
        ug_p=75, lg_p=25, lc_p=False,
        wo_data_protected=False,
        gn_sc_path=GN_SC_FILE_PATH,
        int_str=(min, max),
        protect_inner_spont=False,
        rc_sc_path=RC_SC_FILE_PATH,
        envcond=ENVCOND,
        medium_path=MEDIUM_PATH,
        original_task_path=CONSEN_TASK_PATH,
        tsk_json_path=CONSEN_TSK_JSON_PATH,
        essen_tsk_json_path=ESSEN_TSK_JSON_PATH,
        path_rc_wgts=WEIGHTS_FILE_PATH,
        reconst_fld=RECONST_FLD,
        reconst_file_path=RECONST_FILE_PATH,
        flx_file_path=FLX_FILE_PATH,
        algo=algo,
        constr_ex=False,
        path_tsk_keep=TSK_KEEP_FILE_PATH,
        simul='pFBA',
        md_info_pth=EXTRA_PATH,
        md_info_sheet='DNAtot_coef',
        tss_spc_md_fld_path=TSS_SPC_MD_FLD_PATH)

    ##### don't protect cell type specific tasks ####
    get_rc_models(generic_path=GN_MD_PTH, m_name=m_name, obj_id=obj_id, cl_path=NCI60_PTH,
        transcriptomics=TRANSCRIPTOMICS,
        ug_p=75, lg_p=25, lc_p=False,
        wo_data_protected=False,
        gn_sc_path=GN_SC_FILE_PATH,
        int_str=(min, max),
        protect_inner_spont=False,
        rc_sc_path=RC_SC_FILE_PATH,
        envcond=ENVCOND,
        medium_path=MEDIUM_PATH,
        original_task_path=CONSEN_TASK_PATH,
        tsk_json_path=CONSEN_TSK_JSON_PATH,
        essen_tsk_json_path=ESSEN_TSK_JSON_PATH,
        path_rc_wgts=WEIGHTS_FILE_PATH,
        reconst_fld=RECONST_FLD,
        reconst_file_path=RECONST_FILE_PATH,
        flx_file_path=FLX_FILE_PATH,
        algo=algo,
        constr_ex=False,
        simul='pFBA',
        md_info_pth=EXTRA_PATH,
        md_info_sheet='DNAtot_coef',
        tss_spc_md_fld_path=TSS_SPC_MD_FLD_PATH)




