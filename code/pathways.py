home_path = '/home/tbarata'
import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])

from pathway_analysis import PathAnalysis

if __name__ == "__main__":
    ## fluxes - cell lines only - use boxplots from here
    # heatmap mean, log10, No norm min-max on rows
    # boxplot mean
    # cell lines
    EC_FLX_FLD = 'results/ecGEM_simul_human1'
    METHLT_FL_RC = 'data/DNAmeth_reactions.json'
    GEN_MD_PTH = 'data/models/Human-GEM.yml'
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    with_tsk = False  # tasks were not applied during reconstruction
    algo = 'init'  # algorithm used in reconstruction was tINIT
    constr_ex = False  # no constraints were applied to exchange fluxes of exometabolites
    obj_id = {'prot_pool_exchange': 1.0}  # objective used in the flux simulation (dict) - it was minimized
    op = 'mean'  # aggregate flux values of each subsystem by mean
    htmp_log = True  # logaritmize values for heatmap
    trsf = 'none'  # do not apply any standartization or normalization to values in heatmap
    per_tissue = False  # run for cell lines
    gr_const = True # it was used constraint in biomass in final simulations
    cl_spec = True  # it was used cell line specific total DNA reaction in final simulations
    PathAnalysis.analyze_pathway_flux(mth_corr_fld_path=EC_FLX_FLD, algo=algo, obj_id=obj_id, constr_ex=constr_ex,
                         gen_md_pth=GEN_MD_PTH, gen_gecko_pth=GEN_GECKO_PTH, methlt_fl_rc=METHLT_FL_RC, with_tsk=with_tsk,
                         htmp_log=htmp_log, per_tissue=per_tissue, op=op, trsf=trsf,
                         cell_width=0.7, cell_height=0.5, fn=12, xlab_siz=17, ylab_siz=16, xlab_rot=90,
                         xcbar=2.2, ycbar=0.65, cbar_width=0.05, cbar_height=0.09, xheat=1, yheat=1, heat_wf=1,
                         heat_hf=1, rden_xf=1.36, cden_xf=1, rden_yf=1, cden_yf=0.92, rden_wf=1.1, cden_wf=1, rden_hf=1,
                         cden_hf=0.55, cbar_lab_siz=16, rclust=False, cclust=True,
                         gr_const=gr_const, cl_spec=cl_spec)

    ## protein usage cell lines only - use boxplots from here
    # heatmap mean, log10, No norm min-max on rows
    # boxplot mean
    # cell lines
    EC_FLX_FLD = 'results/ecGEM_simul_human1'
    METHLT_FL_RC = 'data/DNAmeth_reactions.json'
    GEN_MD_PTH = 'data/models/Human-GEM.yml'
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    with_tsk = False  # tasks were not applied during reconstruction
    algo = 'init'  # algorithm used in reconstruction was tINIT
    constr_ex = False  # no constraints were applied to exchange fluxes of exometabolites
    obj_id = {'prot_pool_exchange': 1.0}  # objective used in the flux simulation (dict) - it was minimized
    op = 'mean'  # aggregate flux values of each subsystem by mean
    htmp_log = True  # logaritmize values for heatmap
    trsf = 'none'  # do not apply any standartization or normalization to values in heatmap
    per_tissue = False  # run for cell lines
    gr_const = True  # it was used constraint in biomass in final simulations
    cl_spec = True  # it was used cell line specific total DNA reaction in final simulations
    PathAnalysis.analyze_pathway_protein(mth_corr_fld_path=EC_FLX_FLD, algo=algo, obj_id=obj_id, constr_ex=constr_ex,
                            gen_md_pth=GEN_MD_PTH, gen_gecko_pth=GEN_GECKO_PTH, methlt_fl_rc=METHLT_FL_RC, with_tsk=with_tsk,
                            htmp_log=htmp_log, per_tissue=per_tissue, op=op, trsf=trsf,
                            cell_width=0.7, cell_height=0.5, fn=12, xlab_siz=17, ylab_siz=16, xlab_rot=90,
                            xcbar=2.2, ycbar=0.7, cbar_width=0.05, cbar_height=0.09, xheat=1, yheat=1, heat_wf=1,
                            heat_hf=1, rden_xf=1.36, cden_xf=1, rden_yf=1, cden_yf=0.92, rden_wf=1.1, cden_wf=1, rden_hf=1,
                            cden_hf=0.55, cbar_lab_siz=16, rclust=False, cclust=True,
                            gr_const=gr_const, cl_spec=cl_spec)

