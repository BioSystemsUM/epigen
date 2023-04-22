home_path = '/home/tbarata'
import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])

import pandas as pd
import matplotlib.pyplot as plt
from cobra.io import load_matlab_model
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.decomposition import PCA
from faker import Factory
from pathway_analysis import PathAnalysis
import test_mds



#### PCA ####
def dopca(data, pth):
    '''
    - does PCA by tissue
    :param data: dataframe where column labels are cell lines with tissue name included,
                 i.e. of the type 'UO31_KIDNEY'.
    :param pth: path to file where to save
    '''
    # prepare data for PCA:
    data.columns = ['_'.join(col.split('_')[1:]) for col in data.columns]
    data = data.T
    data.reset_index(inplace=True)
    data.rename(columns={'index': 'Tissue'}, inplace=True)
    dataval = data.iloc[:, 1:]
    # scale features/columns/reactions/genes to normally distributed values with a mean of 0 and std of 1 (centered to have zero mean):
    x = StandardScaler().fit_transform(dataval)
    # get dataframe with pc (principal component) values for all samples:
    pca = PCA(n_components=2)  # pca with 2 PCs
    pc_best = pca.fit_transform(x)
    pcdf = pd.concat([pd.DataFrame(data=pc_best, columns=['PC1', 'PC2']), data[['Tissue']]], axis=1)
    # get amount of variance each pc holds:
    pc_var = pca.explained_variance_ratio_
    # graph:
    plt.figure(figsize=(12, 7))
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    plt.xlabel('PC1({}%)'.format(round(pc_var[0] * 100, 2)), fontsize=20, labelpad=10)
    plt.ylabel('PC2({}%)'.format(round(pc_var[1] * 100, 2)), fontsize=20, labelpad=0.3)
    fake = Factory.create()
    clasf = set(pcdf['Tissue'])
    clLst = [fake.hex_color() for n in range(0, len(clasf))]
    for g, c in zip(clasf, clLst):
        keep_ind = pcdf['Tissue'] == g
        plt.scatter(pcdf.loc[keep_ind, 'PC1'], pcdf.loc[keep_ind, 'PC2'], color=c, s=50 * 2)
    plt.legend(clasf, bbox_to_anchor=(1.04, 1), loc='upper left', prop={'size': 15})
    plt.tight_layout()
    plt.savefig(pth, dpi=300)
    plt.close('all')

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

    ##################################################
    ## fluxes - tissue (using 31 cell lines)
    # heatmap mean, log10, No norm min-max on rows
    # boxplot mean
    # TISSUE
    EC_FLX_FLD = 'results/ecGEM_simul_human1'
    METHLT_FL_RC = 'data/DNAmeth_reactions.json'
    GEN_MD_PTH = 'data/models/Human-GEM.yml'
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    with_tsk=False # tasks were not applied during reconstruction
    algo = 'init' # algorithm used in reconstruction was tINIT
    constr_ex=False # no constraints were applied to exchange fluxes of exometabolites
    obj_id = {'prot_pool_exchange': 1.0} # objective used in the flux simulation (dict) - it was minimized
    op = 'mean' # aggregate flux values of each subsystem by mean
    htmp_log=True # logaritmize values for heatmap
    trsf='none' # do not apply any standartization or normalization to values in heatmap
    per_tissue = True
    PathAnalysis.analyze_pathway_flux(mth_corr_fld_path=EC_FLX_FLD, algo=algo, obj_id=obj_id, constr_ex=constr_ex, gen_md_pth=GEN_MD_PTH, gen_gecko_pth=GEN_GECKO_PTH,
                         methlt_fl_rc=METHLT_FL_RC, with_tsk=with_tsk, htmp_log=htmp_log, per_tissue=per_tissue, op=op, trsf=trsf,
                         cell_width=0.7, cell_height=0.5, fn=12, xlab_siz=17, ylab_siz=16, xlab_rot=90,
                         xcbar=2.2, ycbar=0.65, cbar_width=0.05, cbar_height=0.09, xheat=1, yheat=1, heat_wf=1, heat_hf=1,
                         rden_xf=1.36, cden_xf=1, rden_yf=1, cden_yf=0.92, rden_wf=1.1, cden_wf=1, rden_hf=1, cden_hf=0.55,
                         cbar_lab_siz=16, rclust=False, cclust=True)

    ## protein usage - tissue (using 31 cell lines)
    # heatmap mean, log10, No norm min-max on rows
    # boxplot mean
    # TISSUE
    # EC_FLX_FLD = 'results/ecGEM_simul_human1'
    # METHLT_FL_RC = 'data/DNAmeth_reactions.json'
    # GEN_MD_PTH = 'data/models/Human-GEM.yml'
    # GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    # with_tsk=False # tasks were not applied during reconstruction
    # algo = 'init' # algorithm used in reconstruction was tINIT
    # constr_ex=False # no constraints were applied to exchange fluxes of exometabolites
    # obj_id = {'prot_pool_exchange': 1.0} # objective used in the flux simulation (dict) - it was minimized
    # op = 'mean' # aggregate flux values of each subsystem by mean
    # htmp_log=True # logaritmize values for heatmap
    # trsf='none' # do not apply any standartization or normalization to values in heatmap
    # per_tissue = True
    # PathAnalysis.analyze_pathway_protein(mth_corr_fld_path=EC_FLX_FLD, algo=algo, obj_id=obj_id, constr_ex=constr_ex, gen_md_pth=GEN_MD_PTH, gen_gecko_pth=GEN_GECKO_PTH,
    #                         methlt_fl_rc=METHLT_FL_RC, with_tsk=with_tsk, htmp_log=htmp_log, per_tissue=per_tissue, op=op, trsf=trsf,
    #                         cell_width=0.7, cell_height=0.5, fn=12, xlab_siz=17, ylab_siz=16, xlab_rot=90,
    #                         xcbar=2.2, ycbar=0.7, cbar_width=0.05, cbar_height=0.09, xheat=1, yheat=1, heat_wf=1, heat_hf=1,
    #                         rden_xf=1.36, cden_xf=1, rden_yf=1, cden_yf=0.92, rden_wf=1.1, cden_wf=1, rden_hf=1, cden_hf=0.55,
    #                         cbar_lab_siz=16, rclust=False, cclust=True)

    ## PCA:
    TRANSCRIPTOMICS = 'data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.txt'
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    EC_FLX_FLD = 'support/ecGEM_simul_human1'
    mth_corr_fld_path=EC_FLX_FLD
    gen_gecko_pth=GEN_GECKO_PTH
    transcriptomics=TRANSCRIPTOMICS
    # data = rec.reaction_scores
    # data = rec.gene_scores
    # data = rec.gene_data
    # get flux dataframe:
    if with_tsk:
        mth_corr_fld_path = os.path.join(mth_corr_fld_path, algo, 'including_tsks')
    else:
        mth_corr_fld_path = os.path.join(mth_corr_fld_path, algo)
    if constr_ex:
        mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'w_flux_constr')
    else:
        mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'no_constr')
    file_lst = test_mds.TestMds.create_file_list(path=mth_corr_fld_path, incl_wrd='fluxes')
    file_lst = [f for f in file_lst if list(obj_id.keys())[0] in f]
    if type(constr_ex) == list:
        file_lst = [f for f in file_lst if
                    '-'.join(constr_ex) in f]  # use only the files corresponding to the constraints we are testing
    file_pth = file_lst[0]
    print(file_pth)
    df = pd.read_csv(file_pth, sep='\t', index_col=0)
    # load generic gecko model:
    ms = load_matlab_model(gen_gecko_pth)
    # load transcriptomics data:
    data = pd.read_csv(transcriptomics, sep = '\t')
    data['gene_id'] = data.gene_id.apply(lambda x: x.split('.')[0])
    data.set_index('gene_id', inplace=True)
    data.drop('transcript_ids', axis=1, inplace=True)
    # select just the genes in the generic gecko model and the cell lines used in reconstructed models:
    transdata = data.loc[set(data.index).intersection(set([g.id for g in ms.genes])), set(df.columns).intersection(set(data.columns))]
    # do PCA for transcriptomics data:
    dopca(data=transdata, pth=os.path.join(mth_corr_fld_path, 'pca_genexp_mdgen.svg'))
    # do PCA for flux data:
    # df.drop('adaptbiomass', axis=0)
    dopca(data=df, pth=os.path.join(mth_corr_fld_path, 'pca_rc_flux.svg'))




    '''
    # get flux of lactate, O2 and glucose:
    ff = df.loc[list(rc_eval) + [el+'_REV' for el in rc_eval.keys()] + ['adaptbiomass'], :]
    for el in rc_eval: ff.loc[el + '_dif'] = np.abs(ff.loc[el] - ff.loc[el+'_REV'])
    ff = ff.loc[ff.index.str.contains('dif')]
    ff = ff.rename({k+'_dif': v for k,v in rc_eval.items()}, axis=0)
    ff = ff.T
    ff['Cell_lines'] = ff.index
    ff = pd.melt(ff, id_vars=['Cell_lines'], value_vars=['glucose', 'lactate', 'O2'], var_name='Metabolites', value_name='Flux')
    # ff['Simulation_type'] = [lt_lst[i]] * ff.shape[0]
    # ff_lst.append(ff)
    rc_eval = {'MAR09034': 'glucose', 'MAR09135': 'lactate', 'MAR09048': 'O2'}
    '''

