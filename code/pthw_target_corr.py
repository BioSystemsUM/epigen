home_path = '/home/tbarata'

import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])

import pandas as pd
from test_mds import TestMds
from scipy.stats import pearsonr, spearmanr
import copy
from pathway_analysis import PathAnalysis
import collections
import numpy as np
from cobra.io import load_matlab_model, load_yaml_model, read_sbml_model

def pthw_corr_meth_bmindp(fld_pth, fl_nm, methlt_dt_fld, smp_info, gen_gecko_pth, div_bm, **kwargs):
    '''
    - obtain the metabolic pathways which flux/protein usage correlated with overall DNA methylation independently of cell growth rate across the different cell lines
    - obtain reactions which flux correlated with overall DNA methylation independently of cell growth rate across the different cell lines
    - obtain enzymes which protein usage correlated with overall DNA methylation independently of cell growth rate across the different cell lines
    :param fld_pth: path to folder containing fluxes (mmol/gDW.h) or protein usage (mg/gDW)
                    for reactions in cell line-specific gecko models.
                    For the table with fluxes, the reactions excluded from the table are: 'draw_prot', 'prot_pool_exchange', and leg reactions corresponding to 'arm_' reactions (cause 'arm_' reactions flux is already representative of all the associated legs).
                    For the table with protein usage, the reactions excluded from the table are: 'draw_prot', 'prot_pool_exchange', 'arm_' reactions (cause those do not have an associated enzyme)
    :param fl_nm: name of file. has to contain either 'flx' or 'prot' depending on whether relates to fluxes or protein usage
    :param methlt_dt_fld: path to folder where file with experimental DNA methylation fluxes is.
    :param smp_info: path to file with cell line info
    :param gen_gecko_pth: path to generic gecko model
    :param div_bm: whether to obtain values independently of cell growth rate or not.
    :param kwargs - fl_nm_bm: path to dataframe with biomass flux. use when analyzing protein usage.
                  - ensembl_symbol_file: path to .json file where to retrieve or where to save (if the file was not created before) the mappings between ensembl gene id and gene symbol.
    '''
    for othr in ['spearman', 'pearson', 'common_butcoefpearson']:
        if not os.path.exists(os.path.join(fld_pth, othr)):
            os.makedirs(os.path.join(fld_pth, othr))
    # get dataframe with fluxes/protein usage for reactions in cell line-specific models
    # even reactions not present in cell specific model are in this dataframe (just have flux/protein usage of 0):
    simul_df = pd.read_csv(os.path.join(fld_pth, fl_nm), index_col=0, sep='\t')
    # obtain list of experimental methylation degree:
    mth_res_lst = TestMds.obtain_methylation(methlt_dt_fld=methlt_dt_fld, smp_info=smp_info)
    col_lst = [[col for col in simul_df.columns if col in exp_df.index] for exp_df in mth_res_lst]
    exp_df_lst = [exp_df[lst] for exp_df, lst in zip(mth_res_lst, col_lst)]
    # join experimental and simulated datasets:
    df = pd.concat([simul_df.T] + exp_df_lst, axis=1)
    d = {'Illumina_450K_methylation_Gene_average': 'Genes',
         'DNA_methylation_ratio': 'Global DNA methylation',
         'TSS1kb': 'Upstream of TSS',
         'cgi_CpG_clusters': 'CpG islands (clusters)',
         'enh_CpG_clusters': 'Enhancers (clusters)',
         'tss_CpG_clusters': 'TSS (clusters)'}
    df.rename(columns=d, inplace=True)
    exp_vars = [d[exp.name] for exp in exp_df_lst]
    sim_vars = [cl for cl in df.columns if (cl not in exp_vars)]
    # get only the simulated fluxes of all reactions (exclude real methylation levels), and all reactions except biomass:
    if 'flx' in fl_nm: ndf2 = df.drop('adaptbiomass', axis=1)
    else: ndf2 = copy.deepcopy(df)
    ndf2=ndf2.loc[:, set(sim_vars)-{'adaptbiomass'}]
    # get global DNA methylation and biomass:
    gb = df['Global DNA methylation'].drop('Subsystem', axis=0) # exclude 'Subsystem'
    if 'flx' in fl_nm:
        biomass=df['adaptbiomass'].drop('Subsystem') # exclude 'Subsystem'
    else:
        #df_w_bm = pd.read_csv(os.path.join(fld_pth, fl_nm_bm), index_col=0, sep='\t')
        df_w_bm = pd.read_csv(os.path.join(fld_pth, kwargs['fl_nm_bm']), index_col=0, sep='\t')
        biomass = df_w_bm.loc['adaptbiomass'].drop('Subsystem')
        ndf2.columns = pd.MultiIndex.from_arrays(ndf2.iloc[0:2].values)
        ndf2.drop(['Enzyme', 'Reaction'], axis=0, inplace=True)
        gb.drop(['Enzyme', 'Reaction'], axis=0, inplace=True)
    # corrected by (divide by) biomass flux or not:
    if div_bm:
        rc_by_bm = ndf2.drop('Subsystem', axis=0).divide(biomass, axis=0) # reaction flux/protein usage divided by biomass
    else:
        rc_by_bm = ndf2.drop('Subsystem', axis=0)
    # aggregate by subsystem:
    sbsm_by_bm_old = rc_by_bm.append(pd.DataFrame(ndf2.loc['Subsystem', :]).T)
    sbsm_by_bm = sbsm_by_bm_old.T.groupby('Subsystem').mean().T # each average substm flux/protein usage divided by biomass
    # calculate coefficient and P-value:
    # with spearman correlation:
    res_spm = sbsm_by_bm.apply(lambda x: spearmanr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
    res_spm.index = ['Coefficient', 'P-value']
    # with pearson correlation:
    res_prs = sbsm_by_bm.apply(lambda x: pearsonr(x, pd.DataFrame(gb)), axis=0) # NaN is when there is no variation
    res_prs.index = ['Coefficient', 'P-value']
    # select metabolic systems/reactions with p-value < 0.05:
    sig_sys_spm = res_spm.loc[:, res_spm.loc['P-value'] < 0.05].T
    sig_sys_prs = res_prs.loc[:, res_prs.loc['P-value'] < 0.05].T
    # sort by correlation coefficient:
    sig_sys_spm.sort_values(by='Coefficient', ascending=False, inplace=True)
    sig_sys_prs.sort_values(by='Coefficient', ascending=False, inplace=True)
    # get rid of weird subsystems:
    def remove_subst(f):
        mett = ['Isolated', 'Miscellaneous', 'Artificial reactions', 'Pool reactions', 'Exchange/demand reactions']
        for x in mett:
            if x in f.index:
                f = f.drop(index=x)
        return f
    sig_sys_spm = remove_subst(sig_sys_spm)
    sig_sys_prs = remove_subst(sig_sys_prs)
    spm_pos= sig_sys_spm.loc[sig_sys_spm['Coefficient'] > 0]
    spm_neg= sig_sys_spm.loc[sig_sys_spm['Coefficient'] < 0]
    prs_pos= sig_sys_prs.loc[sig_sys_prs['Coefficient'] > 0]
    prs_neg= sig_sys_prs.loc[sig_sys_prs['Coefficient'] < 0]
    common_spm_prs_pos = set(spm_pos.index).intersection(prs_pos.index)
    common_spm_prs_neg = set(spm_neg.index).intersection(prs_neg.index)
    sig_sys_common_pos = sig_sys_prs.loc[common_spm_prs_pos] # common to spearman and pearson, but coefficient from pearson
    sig_sys_common_neg = sig_sys_prs.loc[common_spm_prs_neg] # common to spearman and pearson, but coefficient from pearson
    sig_sys_common_pos.sort_values(by='Coefficient', ascending=False, inplace=True)
    sig_sys_common_neg.sort_values(by='Coefficient', ascending=False, inplace=True)
    # save results:
    if 'flx' in fl_nm: el = 'flx'
    elif 'prot' in fl_nm: el = 'prot'
    if div_bm: bmv = '_divbm_'
    else: bmv = ''
    sig_sys_spm.to_csv(os.path.join(fld_pth, 'spearman', f'corr_meth_{el}_substm{bmv}.tsv'), sep='\t')
    sig_sys_prs.to_csv(os.path.join(fld_pth, 'pearson', f'corr_meth_{el}_substm{bmv}.tsv'), sep='\t')
    sig_sys_common_pos.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_{el}_substm_pos{bmv}.tsv'), sep='\t')
    sig_sys_common_neg.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_{el}_substm_neg{bmv}.tsv'), sep='\t')
    if 'flx' in fl_nm:
        ## identify reactions of each subsystem which flux correlates with overall DNA methylation,
        ## independently of cell growth rates:
        # with spearman correlation:
        spm = sbsm_by_bm_old.iloc[:-1,:].apply(lambda x: spearmanr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
        spm.index = ['Coefficient', 'P-value']
        spm.loc['Subsystem'] = sbsm_by_bm_old.loc['Subsystem']
        # with pearson correlation:
        prs = sbsm_by_bm_old.iloc[:-1,:].apply(lambda x: pearsonr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
        prs.index = ['Coefficient', 'P-value']
        prs.loc['Subsystem'] = sbsm_by_bm_old.loc['Subsystem']
        # select metabolic reactions with p-value < 0.05:
        sig_spm = spm.loc[:, spm.loc['P-value'] < 0.05].T
        sig_prs = prs.loc[:, prs.loc['P-value'] < 0.05].T
        # get rid of weird subsystems:
        def remove_subst_cl(f):
            mett = ['Isolated', 'Miscellaneous', 'Artificial reactions', 'Pool reactions', 'Exchange/demand reactions']
            for x in mett:
                if x in f['Subsystem']:
                    f = f.drop(index=x)
            return f
        sig_spm = remove_subst_cl(sig_spm)
        sig_prs = remove_subst_cl(sig_prs)
        # intersect positively/negatively correlated of spearman and pearson methods:
        spm_pos = sig_spm.loc[sig_spm['Coefficient'] > 0]
        spm_neg = sig_spm.loc[sig_spm['Coefficient'] < 0]
        prs_pos = sig_prs.loc[sig_prs['Coefficient'] > 0]
        prs_neg = sig_prs.loc[sig_prs['Coefficient'] < 0]
        common_spm_prs_pos = set(spm_pos.index).intersection(prs_pos.index)
        common_spm_prs_neg = set(spm_neg.index).intersection(prs_neg.index)
        # common to spearman and pearson, but coefficient from pearson:
        sig_common_pos = sig_prs.loc[common_spm_prs_pos]
        sig_common_neg = sig_prs.loc[common_spm_prs_neg]
        sig_common_pos['Coefficient'] = sig_common_pos['Coefficient'].apply(lambda x: x[0])
        sig_common_neg['Coefficient'] = sig_common_neg['Coefficient'].apply(lambda x: x[0])
        sig_common_pos.sort_values(by=['Coefficient'], ascending=False, inplace=True)
        sig_common_neg.sort_values(by=['Coefficient'], ascending=True, inplace=True)
        # add reaction equations:
        ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
        mtid_nm = {mt.id: mt.name for mt in ms.metabolites}
        def spl_fc(x):
            rgs = ' + '.join([mtid_nm[rg.id]+f'[{rg.compartment}]' for rg in ms.reactions.get_by_id(x).reactants if 'prot_' not in rg.id])
            pds = list()
            for pd in ms.reactions.get_by_id(x).products:
                if 'pmet_' in pd.id:
                    for fpd in ms.reactions.get_by_id(pd.id[5:] + 'No1').products:
                        pds.append(mtid_nm[fpd.id] + f'[{fpd.compartment}]')
                else:
                    pds.append(mtid_nm[pd.id] + f'[{pd.compartment}]')
            pds = ' + '.join(pds)
            return ' --> '.join([rgs, pds])
        sig_common_pos['Reaction'] = sig_common_pos.index.to_series().apply(spl_fc)
        sig_common_neg['Reaction'] = sig_common_neg.index.to_series().apply(spl_fc)
        # save results:
        sig_common_pos.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_flx_rc_pos{bmv}.tsv'), sep='\t')
        sig_common_neg.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_flx_rc_neg{bmv}.tsv'), sep='\t')
    else:
        ## identify enzymes which flux correlates with overall DNA methylation,
        ## independently of cell growth rates:
        sbsm_by_bm_old.drop('Subsystem', axis=0, inplace=True)
        sbsm_by_bm_old.loc['Enzyme'] = [enz for enz, rc in sbsm_by_bm_old.columns]
        sbsm_divided_by_bm_old = sbsm_by_bm_old.T.groupby('Enzyme').sum().T  # for each enzyme get sum of ratio of that enzyme usage divided by cell growth across different reactions
        # with spearman correlation:
        spm = sbsm_divided_by_bm_old.apply(lambda x: spearmanr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
        spm.index = ['Coefficient', 'P-value']
        # with pearson correlation:
        prs = sbsm_divided_by_bm_old.apply(lambda x: pearsonr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
        prs.index = ['Coefficient', 'P-value']
        # select metabolic reactions with p-value < 0.05:
        sig_spm = spm.loc[:, spm.loc['P-value'] < 0.05].T
        sig_prs = prs.loc[:, prs.loc['P-value'] < 0.05].T
        # intersect positively/negatively correlated of spearman and pearson methods:
        spm_pos = sig_spm.loc[sig_spm['Coefficient'] > 0]
        spm_neg = sig_spm.loc[sig_spm['Coefficient'] < 0]
        prs_pos = sig_prs.loc[sig_prs['Coefficient'] > 0]
        prs_neg = sig_prs.loc[sig_prs['Coefficient'] < 0]
        common_spm_prs_pos = set(spm_pos.index).intersection(prs_pos.index)
        common_spm_prs_neg = set(spm_neg.index).intersection(prs_neg.index)
        # common to spearman and pearson, but coefficient from pearson:
        sig_common_pos = sig_prs.loc[common_spm_prs_pos]
        sig_common_neg = sig_prs.loc[common_spm_prs_neg]
        sig_common_pos['Coefficient'] = sig_common_pos['Coefficient'].apply(lambda x: x[0])
        sig_common_neg['Coefficient'] = sig_common_neg['Coefficient'].apply(lambda x: x[0])
        sig_common_pos.sort_values(by=['Coefficient'], ascending=False, inplace=True)
        sig_common_neg.sort_values(by=['Coefficient'], ascending=True, inplace=True)
        # add gene corresponding to enzyme:
        ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
        drw_rc = {r.id for r in ms.reactions if 'draw_' in r.id}
        prt_gene_d = {k: el.id for k, v in {mt.id: ms.reactions.get_by_id(rid).genes for rid in drw_rc for mt in ms.reactions.get_by_id(rid).metabolites if 'pool' not in mt.id}.items() for el in v}  # dictionary {enzyme: gene}
        sig_common_pos['Gene_id']=[prt_gene_d[el] for el in sig_common_pos.index]
        sig_common_neg['Gene_id']=[prt_gene_d[el] for el in sig_common_neg.index]
        # symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=ensembl_symbol_file)
        symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=kwargs['ensembl_symbol_file'])
        id_smb_pos = {id: smb for smb, id in symbol_ensembl if id in list(sig_common_pos['Gene_id'])}
        sig_common_pos['Symbol']=[id_smb_pos[el] for el in sig_common_pos['Gene_id']]
        id_smb_neg = {id: smb for smb, id in symbol_ensembl if id in list(sig_common_neg['Gene_id'])}
        sig_common_neg['Symbol'] = [id_smb_neg[el] for el in sig_common_neg['Gene_id']]
        # save results:
        sig_common_pos.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_prot_usage_pos{bmv}.tsv'), sep='\t')
        sig_common_neg.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_prot_usage_neg{bmv}.tsv'), sep='\t')

def experimental_corr_test(methlt_dt_fld, trscp_pth, mth_trscp_corrvalue, envcond, gen_gecko_pth, ensembl_symbol_file,
                           smp_info, gen_md_pth_sbml, gen_md_pth, fld_pth, div_growth):
    '''
    - identify metabolic subsystems over-represented in the list of reactions associated with:
     * genes which expression level correlate (positively/negatively) with overall DNA methylation level
     independently of cell growth rate across different cell lines
    - check whether certain genes are in the list of genes which expression correlates with overall methylation
    :param methlt_dt_fld: path to folder where file with experimental DNA methylation fluxes is.
    :param trscp_pth: path to file with gene expression data
    :param mth_trscp_corrvalue: whether to use direct or inverse correlation between overall methylation and gene expression
    :param envcond: path to dataframe with experimental growth rate
    :param gen_gecko_pth: path to generic gecko model with DNA (de)/methylation reactions, but NO subsystem information
    :param ensembl_symbol_file: path to .json file where to retrieve or where to save (if the file was not created before) the mappings between ensembl gene id and gene symbol.
    :param smp_info: path to file with information on cell lines
    :param gen_md_pth_sbml: path to traditional generic model with DNA methylation and demethylation reactions in sbml format
    :param gen_md_pth: path to generic traditional model withOUT DNA (de)/methylation reactions but with subsystem information, in yaml format
    :param fld_pth: path to folder where to save results
    :param div_growth: whether to divide overall methylation value by the experimental growth rate
    '''
    ### symbol- ensembl gene id in generic model:
    symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=ensembl_symbol_file)
    ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
    mdgenes = {g.id for g in ms.genes}
    # (symbol, ensembl) for all genes in the generic model:
    smbl_id_md = [(smb, id) for smb, id in symbol_ensembl if id in mdgenes]  # 2515 genes
    if len([id for s, id in smbl_id_md]) != len({id for s, id in smbl_id_md}):
        print("lookout! the same ensembl id has more than one symbol")
    elif len([s for s, id in smbl_id_md]) != len({s for s, id in smbl_id_md}):
        print("lookout! the same symbol has more than one ensembl id")
    ### genes which expression correlates with overall DNA methylation:
    # load gene expression matrix:
    trscp = pd.read_csv(trscp_pth, sep='\t')
    trscp.index = trscp['gene_id'].apply(lambda x: x.split('.')[0])
    trscp.drop(['gene_id', 'transcript_ids'], axis=1, inplace=True)
    # get overall DNA methylation degree (load again, cause meanwhile was changed):
    mth_res_lst = TestMds.obtain_methylation(methlt_dt_fld=methlt_dt_fld, smp_info=smp_info)
    for el in mth_res_lst:
        if el.name == 'DNA_methylation_ratio': glb_mth = el
    glb_mth.name = 'Global DNA methylation'
    # if required, divide gene expression and overall methylation by experimental growth rate:
    if div_growth:
        grf = pd.read_excel(envcond, sheet_name='GrowthRate', index_col=0)
        grf['mean'] = grf.mean(axis=1)
        biomass = grf['mean']
        inters = list(set(biomass.index).intersection(glb_mth.index).intersection(trscp.columns))
        glb_mth = glb_mth[inters]
        biomass = biomass[inters]
        trscp = trscp[inters]
        glb_mth = glb_mth.divide(biomass)
    # get correlations:
    cl_ord = list(set(trscp.columns).intersection(set(glb_mth.index)))  # 32 common cell lines to both datasets
    trscp = trscp[cl_ord]
    glb_mth = glb_mth.loc[cl_ord]
    # with spearman method:
    res_trscp_spm = trscp.apply(lambda x: spearmanr(x, glb_mth), axis=1)  # NaN is when there is no variation
    corr_trscp_mth_spm = pd.DataFrame(map(lambda x: (x.correlation, x.pvalue), res_trscp_spm), index=trscp.index, columns=['Coefficient', 'P-value'])
    # select genes with a significant correlation in spearman method:
    corr_trscp_mth_spm_sig = corr_trscp_mth_spm[corr_trscp_mth_spm['P-value'] < 0.05]
    # with pearson method:
    res_trscp_prs = trscp.apply(lambda x: pearsonr(x, glb_mth), axis=1)  # NaN is when there is no variation
    corr_trscp_mth_prs = pd.DataFrame(map(lambda x: (x[0], x[1]), res_trscp_prs), index=trscp.index, columns=['Coefficient', 'P-value'])
    # select genes with a significant correlation in pearson method:
    corr_trscp_mth_prs_sig = corr_trscp_mth_prs[corr_trscp_mth_prs['P-value'] < 0.05]
    # from those above select those with direct or inverse correlation:
    if mth_trscp_corrvalue == 'direct':
        corr_sig_spm_trscp = corr_trscp_mth_spm_sig[corr_trscp_mth_spm_sig['Coefficient'] > 0]
        corr_sig_prs_trscp = corr_trscp_mth_prs_sig[corr_trscp_mth_prs_sig['Coefficient'] > 0]
    elif mth_trscp_corrvalue == 'inverse':
        corr_sig_spm_trscp = corr_trscp_mth_spm_sig[corr_trscp_mth_spm_sig['Coefficient'] < 0]
        corr_sig_prs_trscp = corr_trscp_mth_prs_sig[corr_trscp_mth_prs_sig['Coefficient'] < 0]
    # select genes with significant correlation (direct/inverse) using both spearman and pearson and report their pearson p-values:
    common_sig = list(set(corr_sig_spm_trscp.index).intersection(set(corr_sig_prs_trscp.index)))
    corr_trscp_mth = corr_trscp_mth_prs_sig.loc[common_sig]
    # get column with gene name:
    corr_trscp_mth['gene_name'] = list(map(lambda x: x.split('_')[0], list(corr_trscp_mth.index)))
    # from above-mentioned genes get those with a significant correlation:
    corr_same_sig = corr_trscp_mth[corr_trscp_mth['P-value'] < 0.05]
    # get the significant correlated genes:
    mth_trscp = set(corr_same_sig['gene_name'])
    mth_trscp_id = {id for s, id in smbl_id_md if id in mth_trscp}  # genes in the generic metabolic model only
    if len(mth_trscp) > 0 and len(mth_trscp_id) == 0:
        print('None of the genes identified were metabolic genes (were in the generic model)')
        return
    # Note: many genes lost (len(mth_trscp_id) is 843) are because they are not part of model (not metabolic),
    # not because of conversion symbol to ensembl as len(mdgenes) is 2520 and len(smbl_id_md) is 2515,
    # so just 5 ensembl ids of the model didn't have corresponding symbol
    ## Get Reactions corresponding to those genes:
    md = read_sbml_model(gen_md_pth_sbml)  # generic traditional model with DNA meth reactions
    rc_m = [r.id for r in md.reactions]
    rc_gene_d = {k: [el.id for el in v] for k, v in {rid: md.reactions.get_by_id(rid).genes for rid in rc_m}.items()}
    # dict with gene vs list of reactions for all genes in traditional generic model:
    gene_rc_d = dict()
    for r, gl in rc_gene_d.items():
        for g in gl:
            if g not in gene_rc_d:
                gene_rc_d[g] = [r]
            elif g in gene_rc_d:
                gene_rc_d[g].append(r)
    ### Check over-representation of each metabolic subsystem in the list of reactions
    ### of genes which expression correlates with overall DNA methylation:
    # dict with gene vs list of reactions for aforementioned genes:
    mth_trscp_gene_rc = {g: gene_rc_d[g] for g in list(mth_trscp_id)}
    # dict with gene vs list of subsystems for all genes in traditional generic model:
    m = load_yaml_model(gen_md_pth)  # generic traditional model with subsystems
    rc_sb_d = {r.id: r.subsystem[0] for r in m.reactions}
    extra_rc = {r.id for r in md.reactions} - {r.id for r in m.reactions}
    excp = ['prodDNAtot', 'adaptbiomass']  # artificial reactions - pseudo-reactions
    extra_rc_sb = {r: 'Transport reactions' if 'transp' in r else 'Artificial reactions' if r in excp else 'Dna (de)/methylation' for r in extra_rc}
    rc_sb_d.update(extra_rc_sb)
    gene_sb_d = {g: [rc_sb_d[r] for r in rl] for g, rl in gene_rc_d.items()}
    # dict with gene vs list of subsystems for aforementioned genes:
    mth_trscp_gene_sb = {g: gene_sb_d[g] for g, r in mth_trscp_gene_rc.items()}
    # subsystems associated with those genes:
    mth_trscp_sb = {sb for sb_l in mth_trscp_gene_sb.values() for sb in sb_l}
    # do hypergemotric test:
    mth_trscp_hyp = PathAnalysis.hypergemtest(sample_gene_sb=mth_trscp_gene_sb, population_gene_sb=gene_sb_d, sample_sb=mth_trscp_sb)
    # remove weird subsystems:
    def remove_subst(f):
        mett = ['Isolated', 'Miscellaneous', 'Artificial reactions', 'Pool reactions', 'Exchange/demand reactions']
        for x in mett:
            if x in f.index:
                f = f.drop(index=x)
        return f
    mth_trscp_hyp = remove_subst(mth_trscp_hyp)
    # save result:
    fld1 = os.path.join(fld_pth, 'experimental', 'overall_mth')
    if div_growth: bm = '_div_bm'
    else: bm = ''
    nm1 = 'overall.mth.vs.gene.exp.' + mth_trscp_corrvalue + '.common.pearson.spearman' + bm + '.tsv'
    if not os.path.exists(fld1):
        os.makedirs(fld1)
    final_pth_3 = os.path.join(fld1, nm1)
    mth_trscp_hyp.to_csv(final_pth_3, sep='\t')
    ### check whether certain genes are in the list of genes which expression correlates with overall DNA methylation:
    chk_gn_dct = {'ENSG00000168906': 'MAT2A',
                  'ENSG00000038274': 'MAT2B',
                  'ENSG00000127804': 'METTL16',
                  'ENSG00000139160': 'ETFBKMT/METTL20',
                  'ENSG00000198911': 'SREBF2',
                  'ENSG00000280496': 'ACA43/SNORA17B',
                  'ENSG00000276161': 'ACA43/SNORA17B'}
    have_gn = [(idn, smb) for idn, smb in chk_gn_dct.items() if idn in mth_trscp]
    if have_gn:
        print(f'Of the genes we want to check, those which expression correlates with overall methylation are: {have_gn}')
    else:
        print('None of the genes we want to check, had expression correlating with overall methylation')


def experimental_corr(meth_fl_nm, methlt_dt_fld, smp_info, ensembl_symbol_file, gen_gecko_pth, gen_md_pth,
                      gen_md_pth_sbml, all_prom_same_behavior, mth_mth_corrvalue, fld_pth,
                      envcond, div_growth, mth_exp_corrvalue):
    '''
    - identify metabolic subsystems over-represented in the list of reactions associated with:
     * genes which methylation level correlate (directly/inversely) with overall DNA methylation level across different cell lines
     * those same genes that also are correlated (directly/inversely) with gene expression.
    :param meth_fl_nm: path to file with experimental methylation we want to use. can be from TSS or other genomic regions.
    :param methlt_dt_fld: path to folder where file with experimental DNA methylation fluxes is.
    :param smp_info: path to file with information on cell lines
    :param ensembl_symbol_file: path to .json file where to retrieve or where to save (if the file was not created before) the mappings between ensembl gene id and gene symbol.
    :param gen_gecko_pth: path to generic gecko model with DNA (de)/methylation reactions, but NO subsystem information
    :param gen_md_pth: path to generic traditional model withOUT DNA (de)/methylation reactions but with subsystem information, in yaml format
    :param gen_md_pth_sbml: path to traditional generic model with DNA methylation and demethylation reactions in sbml format
    :param all_prom_same_behavior: whether to use only genes where all promoters are either directly or inversely correlated
                                   or instead consider any gene with at least one promoter with the behavior we want (direct/inversely correlated)
    :param mth_mth_corrvalue: whether to use direct or inverse correlation between overall methylation and gene methylation
    :param fld_pth: path to folder where to save results
    :param envcond: path to dataframe with experimental growth rate
    :param div_growth: whether to divide gene methylation value by the experimental growth rate
    :param mth_exp_corrvalue: whether to use direct or inverse correlation between gene methylation and gene expression
    '''
    ## symbol- ensembl gene id in generic model:
    symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=ensembl_symbol_file)
    ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
    mdgenes = {g.id for g in ms.genes}
    # (symbol, ensembl) for all genes in the generic model:
    smbl_id_md = [(smb, id) for smb, id in symbol_ensembl if id in mdgenes]  # 2515 genes
    if len([id for s, id in smbl_id_md]) != len({id for s, id in smbl_id_md}):
        print("lookout! the same ensembl id has more than one symbol")
    elif len([s for s, id in smbl_id_md]) != len({s for s, id in smbl_id_md}):
        print("lookout! the same symbol has more than one ensembl id")
    ### Genes which methylation level significantly correlate with overall methylation level:
    # get overall DNA methylation degree:
    mth_res_lst = TestMds.obtain_methylation(methlt_dt_fld=methlt_dt_fld, smp_info=smp_info)
    for el in mth_res_lst:
        if el.name == 'DNA_methylation_ratio': glb_mth = el
    glb_mth.name = 'Global DNA methylation'
    # obtain experimental methylation degree of 1000bp upstream of TSS:
    raw_df = pd.read_csv(os.path.join(methlt_dt_fld, meth_fl_nm), sep='\t', index_col=0, na_values=['    NaN', '     NA']).iloc[:-1, :]
    meth_df = raw_df.iloc[:, 2:]
    mean_mth = meth_df.apply(lambda row: row.dropna().median(), axis=1)
    meth_df.apply(lambda row: row.fillna(value=mean_mth[row.name], inplace=True), axis=1)
    # if required, divide gene methylation and overall methylation by experimental growth rate:
    if div_growth:
        grf = pd.read_excel(envcond, sheet_name='GrowthRate', index_col=0)
        grf['mean'] = grf.mean(axis=1)
        biomass = grf['mean']
        inters = list(set(biomass.index).intersection(glb_mth.index).intersection(meth_df.columns))
        glb_mth = glb_mth[inters]
        biomass = biomass[inters]
        meth_df = meth_df[inters]
        glb_mth = glb_mth.divide(biomass)
        meth_df=meth_df.divide(biomass, axis=1)
    # get correlations:
    cl_ord = list(set(meth_df.columns).intersection(set(glb_mth.index)))  # 32 common cell lines to both datasets
    meth_df = meth_df[cl_ord]
    glb_mth = glb_mth.loc[cl_ord]
    # with spearman method:
    res_mth_spm = meth_df.apply(lambda x: spearmanr(x, glb_mth), axis=1)  # NaN is when there is no variation
    corr_meth_mth_spm = pd.DataFrame(map(lambda x: (x.correlation, x.pvalue), res_mth_spm), index=meth_df.index, columns=['Coefficient', 'P-value'])
    # select promoters with a significant correlation in spearman method:
    corr_meth_mth_spm_sig = corr_meth_mth_spm[corr_meth_mth_spm['P-value'] < 0.05]
    # with pearson method:
    res_mth_prs = meth_df.apply(lambda x: pearsonr(x, glb_mth), axis=1)  # NaN is when there is no variation
    corr_meth_mth_prs = pd.DataFrame(map(lambda x: (x[0], x[1]), res_mth_prs), index=meth_df.index, columns=['Coefficient', 'P-value'])
    # select promoters with a significant correlation in pearson method:
    corr_meth_mth_prs_sig = corr_meth_mth_prs[corr_meth_mth_prs['P-value'] < 0.05]
    # select promoters with significant correlation using both spearman and pearson and report their pearson p-values:
    common_sig = list(set(corr_meth_mth_spm_sig.index).intersection(set(corr_meth_mth_prs_sig.index)))
    corr_meth_mth = corr_meth_mth_prs_sig.loc[common_sig]
    # get column with gene name:
    corr_meth_mth['gene_name'] = list(map(lambda x: x.split('_')[0], list(corr_meth_mth.index)))
    if all_prom_same_behavior:
        # select promoters of genes which all promoters have same signal (either are inversely or directly) correlated:
        total_prm = dict(collections.Counter(corr_meth_mth['gene_name']))
        to_keep = list()
        for g in total_prm:
            sbs = corr_meth_mth.loc[corr_meth_mth['gene_name'] == g, 'Coefficient']
            if np.sum(sbs < 0) == total_prm[g] or np.sum(sbs > 0) == total_prm[g]:
                to_keep.append(g)
        corr_same_sig_meth = corr_meth_mth[corr_meth_mth['gene_name'].isin(to_keep)]
    else:
        corr_same_sig_meth = corr_meth_mth
    # from those above select those with direct or inverse correlation:
    if mth_mth_corrvalue == 'direct':
        corr_same_sig_meth = corr_same_sig_meth[corr_same_sig_meth['Coefficient'] > 0]
    elif mth_mth_corrvalue == 'inverse':
        corr_same_sig_meth = corr_same_sig_meth[corr_same_sig_meth['Coefficient'] < 0]
    # get the significant correlated genes:
    mth_mth = set(corr_same_sig_meth['gene_name'])
    mth_mth_id = {id for s, id in smbl_id_md if s in mth_mth} # genes in the generic metabolic model only
    # Note: many genes lost are because they are not part of model (not metabolic),
    # not because of conversion symbol to ensembl as len(mdgenes) is 2520 and len(smbl_id_md) is 2515,
    # so just 5 ensembl ids of the model didn't have corresponding symbol
    ## Get Reactions corresponding to those genes:
    md = read_sbml_model(gen_md_pth_sbml) # generic traditional model with DNA meth reactions
    rc_m = [r.id for r in md.reactions]
    rc_gene_d = {k: [el.id for el in v] for k,v in {rid: md.reactions.get_by_id(rid).genes for rid in rc_m}.items()}
    # dict with gene vs list of reactions for all genes in traditional generic model:
    gene_rc_d = dict()
    for r, gl in rc_gene_d.items():
        for g in gl:
            if g not in gene_rc_d:
                gene_rc_d[g] = [r]
            elif g in gene_rc_d:
                gene_rc_d[g].append(r)
    # dict with gene vs list of reactions for 'mth_mth_id' (genes which methylation correlates with overal methylation):
    mth_mth_gene_rc = {g: gene_rc_d[g] for g in list(mth_mth_id)}
    ## check over-representation of each metabolic subsystem in the list of reactions:
    # dict with gene vs list of subsystems for all genes in traditional generic model:
    m = load_yaml_model(gen_md_pth) # generic traditional model with subsystems
    rc_sb_d = {r.id: r.subsystem[0] for r in m.reactions}
    extra_rc = {r.id for r in md.reactions} - {r.id for r in m.reactions}
    excp = ['prodDNAtot', 'adaptbiomass']  # artificial reactions - pseudo-reactions
    extra_rc_sb = {r: 'Transport reactions' if 'transp' in r else 'Artificial reactions' if r in excp else 'Dna (de)/methylation'for r in extra_rc}
    rc_sb_d.update(extra_rc_sb)
    gene_sb_d={g: [rc_sb_d[r] for r in rl] for g, rl in gene_rc_d.items()}
    # dict with gene vs list of subsystems for 'mth_mth_id' (genes which methylation correlates with overall methylation):
    mth_mth_gene_sb = {g: gene_sb_d[g] for g, r in mth_mth_gene_rc.items()}
    # subsystems associated with 'mth_mth_id' genes:
    mth_mth_sb = {sb for sb_l in mth_mth_gene_sb.values() for sb in sb_l}
    # do hypergeometric test:
    mth_mth_hyp=PathAnalysis.hypergemtest(sample_gene_sb=mth_mth_gene_sb, population_gene_sb=gene_sb_d, sample_sb=mth_mth_sb)
    # remove weird subsystems:
    def remove_subst(f):
        mett = ['Isolated', 'Miscellaneous', 'Artificial reactions', 'Pool reactions', 'Exchange/demand reactions']
        for x in mett:
            if x in f.index:
                f = f.drop(index=x)
        return f
    mth_mth_hyp = remove_subst(mth_mth_hyp)
    # save result:
    if all_prom_same_behavior: prom = 'promoters.all.same'
    else: prom = 'promoter.any'
    if 'TSS1kb' in meth_fl_nm: mth_rg = 'gene_prom' # type of genomic region being methylated that is associated to each gene
    elif 'tss' in meth_fl_nm and 'clusters' in meth_fl_nm: mth_rg = 'tss_clust'
    fld1 = os.path.join(fld_pth, 'experimental', 'gene_overall_mth', mth_rg)
    if div_growth: bm = '_div_bm_'
    else: bm = ''
    nm1 = prom + '_' + mth_mth_corrvalue + '.common.pearson.spearman' + bm +'.tsv'
    if not os.path.exists(fld1):
        os.makedirs(fld1)
    final_pth_1 = os.path.join(fld1, nm1)
    mth_mth_hyp.to_csv(final_pth_1, sep='\t')

    ### Genes which methylation is correlated with gene expression - this is always Pearson Correlation (provided by study):
    mth_exp_corr = pd.read_excel(os.path.join(methlt_dt_fld, 'meth_exp_corr.xlsx'), sheet_name='meth_exp_corr', index_col=0)
    mth_exp_corr['gene_name'] = list(map(lambda x: x.split('_')[0], list(mth_exp_corr.index)))
    # get p-value using n data points and correlation score:
    mth_exp_corr['P-value'] = mth_exp_corr.apply(lambda x: PathAnalysis.r2p(r=x['Pearson_correlation'], N=x['n_datapoints']), axis=1)
    if all_prom_same_behavior:
        # select promoters of genes which all promoters have same signal (either are inversely or directly) correlated:
        total_prm = dict(collections.Counter(mth_exp_corr['gene_name']))
        to_keep = list()
        for g in total_prm:
            sbs = mth_exp_corr.loc[mth_exp_corr['gene_name'] == g, 'Pearson_correlation']
            if np.sum(sbs < 0) == total_prm[g] or np.sum(sbs > 0) == total_prm[g]:
                to_keep.append(g)
        corr_same = mth_exp_corr[mth_exp_corr['gene_name'].isin(to_keep)]
        # from above-mentioned promoters get those with a significant correlation:
        corr_same_sig = corr_same[corr_same['P-value'] < 0.05]
    else:
        # from above-mentioned promoters get those with a significant correlation:
        corr_same_sig = mth_exp_corr[mth_exp_corr['P-value'] < 0.05]
    # from those above select those with inverse or direct correlation:
    if mth_exp_corrvalue == 'direct':
        corr_same_sig = corr_same_sig.loc[corr_same_sig['Pearson_correlation'] > 0]
    elif mth_exp_corrvalue == 'inverse':
        corr_same_sig = corr_same_sig.loc[corr_same_sig['Pearson_correlation'] < 0]
    # get the significant correlated genes:
    mth_exp = set(corr_same_sig['gene_name'])
    mth_exp_id = {id for s, id in smbl_id_md if s in mth_exp}
    # Note: many genes lost are because they are not part of model (not metabolic),
    # not because of conversion symbol to ensembl as len(mdgenes) is 2520 and len(smbl_id_md) is 2515,
    # so just 5 ensembl ids of the model didn't have corresponding symbol
    ### Intersect genes which methylation correlates with overall methylation with genes which methylation correlates with gene expression:
    intersected_id = mth_mth_id.intersection(mth_exp_id) # this is just for genes in the metabolic model
    inters_all = mth_mth.intersection(mth_exp)
    ### Check over-representation of each metabolic subsystem in the list of reactions
    ### of the intersected genes:
    # dict with gene vs list of reactions for aforementioned 'intersected' genes:
    intersected_gene_rc = {g: gene_rc_d[g] for g in list(intersected_id)}
    # dict with gene vs list of subsystems for aforementioned 'intersected' genes:
    intersected_gene_sb = {g: gene_sb_d[g] for g, r in intersected_gene_rc.items()}
    # subsystems associated with 'intersected' genes:
    intersected_sb = {sb for sb_l in intersected_gene_sb.values() for sb in sb_l}
    # do hypergemotric test:
    intersected_hyp=PathAnalysis.hypergemtest(sample_gene_sb=intersected_gene_sb, population_gene_sb=gene_sb_d, sample_sb=intersected_sb)
    # remove weird subsystems:
    intersected_hyp = remove_subst(intersected_hyp)
    # save result:
    fld2 = os.path.join(fld_pth, 'experimental', 'gene_overall_mth_AND_gene_exp_mth', mth_rg)
    if div_growth: bm = '_div_bm_'
    else: bm = ''
    nm2 = prom + '_1st.' + mth_mth_corrvalue + '.common.pearson.spearman' + '_2nd.' + mth_exp_corrvalue + '.pearson'+ bm +'.tsv'
    if not os.path.exists(fld2):
        os.makedirs(fld2)
    final_pth_2 = os.path.join(fld2, nm2)
    intersected_hyp.to_csv(final_pth_2, sep='\t')



'''
    ### genes which expression correlates with overall DNA methylation:
    # lead gene expression matrix:
    # mth_trscp_corrvalue= 'direct'
    trscp_pth ='data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.txt'
    trscp = pd.read_csv(trscp_pth, sep='\t')
    trscp.index = trscp['gene_id'].apply(lambda x: x.split('.')[0])
    trscp.drop(['gene_id', 'transcript_ids'], axis=1, inplace=True)
    # get overall DNA methylation degree (load again, cause meanwhile was changed):
    mth_res_lst = TestMds.obtain_methylation(methlt_dt_fld=methlt_dt_fld, smp_info=smp_info)
    for el in mth_res_lst:
        if el.name == 'DNA_methylation_ratio': glb_mth = el
    glb_mth.name = 'Global DNA methylation'
    # get correlations:
    cl_ord = list(set(trscp.columns).intersection(set(glb_mth.index)))  # 32 common cell lines to both datasets
    trscp = trscp[cl_ord]
    glb_mth = glb_mth.loc[cl_ord]
    # with spearman method:
    res_trscp_spm = trscp.apply(lambda x: spearmanr(x, glb_mth), axis=1)  # NaN is when there is no variation
    corr_meth_trscp_spm = pd.DataFrame(map(lambda x: (x.correlation, x.pvalue), res_trscp_spm), index=trscp.index, columns=['Coefficient', 'P-value'])
    # select promoters with a significant correlation in spearman method:
    # corr_meth_trscp_spm = corr_meth_trscp_spm[corr_meth_trscp_spm['P-value'] < 0.05]
    # from those above select those with direct or inverse correlation:
    if mth_trscp_corrvalue == 'direct':
        corr_same_trscp_spm_pos = corr_meth_trscp_spm[corr_meth_trscp_spm['Coefficient'] > 0]
    elif mth_trscp_corrvalue == 'inverse':
        corr_same_trscp_spm_neg = corr_meth_trscp_spm[corr_meth_trscp_spm['Coefficient'] < 0]
    'ENSG00000198911' in corr_same_trscp_spm_pos.index
    'ENSG00000198911' in corr_same_trscp_spm_neg.index
    corr_same_trscp_spm_neg.loc['ENSG00000198911']
    'ENSG00000198911' in trscp.index
    # with pearson method:
    res_trscp_prs = trscp.apply(lambda x: pearsonr(x, glb_mth), axis=1)  # NaN is when there is no variation
    corr_meth_trscp_prs = pd.DataFrame(map(lambda x: (x[0], x[1]), res_trscp_prs), index=trscp.index, columns=['Coefficient', 'P-value'])
    # select promoters with a significant correlation in pearson method:
    # corr_meth_trscp_prs = corr_meth_trscp_prs[corr_meth_trscp_prs['P-value'] < 0.05]
    # from those above select those with direct or inverse correlation:
    if mth_trscp_corrvalue == 'direct':
        corr_same_trscp_prs = corr_meth_trscp_prs[corr_meth_trscp_prs['Coefficient'] > 0]
    elif mth_trscp_corrvalue == 'inverse':
        corr_same_trscp_prs = corr_meth_trscp_prs[corr_meth_trscp_prs['Coefficient'] < 0]
    # select promoters with significant correlation using both spearman OR pearson and report their pearson p-values:
    common_mth_trscp = list(set(corr_same_trscp_spm.index).union(set(corr_same_trscp_prs.index)))
    # len(common_mth_trscp) # 27809
    # get the significant correlated genes:
    corr_meth_trscp_gnm = list(map(lambda x: [s for s,id in symbol_ensembl if x==id], common_mth_trscp))
    # len(corr_meth_trscp_gnm) # 27809
    mth_trscp = {el for l in corr_meth_trscp_gnm for el in l}
    # len(mth_trscp) # 20078
    'SREBF2' in mth_trscp
    '''
'''
    ###########################################
    # METH_FL_NM = 'CCLE_RRBS_TSS1kb_20181022.txt'
    METH_FL_NM = 'CCLE_RRBS_tss_CpG_clusters_20181022.txt'
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    METHLT_FLD = 'data/methylation'
    ENSEMBL_SYMBOL_FILE = 'support/mappings/ensembl_symbol_human.json'
    FLD_PTH = 'results/ecGEM_simul_human1/init/no_tsks/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'
    methlt_dt_fld = METHLT_FLD
    GEN_MD_PTH='data/models/Human-GEM.yml'
    FLD_PTH='results/ecGEM_simul_human1/init/no_tsks/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    PRT_FL_NM='allprot_rc_substm.tsv' # protein!
    fld_pth = FLD_PTH
    fl_nm = PRT_FL_NM
    gen_gecko_pth = GEN_GECKO_PTH
    gen_md_pth = GEN_MD_PTH
    fld_pth = FLD_PTH
    meth_fl_nm = METH_FL_NM
    smp_info= ACHILLES_SMP_INFO
    ensembl_symbol_file=ENSEMBL_SYMBOL_FILE
    ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
    drw_rc = {r.id for r in ms.reactions if 'draw_' in r.id}
    ## Genes which corresponding reaction flux/protein usage is correlated with overall methylation level:
    if 'flx' in fl_nm:
        # get dataframe with reactions flux in different cell lines:
        flnm = os.path.join('/'.join(fld_pth.split('/')[:-1]), 'prot_pool_exchange1e+00_medonly_FBA_fluxes.tsv')
        df = pd.read_csv(flnm, sep='\t', index_col=0)
        df2 = df.drop({r for r in df.index if 'No' not in r}, axis=0)  # exclude reactions not associated with an enzyme/gene and 'draw' reactions
        # get same table but where each row is a reaction-gene combination
        # (one reaction can be associated with more than one gene and vice-versa:
        rc_gn = pd.DataFrame([(k, el) for k, v in {r: [g.id for g in ms.reactions.get_by_id(r).genes] for r in df2.index}.items() for el in v])
        rc_gn.rename(columns={0: 'reaction', 1: 'gene'}, inplace=True)
        rc_gn.index = rc_gn['reaction']
        rc_gn.drop('reaction', axis=1, inplace=True)
        rc_gn = rc_gn.join(df2, how='left')
    else:
        # get dataframe with protein usage for each reaction in different cell lines:
        prot_use = pd.read_csv(os.path.join(fld_pth, fl_nm), index_col=0, sep='\t')
        prot_use.index=prot_use[['Reaction','Enzyme']] # reactions not associated with an enzyme/gene and 'draw' reactions are already excluded
        prot_use.drop(['Enzyme', 'Reaction', 'Subsystem'], axis=1, inplace=True)
        prt_gene_d = {k: el.id for k, v in {mt.id: ms.reactions.get_by_id(rid).genes for rid in drw_rc for mt in ms.reactions.get_by_id(rid).metabolites if 'pool' not in mt.id}.items() for el in v} # dictionary {enzyme: gene}
        # get same table but where each row is a reaction-gene combination
        # (one reaction can be associated with more than one gene and vice-versa:
        prot_use.insert(loc=0, column='gene', value=[prt_gene_d[e] for r, e in prot_use.index])
        rc_gn=prot_use
    # get overall DNA methylation degree:
    mth_res_lst = TestMds.obtain_methylation(methlt_dt_fld=methlt_dt_fld, smp_info=smp_info)
    for el in mth_res_lst:
        if el.name == 'DNA_methylation_ratio': glb_mth = el
    glb_mth.name = 'Global DNA methylation'
    # correlate each gene-reaction combination with global DNA methylation:
    ncl_ord = list(set(rc_gn.iloc[:, 1:].columns).intersection(set(glb_mth.index)))  # 30 common
    rc_gnval = rc_gn.iloc[:, 1:]
    rc_gnval = rc_gnval[ncl_ord]
    glb_mth = glb_mth[ncl_ord]
    if enzflx_mth_corrtype == 'spearman':
        rc_gn_corr = rc_gnval.apply(lambda x: spearmanr(x, glb_mth), axis=1)  # NaN is when there is no variation
        corr_rc_gn = pd.DataFrame(map(lambda x: (x.correlation, x.pvalue), rc_gn_corr), index=rc_gnval.index, columns=['Coefficient', 'P-value'])
    elif enzflx_mth_corrtype == 'pearson':
        rc_gn_corr = rc_gnval.apply(lambda x: pearsonr(x, glb_mth), axis=1)  # NaN is when there is no variation
        corr_rc_gn = pd.DataFrame(map(lambda x: (x[0], x[1]), rc_gn_corr), index=rc_gnval.index, columns=['Coefficient', 'P-value'])
    corr_rc_gn['gene'] = rc_gn['gene']
    # get reaction-gene combinations with a significant correlation:
    corr_sig_rc_gn = corr_rc_gn[corr_rc_gn['P-value'] < 0.05]
    # from those above select those with direct or inverse correlation:
    if enzflx_bm_corrvalue == 'direct':
        corr_sig_rc_gn = corr_sig_rc_gn[corr_sig_rc_gn['Coefficient'] > 0]
    elif enzflx_bm_corrvalue == 'inverse':
        corr_sig_rc_gn = corr_sig_rc_gn[corr_sig_rc_gn['Coefficient'] < 0]
    # symbol- ensembl gene id in generic model:
    symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=ensembl_symbol_file)
    ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
    mdgenes = {g.id for g in ms.genes}
    # (symbol, ensembl) for all genes in the generic model:
    smbl_id_md = [(smb, id) for smb, id in symbol_ensembl if id in mdgenes]  # 2515 genes
    if len([id for s, id in smbl_id_md]) != len({id for s, id in smbl_id_md}):
        print("lookout! the same ensembl id has more than one symbol")
    elif len([s for s, id in smbl_id_md]) != len({s for s, id in smbl_id_md}):
        print("lookout! the same symbol has more than one ensembl id")
    ## Which of genes with corresponding reaction flux is correlated with overall methylation level are from fatty acid oxidation:
    m = load_yaml_model(gen_md_pth)  # generic traditional model
    sb_d = {r.id: r.subsystem[0] for r in m.reactions}
    if 'flx' in fl_nm:
        corr_sig_rc_gn_sb = [sb for r in corr_sig_rc_gn.index for id, sb in sb_d.items() if id in r]
    else:
        corr_sig_rc_gn_sb = [sb for r in corr_sig_rc_gn.index for id, sb in sb_d.items() if id in r[0]]
    corr_sig_rc_gn['subsystem'] = corr_sig_rc_gn_sb
    smbl_id_md_dct = {v:k for k,v in smbl_id_md}
    corr_sig_rc_gn['gene_name'] = [smbl_id_md_dct[id] for id in list(corr_sig_rc_gn['gene'])]
    fac = set(corr_sig_rc_gn.loc[corr_sig_rc_gn['subsystem'] == 'Fatty acid oxidation', 'gene_name'])
    # {'ACAA2', 'ACADS', 'ECHS1', 'ELOVL6', 'HADH', 'MECR', 'SLC27A2'}
    fac.intersection(mth_exp)
    # {'ACAA2', 'ECHS1', 'ELOVL6', 'SLC27A2'}
    fac.intersection(mth_mth)
    ## No!!!!Check which genes have individual methylation inversely correlated with flux in fatty acid oxidation:
    ## !!!!try force demethylation with an objective, test weights to choose value from plateau
    ## !!!!check if with demethylation there is less production of aKG (less tet activity) with increase in overall methylation
    ## !!!and which pathways follow the same trend with that diference.
    ## This frist!!!is DNMT more/less methylated and more/less expressed with increase in overal methylation?
    ## !!histone acetylation with NAD+/nicotinamide...

    ## Genes which methylation level significantly correlates with overall methylation level:
    # obtain experimental methylation degree of 1000bp upstream of TSS:
    raw_df = pd.read_csv(os.path.join(methlt_dt_fld, meth_fl_nm), sep='\t', index_col=0, na_values=['    NaN', '     NA']).iloc[:-1, :]
    meth_df = raw_df.iloc[:, 2:]
    mean_mth = meth_df.apply(lambda row: row.dropna().median(), axis=1)
    meth_df.apply(lambda row: row.fillna(value=mean_mth[row.name], inplace=True), axis=1)
    # get correlations:
    cl_ord = list(set(meth_df.columns).intersection(set(glb_mth.index)))  # 32 common cell lines to both datasets
    meth_df = meth_df[cl_ord]
    glb_mth = glb_mth.loc[cl_ord]
    if mth_bm_corrtype == 'spearman':
        res_mth = meth_df.apply(lambda x: spearmanr(x, glb_mth), axis=1)  # NaN is when there is no variation
        corr_meth_mth = pd.DataFrame(map(lambda x: (x.correlation, x.pvalue), res_mth), index=meth_df.index, columns=['Coefficient', 'P-value'])
    elif mth_bm_corrtype == 'pearson':
        res_mth = meth_df.apply(lambda x: pearsonr(x, glb_mth), axis=1)  # NaN is when there is no variation
        corr_meth_mth = pd.DataFrame(map(lambda x: (x[0], x[1]), res_mth), index=meth_df.index, columns=['Coefficient', 'P-value'])
    # get column with gene name:
    corr_meth_mth['gene_name'] = list(map(lambda x: x.split('_')[0], list(corr_meth_mth.index)))
    if all_prom_same_behavior:
        # select promoters of genes which all promoters have same signal (either are inversely or directly) correlated:
        total_prm = dict(collections.Counter(corr_meth_mth['gene_name']))
        to_keep = list()
        for g in total_prm:
            sbs = corr_meth_mth.loc[corr_meth_mth['gene_name'] == g, 'Coefficient']
            if np.sum(sbs < 0) == total_prm[g] or np.sum(sbs > 0) == total_prm[g]:
                to_keep.append(g)
        corr_same_meth = corr_meth_mth[corr_meth_mth['gene_name'].isin(to_keep)]
        # from above-mentioned promoters get those with a significant correlation:
        corr_same_sig_meth = corr_same_meth[corr_same_meth['P-value'] < 0.05]
    else:
        # from above-mentioned promoters get those with a significant correlation:
        corr_same_sig_meth = corr_meth_mth[corr_meth_mth['P-value'] < 0.05]
    # from those above select those with direct or inverse correlation:
    if mth_bm_corrvalue == 'direct':
        corr_same_sig_meth = corr_same_sig_meth[corr_same_sig_meth['Coefficient'] > 0]
    elif mth_bm_corrvalue == 'inverse':
        corr_same_sig_meth = corr_same_sig_meth[corr_same_sig_meth['Coefficient'] < 0]
    # get the significant correlated genes:
    mth_mth = set(corr_same_sig_meth['gene_name'])
    mth_mth_id = {id for s, id in smbl_id_md if s in mth_mth}
    # Note: many genes lost are because they are not part of model (not metabolic),
    # not because of conversion symbol to ensembl as len(mdgenes) is 2520 and len(2515),
    # so just 5 ensembl ids of the model didn't have corresponding symbol

    ## Genes which methylation is correlated with gene expression - this is always Pearson Correlation (provided by study):
    mth_exp_corr = pd.read_excel(os.path.join(methlt_dt_fld, 'meth_exp_corr.xlsx'), sheet_name='meth_exp_corr', index_col=0)
    mth_exp_corr['gene_name'] = list(map(lambda x: x.split('_')[0], list(mth_exp_corr.index)))
    # get p-value using n data points and correlation score:
    mth_exp_corr['pvalue'] = mth_exp_corr.apply(lambda x: PathAnalysis.r2p(r=x['Pearson_correlation'], N=x['n_datapoints']), axis=1)
    # select promoters of genes which all promoters have same signal (either are inversely or directly) correlated:
    total_prm = dict(collections.Counter(mth_exp_corr['gene_name']))
    to_keep = list()
    for g in total_prm:
        sbs = mth_exp_corr.loc[mth_exp_corr['gene_name'] == g, 'Pearson_correlation']
        if np.sum(sbs < 0) == total_prm[g] or np.sum(sbs > 0) == total_prm[g]:
            to_keep.append(g)
    corr_same = mth_exp_corr[mth_exp_corr['gene_name'].isin(to_keep)]
    # from above-mentioned promoters get those with a significant correlation:
    corr_same_sig = corr_same[corr_same['pvalue'] < 0.05]
    # from those above select those with inverse or direct correlation:
    if mth_exp_corrvalue == 'inverse':
        corr_same_sig = corr_same_sig.loc[corr_same_sig['Pearson_correlation'] < 0]
    elif mth_exp_corrvalue == 'direct':
        corr_same_sig = corr_same_sig.loc[corr_same_sig['Pearson_correlation'] > 0]
    # get the significant correlated genes:
    mth_exp = set(corr_same_sig['gene_name'])
    # convert to ensembl id present in the model:
    mth_exp_id = {id for s, id in smbl_id_md if s in mth_exp}
    # Note: many genes lost are because they are not part of model (not metabolic),
    # not because of conversion symbol to ensembl as len(mdgenes) is 2520 and len(2515),
    # so just 5 ensembl ids of the model didn't have corresponding symbol



    ## intersect above table with gene lists above:
    intr = mth_mth_id.intersection(mth_exp_id).intersection(set(corr_sig_rc_gn['gene']))
    # {'ENSG00000059377', 'ENSG00000138029'}
    mth_exp_id.intersection(set(corr_sig_rc_gn['gene']))
    # {'ENSG00000059377',
    #  'ENSG00000113790',
    #  'ENSG00000127884',
    #  'ENSG00000134333',
    #  'ENSG00000138029',
    #  'ENSG00000140284',
    #  'ENSG00000167315',
    #  'ENSG00000170522'}
    mth_mth_id.intersection(set(corr_sig_rc_gn['gene']))
    # {'ENSG00000059377', 'ENSG00000138029', 'ENSG00000173614'}
'''
if '__name__' == '__main__':
    ### Pathways which flux/protein usage correlates with overall methylation
    ## independently of growth rate:
    FLD_PTH='results/ecGEM_simul_human1/init/no_tsks/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    FLX_FL_NM='allflx_rc_substm.tsv'
    PRT_FL_NM='allprot_rc_substm.tsv'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    ENSEMBL_SYMBOL_FILE = 'support/mappings/ensembl_symbol_human.json'
    div_bm=True # obtain results independently of growth rate
    # pathways which fluxes in average significantly correlates (directly/inversely) with DNA methylation level/growth rate:
    pthw_corr_meth_bmindp(fld_pth=FLD_PTH, fl_nm=FLX_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO, gen_gecko_pth=GEN_GECKO_PTH, div_bm=div_bm)
    # pathways which protein usage in average significantly correlates (directly/inversely) with DNA methylation level/growth rate:
    pthw_corr_meth_bmindp(fld_pth=FLD_PTH, fl_nm=PRT_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO, fl_nm_bm=FLX_FL_NM, gen_gecko_pth=GEN_GECKO_PTH, div_bm=div_bm, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE)
    
    ### Pathways over-represented in genes which gene expression correlates with overall DNA methylation:
    ## independently of growth rate:
    METHLT_FLD='data/methylation'
    TRSCP_PTH='data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.txt'
    ENVCOND='data/constraints/1218595databases1_corrected_further_cal.xls'
    GEN_GECKO_PTH='support/models/prodDNAtot/ecModel_batch.mat'
    ENSEMBL_SYMBOL_FILE='support/mappings/ensembl_symbol_human.json'
    ACHILLES_SMP_INFO='data/sample_info.csv'
    GEN_MD_PTH_SBML='support/models/prodDNAtot.xml'  # traditional generic model with DNA meth and demethylation reactions
    GEN_MD_PTH='data/models/Human-GEM.yml'
    FLD_PTH='results/ecGEM_simul_human1/init/no_tsks/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    DIV_GROWTH=True
    mth_trscp_corrvalue='direct'
    experimental_corr_test(methlt_dt_fld=METHLT_FLD, trscp_pth=TRSCP_PTH, mth_trscp_corrvalue=mth_trscp_corrvalue, envcond=ENVCOND,
                           gen_gecko_pth=GEN_GECKO_PTH, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, smp_info=ACHILLES_SMP_INFO,
                           gen_md_pth_sbml=GEN_MD_PTH_SBML, gen_md_pth=GEN_MD_PTH, fld_pth=FLD_PTH, div_growth=DIV_GROWTH)
    mth_trscp_corrvalue='inverse'
    experimental_corr_test(methlt_dt_fld=METHLT_FLD, trscp_pth=TRSCP_PTH, mth_trscp_corrvalue=mth_trscp_corrvalue, envcond=ENVCOND,
                           gen_gecko_pth=GEN_GECKO_PTH, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, smp_info=ACHILLES_SMP_INFO,
                           gen_md_pth_sbml=GEN_MD_PTH_SBML, gen_md_pth=GEN_MD_PTH, fld_pth=FLD_PTH, div_growth=DIV_GROWTH)

    ### Pathways over-represented in genes which methylation correlates directly with overall methylation and inversely with gene expression:
    ## gene promoter:
    METH_FL_NM='CCLE_RRBS_TSS1kb_20181022.txt'
    METHLT_FLD='data/methylation'
    ACHILLES_SMP_INFO='data/sample_info.csv'
    ENSEMBL_SYMBOL_FILE='support/mappings/ensembl_symbol_human.json'
    GEN_GECKO_PTH='support/models/prodDNAtot/ecModel_batch.mat'
    GEN_MD_PTH='data/models/Human-GEM.yml'
    GEN_MD_PTH_SBML='support/models/prodDNAtot.xml'  # traditional generic model with DNA meth and demethylation reactions
    FLD_PTH='results/ecGEM_simul_human1/init/no_tsks/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    FL_NM='allflx_rc_substm.tsv'
    ENVCOND='data/constraints/1218595databases1_corrected_further_cal.xls'
    all_prom_same_behavior=False  # whether to use only genes where all promoters are either directly or inversely correlated
                                    # or instead consider any gene with at least one promoter with the behavior we want (direct/inversely correlated)
    '''
    # not correcting for biomass influence:
    mth_mth_corrvalue='direct'  # whether to use direct or inverse correlation between overall methylation and gene methylation
    div_growth=False  # whether to divide by the experimental growth rate
    mth_exp_corrvalue='inverse'  # whether to use direct or inverse correlation between gene methylation and gene expression
    experimental_corr(meth_fl_nm=METH_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, gen_md_pth=GEN_MD_PTH,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, all_prom_same_behavior=all_prom_same_behavior, mth_mth_corrvalue=mth_mth_corrvalue, fld_pth=FLD_PTH,
                      envcond=ENVCOND, div_growth=div_growth, mth_exp_corrvalue=mth_exp_corrvalue)
    mth_exp_corrvalue = 'direct'  # whether to use direct or inverse correlation between gene methylation and gene expression
    experimental_corr(meth_fl_nm=METH_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO,
                      ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, gen_md_pth=GEN_MD_PTH,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, all_prom_same_behavior=all_prom_same_behavior,
                      mth_mth_corrvalue=mth_mth_corrvalue, fld_pth=FLD_PTH,
                      envcond=ENVCOND, div_growth=div_growth, mth_exp_corrvalue=mth_exp_corrvalue)
    '''
    # correcting for biomass influence:
    mth_mth_corrvalue = 'direct'  # whether to use direct or inverse correlation between overall methylation and gene methylation
    div_growth = True  # whether to divide by the experimental growth rate
    mth_exp_corrvalue = 'inverse'  # whether to use direct or inverse correlation between gene methylation and gene expression
    experimental_corr(meth_fl_nm=METH_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO,
                      ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, gen_md_pth=GEN_MD_PTH,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, all_prom_same_behavior=all_prom_same_behavior,
                      mth_mth_corrvalue=mth_mth_corrvalue, fld_pth=FLD_PTH,
                      envcond=ENVCOND, div_growth=div_growth, mth_exp_corrvalue=mth_exp_corrvalue)
    mth_exp_corrvalue = 'direct'  # whether to use direct or inverse correlation between gene methylation and gene expression
    experimental_corr(meth_fl_nm=METH_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO,
                      ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, gen_md_pth=GEN_MD_PTH,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, all_prom_same_behavior=all_prom_same_behavior,
                      mth_mth_corrvalue=mth_mth_corrvalue, fld_pth=FLD_PTH,
                      envcond=ENVCOND, div_growth=div_growth, mth_exp_corrvalue=mth_exp_corrvalue)
    ## Pathways over-represented in genes which methylation correlates directly with overall methylation and inversely with gene expression:
    ## TSS cluster:
    METH_FL_NM = 'CCLE_RRBS_tss_CpG_clusters_20181022.txt'  # Here is the difference from the above
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'
    ENSEMBL_SYMBOL_FILE = 'support/mappings/ensembl_symbol_human.json'
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    GEN_MD_PTH = 'data/models/Human-GEM.yml'
    GEN_MD_PTH_SBML = 'support/models/prodDNAtot.xml'  # traditional generic model with DNA meth and demethylation reactions
    FLD_PTH = 'results/ecGEM_simul_human1/init/no_tsks/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    FL_NM = 'allflx_rc_substm.tsv'
    ENVCOND = 'data/constraints/1218595databases1_corrected_further_cal.xls'
    all_prom_same_behavior = False  # whether to use only genes where all promoters are either directly or inversely correlated
    # or instead consider any gene with at least one promoter with the behavior we want (direct/inversely correlated)
    # correcting for biomass influence:
    mth_mth_corrvalue = 'direct'  # whether to use direct or inverse correlation between overall methylation and gene methylation
    div_growth = True  # whether to divide by the experimental growth rate
    mth_exp_corrvalue = 'inverse'  # whether to use direct or inverse correlation between gene methylation and gene expression
    experimental_corr(meth_fl_nm=METH_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO,
                      ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, gen_md_pth=GEN_MD_PTH,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, all_prom_same_behavior=all_prom_same_behavior,
                      mth_mth_corrvalue=mth_mth_corrvalue, fld_pth=FLD_PTH,
                      envcond=ENVCOND, div_growth=div_growth, mth_exp_corrvalue=mth_exp_corrvalue)
    mth_exp_corrvalue = 'direct'  # whether to use direct or inverse correlation between gene methylation and gene expression
    experimental_corr(meth_fl_nm=METH_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO,
                      ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, gen_md_pth=GEN_MD_PTH,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, all_prom_same_behavior=all_prom_same_behavior,
                      mth_mth_corrvalue=mth_mth_corrvalue, fld_pth=FLD_PTH,
                      envcond=ENVCOND, div_growth=div_growth, mth_exp_corrvalue=mth_exp_corrvalue)

    ## Genes:
    # enzyme usage
    # - enzyme vd growth rate: spearman, inverse
    # - methylation vs gene expression: pearson (fixed), inverse
    # - methylation vs growth rate: spearman, direct
    GEN_GECKO_PTH='support/models/prodDNAtot/ecModel_batch.mat'
    METH_FL_NM='CCLE_RRBS_TSS1kb_20181022.txt'
    METHLT_FLD = 'data/methylation'
    FLD_PTH='results/ecGEM_simul_human1/init/no_tsks/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    ENSEMBL_SYMBOL_FILE='support/mappings/ensembl_symbol_human.json'
    enzflx='enzyme' # either 'enzyme' usage or reaction 'flux'
    enzflx_bm_corrtype='spearman' # method of correlation between enzyme usage/flux and cell growth rate
    enzflx_bm_corrvalue='inverse' # either 'inverse' or 'direct' correlation between enzyme usage/flux and cell growth rate
    mth_exp_corrvalue='inverse' # either 'inverse' or 'direct' correlation between gene methylation and gene expression
    mth_bm_corrtype='spearman' # method of correlation between gene methylation and cell growth rate
    mth_bm_corrvalue='direct' # either 'inverse' or 'direct' correlation between gene methylation and cell growth rate
    common_correlated_genes(gen_gecko_pth=GEN_GECKO_PTH, enzflx=enzflx, fld_pth=FLD_PTH,
                            enzflx_bm_corrtype=enzflx_bm_corrtype, enzflx_bm_corrvalue=enzflx_bm_corrvalue, methlt_dt_fld=METHLT_FLD,
                            ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, mth_exp_corrvalue=mth_exp_corrvalue, mth_bm_corrtype=mth_bm_corrtype,
                            mth_bm_corrvalue=mth_bm_corrvalue, meth_fl_nm=METH_FL_NM)

    # reaction flux
    # - flux vd growth rate: spearman, inverse
    # - methylation vs gene expression: pearson (fixed), inverse
    # - methylation vs growth rate: spearman, direct
    enzflx = 'flux'  # either 'enzyme' usage or reaction 'flux'
    enzflx_bm_corrtype = 'spearman'  # method of correlation between enzyme usage/flux and cell growth rate
    enzflx_bm_corrvalue = 'inverse'  # either 'inverse' or 'direct' correlation between enzyme usage/flux and cell growth rate
    mth_exp_corrvalue = 'inverse'  # either 'inverse' or 'direct' correlation between gene methylation and gene expression
    mth_bm_corrtype = 'spearman'  # method of correlation between gene methylation and cell growth rate
    mth_bm_corrvalue = 'direct'  # either 'inverse' or 'direct' correlation between gene methylation and cell growth rate
    common_correlated_genes(gen_gecko_pth=GEN_GECKO_PTH, enzflx=enzflx, fld_pth=FLD_PTH,
                            enzflx_bm_corrtype=enzflx_bm_corrtype, enzflx_bm_corrvalue=enzflx_bm_corrvalue,
                            methlt_dt_fld=METHLT_FLD, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE,
                            mth_exp_corrvalue=mth_exp_corrvalue, mth_bm_corrtype=mth_bm_corrtype,
                            mth_bm_corrvalue=mth_bm_corrvalue, meth_fl_nm=METH_FL_NM)


####################################
def common_correlated_genes(gen_gecko_pth, enzflx, fld_pth, enzflx_bm_corrtype, enzflx_bm_corrvalue,
                            methlt_dt_fld, ensembl_symbol_file, mth_exp_corrvalue, mth_bm_corrtype,
                            mth_bm_corrvalue, meth_fl_nm):
    '''
    - gets genes which enzyme usage/ associated reaction flux is correlated (either directly/inversely) with biomass
    - gets genes which methylation is correlated (either directly/inversely) with their transcription
    - gets genes which methylation is correlated (either directly/inversely) with cell growth rate
    - then intersects the three aforementioned groups of genes.
    :param gen_gecko_pth: path to generic gecko model
    :param enzflx: either 'enzyme' or 'flux' if wanting to use the enzyme usage (mg/gDW) or the reaction fluxes (mmol/(gDW.h))
    :param fld_pth: path to folder containing fluxes (mmol/gDW.h) or protein usage (mg/gDW)
                    for reactions in cell line-specific gecko models.
                    For the table with fluxes, the reactions excluded from the table are: 'draw_prot', 'prot_pool_exchange', and leg reactions corresponding to 'arm_' reactions (cause 'arm_' reactions flux is already representative of all the associated legs).
                    For the table with protein usage, the reactions excluded from the table are: 'draw_prot', 'prot_pool_exchange', 'arm_' reactions (cause those do not have an associated enzyme)

    :param enzflx_bm_corrtype: either 'spearman' or 'pearson'. type of correlation method used to correlate enzyme usage/reaction flux with cell growth rate
    :param enzflx_bm_corrvalue: either 'inverse' or 'direct'. to indicate if we want direct or indirect correlation between enzyme usage/reaction flux and cell growth rate
    :param methlt_dt_fld: path to folder where file with experimental DNA methylation fluxes is.
    :param ensembl_symbol_file: path to .json file where to retrieve the mappings or where to save the mappings (between ensembl gene id and gene symbol) if the file does not exist yet.
    :param mth_exp_corrvalue: either 'inverse' or 'direct'. to indicate if we want direct or indirect correlation between gene methylation and gene expression
    :param mth_bm_corrtype: either 'spearman' or 'pearson'. type of correlation method used to correlate gene methylation and cell growth rate.
    :param mth_bm_corrvalue: either 'inverse' or 'direct'. to indicate if we want direct or indirect correlation between gene methylation and cell growth rate.
    :param meth_fl_nm: path to file with experimental methylation we want to use. can be from TSS or other genomic regions.
    '''
    # Note: when spearman method gives at least some significant correlation coefficients while the Pearson method does Not,
    # that means that there is no linear correlation but exists a monotonic relationship.

    # get dictionary {enzyme: gene}:
    ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
    drw_rc = {r.id for r in ms.reactions if 'draw_' in r.id}
    prt_gene_d = {k: [el.id] for k,v in {mt.id: ms.reactions.get_by_id(rid).genes for rid in drw_rc for mt in ms.reactions.get_by_id(rid).metabolites if 'pool' not in mt.id}.items() for el in v}

    # get cell growth rates:
    flnm = os.path.join('/'.join(fld_pth.split('/')[:-1]), 'prot_pool_exchange1e+00_medonly_FBA_fluxes.tsv')
    df = pd.read_csv(flnm, sep='\t', index_col=0)
    biomass = df.loc['adaptbiomass']

    if enzflx == 'enzyme':
        ## Genes which corresponding Enzyme use is correlated with biomass:
        # get table with gene, p-value and correlation coefficient of enzymes that correlate significantly (direct/reverse) with biomass:
        enz_crr = pd.read_csv(os.path.join(fld_pth, enzflx_bm_corrtype, 'corr_biomass_prot_rc.tsv'), sep='\t', index_col=0)
        joined = enz_crr.join(pd.DataFrame(prt_gene_d).T)
        joined.rename(columns={0:'Gene'}, inplace=True)
        # from those enzymes get those correlating:
        if enzflx_bm_corrvalue == 'inverse':
            dff = joined.loc[joined['Coefficient'] < 0, :]
        elif enzflx_bm_corrvalue == 'direct':
            dff = joined.loc[joined['Coefficient'] > 0, :]
        ensembl_symbol = PathAnalysis.get_ensembl_mappings(fl_json=ensembl_symbol_file)
        enzflx_bm = set(map(lambda x: ensembl_symbol[x], dff.loc[:, 'Gene']))
    elif enzflx == 'flux':
        ## Genes which corresponding reaction flux is correlated with biomass:
        # get dataframe with reactions flux in different cell lines:
        df2 = df.drop({r for r in df.index if 'No' not in r}, axis=0) # exclude reactions not associated with an enzyme/gene and 'draw' reactions
        # get same table but where each row is a reaction-gene combination
        # (one reaction can be associated with more than one gene and vice-versa:
        rc_gn = pd.DataFrame([(k, el) for k, v in {r: [g.id for g in ms.reactions.get_by_id(r).genes] for r in df2.index}.items() for el in v])
        rc_gn.rename(columns={0: 'reaction', 1: 'gene'}, inplace=True)
        rc_gn.index = rc_gn['reaction']
        rc_gn.drop('reaction', axis=1, inplace=True)
        rc_gn = rc_gn.join(df2, how='left')
        # correlate each gene-reaction combination with biomass:
        if enzflx_bm_corrtype == 'spearman':
            rc_gn_corr = rc_gn.iloc[:, 1:].apply(lambda x: spearmanr(x, biomass), axis=1) # NaN is when there is no variation
        elif enzflx_bm_corrtype == 'pearson':
            rc_gn_corr = rc_gn.iloc[:, 1:].apply(lambda x: pearsonr(x, biomass), axis=1)  # NaN is when there is no variation
        corr_rc_gn = pd.DataFrame(map(lambda x: (x.correlation, x.pvalue), rc_gn_corr), index=rc_gn.index, columns=['Coefficient', 'P-value'])
        corr_rc_gn = pd.concat([rc_gn.loc[:, 'gene'].to_frame(), corr_rc_gn], axis=1)
        # get reaction-gene combinations with a significant correlation:
        corr_sig_rc_gn = corr_rc_gn[corr_rc_gn['P-value'] < 0.05]
        # from those above select those with direct or inverse correlation:
        if enzflx_bm_corrvalue == 'direct':
            corr_sig_rc_gn = corr_sig_rc_gn[corr_sig_rc_gn['Coefficient'] > 0]
        elif enzflx_bm_corrvalue == 'inverse':
            corr_sig_rc_gn = corr_sig_rc_gn[corr_sig_rc_gn['Coefficient'] < 0]
        # get the significant correlated genes:
        enzflx_bm = set(corr_sig_rc_gn['gene'])

    ## Genes which methylation is correlated with gene expression - this is always Pearson Correlation (provided by study):
    # get genes that are significantly inversely correlated with gene expression:
    mth_exp_corr = pd.read_excel(os.path.join(methlt_dt_fld, 'meth_exp_corr.xlsx'), sheet_name='meth_exp_corr', index_col=0)
    mth_exp_corr['gene_name'] = list(map(lambda x: x.split('_')[0], list(mth_exp_corr.index)))
    # get p-value using n data points and correlation score:
    mth_exp_corr['pvalue'] = mth_exp_corr.apply(lambda x: PathAnalysis.r2p(r=x['Pearson_correlation'], N=x['n_datapoints']), axis=1)
    # select promoters of genes which all promoters have same signal (either are inversely or directly) correlated:
    total_prm = dict(collections.Counter(mth_exp_corr['gene_name']))
    to_keep = list()
    for g in total_prm:
        sbs = mth_exp_corr.loc[mth_exp_corr['gene_name']==g, 'Pearson_correlation']
        if np.sum(sbs < 0) == total_prm[g] or np.sum(sbs > 0) == total_prm[g]:
            to_keep.append(g)
    corr_same = mth_exp_corr[mth_exp_corr['gene_name'].isin(to_keep)]
    # from above-mentioned promoters get those with a significant correlation:
    corr_same_sig = corr_same[corr_same['pvalue'] < 0.05]
    # from those above select those with inverse or direct correlation:
    if mth_exp_corrvalue == 'inverse':
        corr_same_sig = corr_same_sig[corr_same_sig['Pearson_correlation'] < 0]
    elif mth_exp_corrvalue == 'direct':
        corr_same_sig = corr_same_sig[corr_same_sig['Pearson_correlation'] > 0]
    # get the significant correlated genes:
    mth_exp = set(corr_same_sig['gene_name'])

    ## Genes which methylation level significantly correlates with biomass:
    # obtain experimental methylation degree of 1000bp upstream of TSS:
    raw_df = pd.read_csv(os.path.join(methlt_dt_fld, meth_fl_nm), sep='\t', index_col=0, na_values=['    NaN', '     NA']).iloc[:-1, :]
    meth_df = raw_df.iloc[:, 2:]
    mean_mth = meth_df.apply(lambda row: row.dropna().median(), axis=1).to_dict()
    meth_df.apply(lambda row: row.fillna(value=mean_mth[row.name], inplace=True), axis=1)
    # get spearman correlations:
    cl_ord = set(meth_df.columns).intersection(set(biomass.index)) # 30 common cell lines to both datasets
    meth_df = meth_df[cl_ord]
    biomass = biomass.loc[cl_ord]
    if mth_bm_corrtype == 'spearman':
        res_bm = meth_df.apply(lambda x: spearmanr(x, biomass), axis=1) # NaN is when there is no variation
    elif mth_bm_corrtype == 'pearson':
        res_bm = meth_df.apply(lambda x: pearsonr(x, biomass), axis=1)  # NaN is when there is no variation
    corr_meth_bm = pd.DataFrame(map(lambda x: (x.correlation, x.pvalue), res_bm), index=meth_df.index, columns=['Coefficient', 'P-value'])
    # get column with gene name:
    corr_meth_bm['gene_name'] = list(map(lambda x: x.split('_')[0], list(corr_meth_bm.index)))
    # select promoters of genes which all promoters have same signal (either are inversely or directly) correlated:
    total_prm = dict(collections.Counter(corr_meth_bm['gene_name']))
    to_keep = list()
    for g in total_prm:
        sbs = corr_meth_bm.loc[corr_meth_bm['gene_name']==g, 'Coefficient']
        if np.sum(sbs < 0) == total_prm[g] or np.sum(sbs > 0) == total_prm[g]:
            to_keep.append(g)
    corr_same_meth_bm = corr_meth_bm[corr_meth_bm['gene_name'].isin(to_keep)]
    # from above-mentioned promoters get those with a significant correlation:
    corr_same_sig_meth_bm = corr_same_meth_bm[corr_same_meth_bm['P-value'] < 0.05]
    # from those above select those with direct or inverse correlation:
    if mth_bm_corrvalue == 'direct':
        corr_same_sig_meth_bm = corr_same_sig_meth_bm[corr_same_sig_meth_bm['Coefficient'] > 0]
    elif mth_bm_corrvalue == 'inverse':
        corr_same_sig_meth_bm = corr_same_sig_meth_bm[corr_same_sig_meth_bm['Coefficient'] < 0]
    # get the significant correlated genes:
    mth_bm = set(corr_same_sig_meth_bm['gene_name'])

    ## intersect all groups of genes and save:
    intsec = enzflx_bm.intersection(mth_exp).intersection(mth_bm)
    ress = pd.DataFrame(intsec)
    a = os.path.join(fld_pth, 'corr_genes')
    b = os.path.join(a, f'{enzflx[:3]}-bm_{enzflx_bm_corrtype[:4]}_{enzflx_bm_corrvalue[:3]}')
    c = os.path.join(b, f'mth-exp_{mth_exp_corrvalue[:3]}')
    d = os.path.join(c, f'mth-bm_{mth_bm_corrtype[:4]}_{mth_bm_corrvalue}.tsv')
    for sb in [a, b, c]:
        if not os.path.exists(sb):
            os.makedirs(sb)
    ress.to_csv(d, sep='\t')
    print(intsec)