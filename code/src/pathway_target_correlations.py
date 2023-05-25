import os
import pandas as pd
from test_mds import TestMds
from scipy.stats import pearsonr, spearmanr
import copy
from pathway_analysis import PathAnalysis
import numpy as np
from cobra.io import load_matlab_model, load_yaml_model, read_sbml_model

class DataCorr:
    '''class with methods to analyze correlations between different data types (flux, protein usage, gene expression, gene methylation)'''
    @staticmethod
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
        # get only the simulated fluxes of all reactions (exclude real methylation levels) except biomass:
        if 'flx' in fl_nm:
            ndf2 = df.drop('adaptbiomass', axis=1)
        else:
            ndf2 = copy.deepcopy(df)
        ndf2 = ndf2.loc[:, set(sim_vars) - {'adaptbiomass'}]
        # get global DNA methylation and biomass:
        gb = df['Global DNA methylation'].drop('Subsystem', axis=0)  # exclude 'Subsystem'
        if 'flx' in fl_nm:
            biomass = df['adaptbiomass'].drop('Subsystem')  # exclude 'Subsystem'
        else:
            # df_w_bm = pd.read_csv(os.path.join(fld_pth, fl_nm_bm), index_col=0, sep='\t')
            df_w_bm = pd.read_csv(os.path.join(fld_pth, kwargs['fl_nm_bm']), index_col=0, sep='\t')
            biomass = df_w_bm.loc['adaptbiomass'].drop('Subsystem')
            ndf2.columns = pd.MultiIndex.from_arrays(ndf2.iloc[0:2].values)
            ndf2.drop(['Enzyme', 'Reaction'], axis=0, inplace=True)
            gb.drop(['Enzyme', 'Reaction'], axis=0, inplace=True)
        # corrected by (divide by) biomass flux or not:
        if div_bm:
            rc_by_bm = ndf2.drop('Subsystem', axis=0).divide(biomass, axis=0)  # reaction flux/protein usage divided by biomass
        else:
            rc_by_bm = ndf2.drop('Subsystem', axis=0)
        # aggregate by subsystem:
        sbsm_by_bm_old = rc_by_bm.append(pd.DataFrame(ndf2.loc['Subsystem', :]).T)
        sbsm_by_bm = sbsm_by_bm_old.T.groupby('Subsystem').mean().T  # each average substm flux/protein usage divided by biomass
        # calculate coefficient and P-value:
        # with spearman correlation:
        res_spm = sbsm_by_bm.apply(lambda x: spearmanr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
        res_spm.index = ['Coefficient', 'P-value']
        # with pearson correlation:
        res_prs = sbsm_by_bm.apply(lambda x: pearsonr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
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
        spm_pos = sig_sys_spm.loc[sig_sys_spm['Coefficient'] > 0]
        spm_neg = sig_sys_spm.loc[sig_sys_spm['Coefficient'] < 0]
        prs_pos = sig_sys_prs.loc[sig_sys_prs['Coefficient'] > 0]
        prs_neg = sig_sys_prs.loc[sig_sys_prs['Coefficient'] < 0]
        common_spm_prs_pos = set(spm_pos.index).intersection(prs_pos.index)
        common_spm_prs_neg = set(spm_neg.index).intersection(prs_neg.index)
        sig_sys_common_pos = sig_sys_prs.loc[
            common_spm_prs_pos]  # common to spearman and pearson, but coefficient from pearson
        sig_sys_common_neg = sig_sys_prs.loc[
            common_spm_prs_neg]  # common to spearman and pearson, but coefficient from pearson
        sig_sys_common_pos.sort_values(by='Coefficient', ascending=False, inplace=True)
        sig_sys_common_neg.sort_values(by='Coefficient', ascending=False, inplace=True)
        # save results:
        if 'flx' in fl_nm:
            el = 'flx'
        elif 'prot' in fl_nm:
            el = 'prot'
        if div_bm:
            bmv = '_divbm_'
        else:
            bmv = ''
        sig_sys_spm.to_csv(os.path.join(fld_pth, 'spearman', f'corr_meth_{el}_substm{bmv}.tsv'), sep='\t')
        sig_sys_prs.to_csv(os.path.join(fld_pth, 'pearson', f'corr_meth_{el}_substm{bmv}.tsv'), sep='\t')
        sig_sys_common_pos.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_{el}_substm_pos{bmv}.tsv'), sep='\t')
        sig_sys_common_neg.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_{el}_substm_neg{bmv}.tsv'), sep='\t')
        if 'flx' in fl_nm:
            ## identify reactions of each subsystem which flux correlates with overall DNA methylation,
            ## independently of cell growth rates:
            # with spearman correlation:
            spm = sbsm_by_bm_old.iloc[:-1, :].apply(lambda x: spearmanr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
            spm.index = ['Coefficient', 'P-value']
            spm.loc['Subsystem'] = sbsm_by_bm_old.loc['Subsystem']
            # with pearson correlation:
            prs = sbsm_by_bm_old.iloc[:-1, :].apply(lambda x: pearsonr(x, pd.DataFrame(gb)), axis=0)  # NaN is when there is no variation
            prs.index = ['Coefficient', 'P-value']
            prs.loc['Subsystem'] = sbsm_by_bm_old.loc['Subsystem']
            # select metabolic reactions with p-value < 0.05:
            sig_spm = spm.loc[:, spm.loc['P-value'] < 0.05].T
            sig_prs = prs.loc[:, prs.loc['P-value'] < 0.05].T

            # get rid of weird subsystems:
            def remove_subst_cl(f):
                mett = ['Isolated', 'Miscellaneous', 'Artificial reactions']
                for x in mett:
                    if sum(f['Subsystem'] == x) != 0:
                        f = f.drop(f.loc[f['Subsystem'] == x].index, axis=0)
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
                rgs = ' + '.join(
                    [mtid_nm[rg.id] + f'[{rg.compartment}]' for rg in ms.reactions.get_by_id(x).reactants if
                     'prot_' not in rg.id])
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
            sig_common_pos['Gene_id'] = [prt_gene_d[el] for el in sig_common_pos.index]
            sig_common_neg['Gene_id'] = [prt_gene_d[el] for el in sig_common_neg.index]
            # symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=ensembl_symbol_file)
            symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=kwargs['ensembl_symbol_file'])
            id_smb_pos = {id: smb for smb, id in symbol_ensembl if id in list(sig_common_pos['Gene_id'])}
            sig_common_pos['Symbol'] = [id_smb_pos[el] for el in sig_common_pos['Gene_id']]
            id_smb_neg = {id: smb for smb, id in symbol_ensembl if id in list(sig_common_neg['Gene_id'])}
            sig_common_neg['Symbol'] = [id_smb_neg[el] for el in sig_common_neg['Gene_id']]
            # save results:
            sig_common_pos.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_prot_usage_pos{bmv}.tsv'), sep='\t')
            sig_common_neg.to_csv(os.path.join(fld_pth, 'common_butcoefpearson', f'corr_meth_prot_usage_neg{bmv}.tsv'), sep='\t')

    @staticmethod
    def experimental_corr(methlt_dt_fld, trscp_pth, mth_trscp_corrvalue, envcond, gen_gecko_pth, ensembl_symbol_file,
                          smp_info, gen_md_pth_sbml, gen_md_pth, fld_pth, div_growth):
        '''
        - identify metabolic subsystems over-represented in the list of reactions associated with:
         * genes which expression level correlate (positively/negatively) with overall DNA methylation level
         (if required independently of cell growth rate) across different cell lines
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
        # get overall DNA methylation degree:
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
        if div_growth:
            bm = '_div_bm'
        else:
            bm = ''
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

    @staticmethod
    def mth_bm_trgt(ensembl_symbol_file, gen_gecko_pth, trscp_pth, meth_fl, fld_pth, fl_nm, **kwargs):
        '''
        - genes whose promoter methylation significantly correlated (either positively/negatively) with its transcription
        and with the cell growth rate across the different cell lines
        is intersected with either the group of genes associated with reactions whose flux,
        or the genes associated with enzymes whose protein usage,
        significantly correlated (either positively/negatively) with the cell growth rate.
        :param ensembl_symbol_file: path to .json file where to retrieve or where to save (if the file was not created before) the mappings between ensembl gene id and gene symbol.
        :param gen_gecko_pth: path to generic gecko model with DNA (de)/methylation reactions, but NO subsystem information
        :param trscp_pth: path to file with gene expression data
        :param meth_fl: path to file with experimental methylation we want to use.
        :param fld_pth: path to folder containing fluxes (mmol/gDW.h) or protein usage (mg/gDW)
                        for reactions in cell line-specific gecko models.
                        For the table with fluxes, the reactions excluded from the table are: 'draw_prot', 'prot_pool_exchange', and leg reactions corresponding to 'arm_' reactions (cause 'arm_' reactions flux is already representative of all the associated legs).
                        For the table with protein usage, the reactions excluded from the table are: 'draw_prot', 'prot_pool_exchange', 'arm_' reactions (cause those do not have an associated enzyme)
        :param fl_nm: name of file. has to contain either 'flx' or 'prot' depending on whether relates to fluxes or protein usage
        :param kwargs - fl_nm_bm: path to dataframe with biomass flux. use when analyzing protein usage.
        '''

        ### genes which expression correlates with its methylation:
        # load gene promoter methylation:
        raw_meth_df = pd.read_csv(meth_fl, sep='\t', index_col=0, na_values=['    NaN', '     NA']).iloc[:-1, :]
        meth_df = raw_meth_df.iloc[:, 2:]
        mean_mth = meth_df.apply(lambda row: row.dropna().median(), axis=1).to_dict()
        meth_df.apply(lambda row: row.fillna(value=mean_mth[row.name], inplace=True), axis=1)
        meth_df = meth_df.drop('ITFG2-AS1_12_2922090_2923090', axis=0)  # this gene encoding a ncRNA has more than one possible ensembl gene id
        # load gene expression matrix:
        trscp = pd.read_csv(trscp_pth, sep='\t')
        trscp.index = trscp['gene_id'].apply(lambda x: x.split('.')[0])
        trscp.drop(['gene_id', 'transcript_ids'], axis=1, inplace=True)
        # get gene expression symbols:
        symbol_ensembl = PathAnalysis.get_ensembl_mappings(fl_json=ensembl_symbol_file)
        ens_symb_trscp = [(gid, s) for gid in trscp.index for s, ens in symbol_ensembl if gid == ens]
        ens_sym_exp = dict()
        for ens, sym in ens_symb_trscp:
            if ens in ens_sym_exp:
                ens_sym_exp[ens].append(sym)
            else:
                ens_sym_exp[ens] = [sym]
        trscp['gene_symbol'] = [ens_sym_exp[tr] if tr in ens_sym_exp else [''] for tr in trscp.index]
        # calculate coefficient and P-value:
        cl_ord = list(set(meth_df.columns).intersection(set(trscp.columns)))  # 836 common cell lines to both datasets
        meth_df = meth_df[cl_ord]
        exp_df = trscp[cl_ord + ['gene_symbol']]

        def get_row_spearman(x):
            print(x.name)
            bool_sl = [x.name.split('_')[0] in el for el in exp_df['gene_symbol']]
            if sum(bool_sl):
                return spearmanr(x, exp_df.iloc[np.array(bool_sl), :-1].T)
            else:
                return (float('nan'), float('nan'))

        def get_row_pearson(x):
            print(x.name)
            bool_sl = [x.name.split('_')[0] in el for el in exp_df['gene_symbol']]
            return pearsonr(x, exp_df.iloc[np.array(bool_sl), :-1].T)

        # with spearman correlation:
        meth_exp_spm_t = [get_row_spearman(row) for _, row in meth_df.iterrows()]  # NaN is when there is no variation
        meth_exp_spm = pd.DataFrame(np.array(meth_exp_spm_t).T)
        meth_exp_spm.columns = meth_df.index
        meth_exp_spm.index = ['Coefficient', 'P-value']
        # with pearson correlation:
        meth_exp_prs_t = [get_row_pearson(row) for _, row in meth_df.iterrows()]  # NaN is when there is no variation
        meth_exp_prs = pd.DataFrame(np.array(meth_exp_prs_t).T)
        meth_exp_prs.columns = meth_df.index
        meth_exp_prs.index = ['Coefficient', 'P-value']
        meth_exp_prs = meth_exp_prs.astype(float)
        # select with p-value < 0.05:
        meth_exp_spm = meth_exp_spm.loc[:, meth_exp_spm.loc['P-value'] < 0.05].T
        meth_exp_prs = meth_exp_prs.loc[:, meth_exp_prs.loc['P-value'] < 0.05].T
        # add column with gene name:
        meth_exp_spm['Gene'] = meth_exp_spm.index.map(lambda x: x.split('_')[0])
        meth_exp_prs['Gene'] = meth_exp_prs.index.map(lambda x: x.split('_')[0])
        # intersect positively/negatively correlated of spearman and pearson methods:
        meth_exp_spm_pos = meth_exp_spm.loc[meth_exp_spm['Coefficient'] > 0]
        meth_exp_spm_neg = meth_exp_spm.loc[meth_exp_spm['Coefficient'] < 0]
        meth_exp_prs_pos = meth_exp_prs.loc[meth_exp_prs['Coefficient'] > 0]
        meth_exp_prs_neg = meth_exp_prs.loc[meth_exp_prs['Coefficient'] < 0]
        common_meth_exp_spm_prs_pos = set(meth_exp_spm_pos['Gene']).intersection(meth_exp_prs_pos['Gene'])
        common_meth_exp_spm_prs_neg = set(meth_exp_spm_neg['Gene']).intersection(meth_exp_prs_neg['Gene'])
        len(common_meth_exp_spm_prs_pos)  # 1066
        len(common_meth_exp_spm_prs_neg)  # 5790
        # Note: any gene with at least one promoter which methylation correlates with gene expression was selected.

        #### reactions/proteins which flux/protein usage correlates with cell growth rate:
        # get dataframe with fluxes/protein usage for reactions in cell line-specific models
        # even reactions not present in cell specific model are in this dataframe (just have flux/protein usage of 0):
        simul_df = pd.read_csv(os.path.join(fld_pth, fl_nm), index_col=0, sep='\t').T
        #  get only the simulated fluxes of all reactions except biomass:
        if 'flx' in fl_nm:
            ndf2 = simul_df.drop('adaptbiomass', axis=1)
        else:
            ndf2 = copy.deepcopy(simul_df)
            ndf2 = ndf2.T.groupby('Enzyme').sum().T
            # ndf2=ndf2.loc[:, set(sim_vars)-{'adaptbiomass'}]
        # get biomass:
        if 'flx' in fl_nm:
            biomass = simul_df['adaptbiomass'].drop('Subsystem')  # exclude 'Subsystem'
        else:

            # df_w_bm = pd.read_csv(os.path.join(fld_pth, fl_nm_bm), index_col=0, sep='\t')
            df_w_bm = pd.read_csv(os.path.join(fld_pth, kwargs['fl_nm_bm']), index_col=0, sep='\t')
            biomass = df_w_bm.loc['adaptbiomass'].drop('Subsystem')
            ndf2.drop('Reaction', axis=0, inplace=True)
        ndf2 = ndf2.loc[biomass.index]
        if 'flx' in fl_nm:
            # with spearman correlation:
            flx_cg_spm = ndf2.apply(lambda x: spearmanr(x, biomass), axis=0)  # NaN is when there is no variation
            flx_cg_spm.index = ['Coefficient', 'P-value']
            flx_cg_spm.loc['Subsystem'] = simul_df.drop('adaptbiomass', axis=1).loc['Subsystem'].values
            # with pearson correlation:
            flx_cg_prs = ndf2.apply(lambda x: pearsonr(x, biomass), axis=0)  # NaN is when there is no variation
            flx_cg_prs.index = ['Coefficient', 'P-value']
            flx_cg_prs.loc['Subsystem'] = simul_df.drop('adaptbiomass', axis=1).loc['Subsystem'].values
        else:
            # with spearman correlation:
            flx_cg_spm = ndf2.apply(lambda x: spearmanr(x, biomass), axis=0)  # NaN is when there is no variation
            flx_cg_spm.index = ['Coefficient', 'P-value']
            # with pearson correlation:
            flx_cg_prs = ndf2.apply(lambda x: pearsonr(x, biomass), axis=0)  # NaN is when there is no variation
            flx_cg_prs.index = ['Coefficient', 'P-value']
        # select with p-value < 0.05:
        flx_cg_spm = flx_cg_spm.loc[:, flx_cg_spm.loc['P-value'] < 0.05].T  # 605 reactions flux / 282 enzymes prot usage
        flx_cg_prs = flx_cg_prs.loc[:, flx_cg_prs.loc['P-value'] < 0.05].T  # 720 reactions flux / 289 enzymes prot usage
        # intersect positively/negatively correlated of spearman and pearson methods:
        flx_cg_spm_pos = flx_cg_spm.loc[flx_cg_spm['Coefficient'] > 0]
        flx_cg_spm_neg = flx_cg_spm.loc[flx_cg_spm['Coefficient'] < 0]
        flx_cg_prs_pos = flx_cg_prs.loc[flx_cg_prs['Coefficient'] > 0]
        flx_cg_prs_neg = flx_cg_prs.loc[flx_cg_prs['Coefficient'] < 0]
        common_flx_cg_spm_prs_pos = set(flx_cg_spm_pos.index).intersection(flx_cg_prs_pos.index)  # 286 enzymes prot usage
        common_flx_cg_spm_prs_neg = set(flx_cg_spm_neg.index).intersection(flx_cg_prs_neg.index)  # 0 enzymes prot usage
        ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
        if 'prot' in fl_nm:
            # convert enzymes to corresponding draw reactions:
            common_flx_cg_spm_prs_pos = {r.id for mt in common_flx_cg_spm_prs_pos for r in ms.metabolites.get_by_id(mt).reactions if 'draw_' in r.id}
            common_flx_cg_spm_prs_neg = {r.id for mt in common_flx_cg_spm_prs_neg for r in ms.metabolites.get_by_id(mt).reactions if 'draw_' in r.id}
        # associate corresponding genes:
        # for a gene with a correlation between mehylation-expression and methylation-growth rate,
        # associated with a reaction which flux is correlated with growth rate,
        # if that associated reaction has gene-reaction-rule with OR operators (the gene is a isozyme of a 'arm_' reaction):
        # - the gene IS selected, because methylation of any of the genes/isozymes could be affecting cell growth
        # if that associated reaction has gene-reaction-rule with AND operators (the gene is a enzyme of a complex):
        # - the gene IS selected, because methylation of any of the genes/complex subunits could be affecting cell growth
        r_gd = {r.id: [g.id for g in r.genes] for r in ms.reactions}
        common_flx_cg_spm_prs_pos = {sbel for el in [r_gd[r] for r in common_flx_cg_spm_prs_pos] for sbel in el}
        common_flx_cg_spm_prs_neg = {sbel for el in [r_gd[r] for r in common_flx_cg_spm_prs_neg] for sbel in el}

        # convert from ensembl id to gene symbols:
        def ens_smb(st):
            '''
                        - to convert ensembl to gene symbol. Note that the same ensembl might have more than one symbol
                        and vice-versa. that's why this function is needed
                        :param st: set with ensemb ids
                        :return: set with gene symbols
                        '''
            ens_sym_cg = dict()
            for gid in st:
                for s, ens in symbol_ensembl:
                    if gid == ens:
                        if ens in ens_sym_cg:
                            ens_sym_cg[ens].append(s)
                        else:
                            ens_sym_cg[ens] = [s]
            return {sel for el in st if el in ens_sym_cg for sel in ens_sym_cg[el]}

        common_flx_cg_spm_prs_pos = ens_smb(common_flx_cg_spm_prs_pos)
        common_flx_cg_spm_prs_neg = ens_smb(common_flx_cg_spm_prs_neg)
        len(common_flx_cg_spm_prs_pos)  # 619 genes of reactions flux / 262 genes of enzymes prot usage
        len(common_flx_cg_spm_prs_neg)  # 2 genes of reactions flux / 0 genes of enzymes prot usage

        #### genes which methylation correlates with growth rate:
        colint = list(set(meth_df.columns).intersection(set(biomass.index)))
        meth_sb = meth_df[colint]
        biomass_sb = biomass.loc[colint]
        # with spearman correlation:
        meth_cg_spm = pd.DataFrame()
        meth_cg_spm['Coefficient'], meth_cg_spm['P-value'] = meth_sb.apply(lambda x: spearmanr(x, pd.DataFrame(biomass_sb)), axis=1).str  # NaN is when there is no variation
        meth_cg_spm = meth_cg_spm.T
        # with pearson correlation:
        meth_cg_prs = pd.DataFrame()
        meth_cg_prs['Coefficient'], meth_cg_prs['P-value'] = meth_sb.apply(lambda x: pearsonr(x, pd.DataFrame(biomass_sb)), axis=1).str  # NaN is when there is no variation
        meth_cg_prs = meth_cg_prs.T.astype(float)
        # select with p-value < 0.05:
        meth_cg_spm = meth_cg_spm.loc[:, meth_cg_spm.loc['P-value'] < 0.05].T  #
        meth_cg_prs = meth_cg_prs.loc[:, meth_cg_prs.loc['P-value'] < 0.05].T  #
        # add gene symbol:
        meth_cg_spm['Gene'] = meth_cg_spm.index.to_series().apply(lambda x: x.split('_')[0])
        meth_cg_prs['Gene'] = meth_cg_prs.index.to_series().apply(lambda x: x.split('_')[0])
        # intersect positively/negatively correlated of spearman and pearson methods:
        meth_cg_spm_pos = meth_cg_spm.loc[meth_cg_spm['Coefficient'] > 0]
        meth_cg_spm_neg = meth_cg_spm.loc[meth_cg_spm['Coefficient'] < 0]
        meth_cg_prs_pos = meth_cg_prs.loc[meth_cg_prs['Coefficient'] > 0]
        meth_cg_prs_neg = meth_cg_prs.loc[meth_cg_prs['Coefficient'] < 0]
        common_meth_cg_spm_prs_pos = set(meth_cg_spm_pos['Gene']).intersection(meth_cg_prs_pos['Gene'])
        common_meth_cg_spm_prs_neg = set(meth_cg_spm_neg['Gene']).intersection(meth_cg_prs_neg['Gene'])
        len(common_meth_cg_spm_prs_pos)  # 1071
        len(common_meth_cg_spm_prs_neg)  # 26
        # Note: any gene with at least one promoter which methylation correlates with growth rate was selected.

        ### Intersect all and save:
        if 'flx' in fl_nm:
            vr = 'rc_flx'
        else:
            vr = 'prt_usg'
        ffld = os.path.join(fld_pth, 'meth_biom_trgt', f'{vr}')
        if not os.path.exists(ffld):
            os.makedirs(ffld)
        al = {'methexppos': common_meth_exp_spm_prs_pos, 'methexpneg': common_meth_exp_spm_prs_neg}
        bl = {'flxcgpos': common_flx_cg_spm_prs_pos, 'flxcgneg': common_flx_cg_spm_prs_neg}
        cl = {'methcgpos': common_meth_cg_spm_prs_pos, 'methcgneg': common_meth_cg_spm_prs_neg}
        for na, a in al.items():
            for nb, b in bl.items():
                for nc, c in cl.items():
                    ffv = a.intersection(b).intersection(c)
                    nm = f'{na}_{nb}_{nc}'
                    pd.DataFrame(ffv).to_csv(os.path.join(ffld, f'{nm}.tsv'), sep='\t', header=False, index=False)
