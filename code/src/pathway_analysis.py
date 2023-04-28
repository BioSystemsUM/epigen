import os
from graphs import Graphs
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import pandas as pd
import numpy as np
from test_mds import TestMds
import json
from cobra.io import load_matlab_model, load_yaml_model
import re
import math
from scipy.stats import t
import biomart
import json
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


class PathAnalysis:
    '''class with methods to analyze flux and protein usage in different pathways or reactions/enzymes'''
    @staticmethod
    def clust_box_plts(f, nm_clst, nm_bx, nm_bx_2, ylabel, htmp_log, trsf, mth_corr_fld_path,
                       cell_width, cell_height, fn, xlab_siz, ylab_siz, xlab_rot,
                       xcbar, ycbar, cbar_width, cbar_height, xheat, yheat, heat_wf, heat_hf,
                       rden_xf, cden_xf, rden_yf, cden_yf, rden_wf, cden_wf, rden_hf, cden_hf,
                       cbar_lab_siz, rclust, cclust):
        '''
        - creates boxplots and clustermaps for function 'analyze_pathway_flux'
        :param f: dataframe
        :param nm_clst: name of file where o save clustermap
        :param nm_bx: name of file where to save boxplots
        :param nm_bx_2: name of file where to save boxplots with top most flux across all cell lines/tissues
        :param ylabel: name of vertical axis in boxplot
        :param htmp_log: bool indicating whether to logarithmize values for the heatmap
        :param trsf: which of 'min_max' normalization, 'std_scaling'  (i.e. standard scaling) and 'none' to apply to rows of heatmap
        :param mth_corr_fld_path: path to analysis folder
        # other params used in applyclustermap
        '''
        # exclude weird subsystems:
        mett = ['Isolated', 'Miscellaneous', 'Artificial reactions']
        for x in mett:
            if x in f.index:
                f = f.drop(index=x)
        plt_d = f.copy(deep=True)
        plt_r = f.copy(deep=True)
        # logaritmize if required:
        if htmp_log:
            f.replace(to_replace=0, value=np.sort(np.unique(f))[1], inplace=True)  # replace 0 by second lowest value because log(0) is indefined
            f = np.log10(f)
        # subset data with the top 20 subsystems with highest variance between tissues or cell lines:
        top = np.var(f, axis=1).sort_values(ascending=False)[:20].index
        f = f.loc[top]
        # apply min-max normalization to the rows/subsystems if required:
        if trsf == 'min_max':
            f = f.T
            scaler = MinMaxScaler()
            scaled = scaler.fit_transform(f)
            f = pd.DataFrame(scaled, columns=f.columns, index=f.index)
            f = f.T
            f.index.name = ' '
        # apply standard scaling to the rows/subsystems if required:
        elif trsf == 'std_scaling':
            f = f.T
            scaler = StandardScaler()
            scaled = scaler.fit_transform(f)
            f = pd.DataFrame(scaled, columns=f.columns, index=f.index)
            f = f.T
            f.index.name = ' '
        pth = os.path.join(mth_corr_fld_path, nm_clst)
        f.columns = map(lambda x: x.lower().replace('_', ' '), f.columns)
        Graphs.applyclustermap(data=f, file_pth=pth, cell_width=cell_width, cell_height=cell_height, fn=fn,
                               xlab_siz=xlab_siz, ylab_siz=ylab_siz, xlab_rot=xlab_rot,
                               xcbar=xcbar, ycbar=ycbar, cbar_width=cbar_width, cbar_height=cbar_height, xheat=xheat,
                               yheat=yheat, heat_wf=heat_wf, heat_hf=heat_hf,
                               rden_xf=rden_xf, cden_xf=cden_xf, rden_yf=rden_yf, cden_yf=cden_yf, rden_wf=rden_wf,
                               cden_wf=cden_wf, rden_hf=rden_hf, cden_hf=cden_hf,
                               cbar_lab_siz=cbar_lab_siz, rclust=rclust, cclust=cclust)
        # get flux in specific pathways for all models or tissues:
        plt_d['Subsystem'] = plt_d.index
        subss = ['Glycolysis / Gluconeogenesis', 'Oxidative phosphorylation', 'Pentose phosphate pathway',
                 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism', 'Folate metabolism',
                 'Cysteine and methionine metabolism', 'Dna (de)/methylation']
        plt_d = plt_d.loc[subss]
        plt_n = pd.melt(plt_d, id_vars=['Subsystem'], value_vars=list(plt_d.columns[:-1]), value_name=ylabel)
        plt_n.index = plt_n['Subsystem']
        rnm_n = {subss[0]: 'Glycolysis or\ngluconeogenesis', subss[1]: 'Oxidative\nphosphorylation',
                 subss[2]: 'Pentose\nphosphate\npathway',
                 subss[3]: 'TCA and glyoxylate\nor dicarboxylate\nmetabolism', subss[4]: 'Folate\nmetabolism',
                 subss[5]: 'Cysteine\nand methionine\nmetabolism', subss[6]: 'DNA methylation\nor demethylation'}
        plt_n['Subsystem'] = plt_n['Subsystem'].apply(lambda x: rnm_n[x])
        pth = os.path.join(mth_corr_fld_path, nm_bx)
        Graphs.plotbxplt(data=plt_n, pth=pth, ylabel=ylabel)
        # get flux in top 5 pathways with most flux for all cell lines/tissues:
        topr = plt_r.median(axis=1).sort_values(ascending=False)[:5].index
        plt_r = plt_r.loc[topr]
        plt_r['Subsystem'] = plt_r.index
        plt_r = pd.melt(plt_r, id_vars=['Subsystem'], value_vars=list(plt_r.columns[:-1]), value_name=ylabel)
        plt_r.index = plt_r['Subsystem']
        rnm_r = {'Glycolysis / Gluconeogenesis': 'Glycolysis or\ngluconeogenesis',
                 'Fructose and mannose metabolism': 'Fructose and\nmannose metabolism',
                 'Cysteine and methionine metabolism': 'Cysteine\nand methionine\nmetabolism',
                 'Oxidative phosphorylation': 'Oxidative\nphosphorylation',
                 'Purine metabolism': 'Purine\nmetabolism', 'Cholesterol biosynthesis 2': 'Cholesterol\nbiosynthesis 2',
                 'Glycerophospholipid metabolism': 'Glycerophospholipid\nmetabolism',
                 'Cholesterol biosynthesis 1 (Bloch pathway)': 'Cholesterol\nbiosynthesis 1\n(Bloch pathway)',
                 'Aminoacyl-tRNA biosynthesis': 'Aminoacyl-tRNA\nbiosynthesis',
                 'Cholesterol metabolism': 'Cholesterol\nmetabolism',
                 'Fatty acid biosynthesis (even-chain)': 'Fatty acid\nbiosynthesis\n(even-chain)',
                 'Acylglycerides metabolism': 'Acylglycerides\nmetabolism'}
        plt_r['Subsystem'] = plt_r['Subsystem'].apply(lambda x: rnm_r[x])
        pth = os.path.join(mth_corr_fld_path, nm_bx_2)
        Graphs.plotbxplt(data=plt_r, pth=pth, ylabel=ylabel)

    @staticmethod
    def analyze_pathway_flux(mth_corr_fld_path, algo, obj_id, constr_ex, gen_md_pth, gen_gecko_pth,
                             methlt_fl_rc, with_tsk, htmp_log, per_tissue, op, trsf,
                             cell_width, cell_height, fn, xlab_siz, ylab_siz,
                             xlab_rot, xcbar, ycbar, cbar_width, cbar_height, xheat,
                             yheat, heat_wf, heat_hf, rden_xf, cden_xf, rden_yf, cden_yf, rden_wf,
                             cden_wf, rden_hf, cden_hf, cbar_lab_siz, rclust, cclust,
                             gr_const, cl_spec):
        '''
        - creates clustermaps with mean flux values of each subsystem (top 20 with highest variance) for tissues/cell lines
        - creates boxplot with flux values of specific subsystems from all models when 'per_tissue' is False or from all tissues when 'per_tissue' is True
        :param mth_corr_fld_path: path to analysis folder
        :param algo: algorithm used during model reconstruction
        :param obj_id: objective used in the flux simulation (dict)
        :param constr_ex: either a list with exometabolites which fluxes are constraint or 'False' indicating no constraints are applied to exometabolite fluxes
        :param gen_md_pth: path to generic traditional GSMM (non-GECKO)
        :param gen_gecko_pth: path to generic gecko model
        :param methlt_fl_rc: path to .json file with added reactions related with DNA methylation
        :param with_tsk: bool indicating if cell-specific tasks were included during reconstruction
        :param htmp_log: bool indicating whether to logarithmize values for the heatmap
        :param per_tissue: bool indicating whether to aggregate values by tissue or not (keep individual cell lines)
        :param op: bool indicating whether to aggregate values of same subsystem by mean or sum
        :param trsf: which of 'min_max' normalization, 'std_scaling'  (i.e. standard scaling) and 'none' to apply to rows of heatmap
        :param gr_const: whether cell growth was constraint with experimental growth rates
        :param cl_spec: bool, whether cell line specific reaction of 'prodDNAtot' was used (True) or the generic was used instead (False).
        # other params used in applyclustermap
        '''
        if with_tsk: mth_corr_fld_path = os.path.join(mth_corr_fld_path, algo, 'including_tsks')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, algo, 'no_tsks')
        if constr_ex: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'w_flux_constr')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'no_flux_constr')
        if gr_const: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'biomass_constr')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'no_biomass_constr')
        if cl_spec: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'cl_spec_DNAtot')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'generic_DNAtot')
        file_lst = TestMds.create_file_list(path=mth_corr_fld_path, incl_wrd='fluxes')
        file_lst = [f for f in file_lst if list(obj_id.keys())[0] in f]
        if type(constr_ex) == list:
            file_lst = [f for f in file_lst if '-'.join(constr_ex) in f]  # use only the files corresponding to the constraints we are testing
        file_pth = file_lst[0]
        df = pd.read_csv(file_pth, sep='\t', index_col=0)
        with open(methlt_fl_rc) as f:
            meth_rc = json.load(f)
        ## get flux of each pathway (sum of flux of all active reactions of a pathway divided by the number of reactions of that pathway in generic model):
        m = load_yaml_model(gen_md_pth)  # generic traditional model
        ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
        sb_d = {r.id: r.subsystem[0] for r in m.reactions}
        gen_gecko_r = [r.id for r in ms.reactions]
        drawpt_ex = [r for r in gen_gecko_r if 'draw_' in r or 'prot' in r]  # list of all draw_protein reactions and the 'prot_pool_exchange' of the generic gecko model
        genr_nodraw = set(gen_gecko_r) - set(drawpt_ex)  # exclude 'draw_prot' reactions and 'prot_pool_exchange'
        one_enz_r = [r for r in genr_nodraw if r.endswith('No1') and 'arm_' + r[:-3] not in genr_nodraw]  # list of remaining reactions with one enzyme only (because they do not have a corresponding 'arm_' reaction)
        arm_r = [r for r in genr_nodraw if 'arm_' in r]  # list of remaining reactions with many enzymes (with 'arm_' reactions)
        no_enz_r = [r for r in genr_nodraw if 'No' not in r and 'arm_' not in r]  # list of remaining reactions with no associated enzyme (non-catalyzed reactions)
        impt_r = set(one_enz_r).union(set(arm_r)).union(set(no_enz_r))  # all reactions except 'draw_prot', 'prot_pool_exchange', and leg reactions corresponding to arm_reactions (i.e. isozymes of arm_reactions)
        excp = ['prodDNAtot', 'adaptbiomass']  # artificial reactions - pseudo-reactions
        mthrc = [r for r in (set(meth_rc) - set(excp)) if 'transp' not in r] + ['MAR08641']  # all DNA methylation reactions except transport and artificial (pseudo-reactions)
        # get dictionary with reaction id vs subsystem in generic gecko model, for all reactions except
        # 'draw_prot', 'prot_pool_exchange', and leg reactions corresponding to arm_reactions (i.e. isozymes of arm_reactions):
        sb_gek_d = {r: 'Transport reactions' if 'transp' in r else 'Artificial reactions' if r in excp else 'Dna (de)/methylation' if [r for mrc in mthrc if mrc in r] else sb_d[r.split('_')[1]] if r.startswith('arm_') else sb_d[r.split('_')[0]] if '_REV' in r else sb_d[r[:-3]] if r.endswith('No1') else sb_d[r] for r in impt_r}
        # from the simulated fluxes keep those of all reactions except 'draw_prot', 'prot_pool_exchange', and leg reactions corresponding to arm_reactions (i.e. isozymes of arm_reactions)
        otf = df.loc[set(df.index).intersection(impt_r)]
        # if True get mean value of each reaction per tissue, otherwise its per cell line:
        if per_tissue:
            otf.columns = ['_'.join(col.split('_')[1:]) for col in otf.columns]
            tis = list()
            for el in otf.columns:
                if el not in tis:
                    tis.append(el)
            lst = list()
            for ts in tis:
                r = otf[ts].median(axis=1)
                r.name = ts
                lst.append(r)
            otf = pd.concat(lst, axis=1)
        # add subsystem to each reaction of dataframe with simulated values:
        otf['Subsystem'] = [sb_gek_d[r] for r in otf.index]
        # subsystems of reactions (except 'draw_prot', 'prot_pool_exchange', and leg reactions corresponding to arm_reactions) of generic model that do not exist in otf dataframe:
        nf = pd.DataFrame.from_dict({k: [v] for k, v in sb_gek_d.items() if k not in otf.index}).T
        nf.columns = ['Subsystem']
        # create dataframe with simulated values and zero fluxes for all reactions
        # in the generic gecko that do not exist in any context-specific model:
        flx = pd.concat([otf, pd.concat([pd.DataFrame(np.zeros((len(nf.index), len(otf.columns[:-1]))), index=nf.index, columns=otf.columns[:-1]), nf], axis=1)])
        fld_allflx = os.path.join(mth_corr_fld_path, 'corr_pth_mth', )
        if not os.path.exists(fld_allflx): os.makedirs(fld_allflx)
        flx.to_csv(os.path.join(fld_allflx, 'allflx_rc_substm.tsv'), sep='\t')
        # use 'groupby' to get the 'mean flux' per each subsystem across different tissues or cell lines:
        if op == 'mean':
            plt_f = flx.groupby('Subsystem').mean()
            ylabel = 'Average Flux (mmol/gDW.h)'
        elif op == 'sum':
            plt_f = flx.groupby('Subsystem').sum()
            ylabel = 'Total Flux (mmol/gDW.h)'
        PathAnalysis.clust_box_plts(f=plt_f, nm_clst='var_heatmap.svg', nm_bx='boxplot.svg', nm_bx_2='boxplot_top5path.svg',
                       ylabel=ylabel, htmp_log=htmp_log, trsf=trsf, mth_corr_fld_path=mth_corr_fld_path,
                       cell_width=cell_width, cell_height=cell_height, fn=fn, xlab_siz=xlab_siz, ylab_siz=ylab_siz,
                       xlab_rot=xlab_rot, xcbar=xcbar, ycbar=ycbar, cbar_width=cbar_width, cbar_height=cbar_height, xheat=xheat,
                       yheat=yheat, heat_wf=heat_wf, heat_hf=heat_hf, rden_xf=rden_xf, cden_xf=cden_xf, rden_yf=rden_yf, cden_yf=cden_yf, rden_wf=rden_wf,
                       cden_wf=cden_wf, rden_hf=rden_hf, cden_hf=cden_hf, cbar_lab_siz=cbar_lab_siz, rclust=rclust, cclust=cclust)
        # flux histogram distribution:
        Graphs.dist_plt(frm=plt_f, fld=mth_corr_fld_path, fnm='flx_distribution.svg')

    @staticmethod
    def analyze_pathway_protein(mth_corr_fld_path, algo, obj_id, constr_ex, gen_md_pth, gen_gecko_pth,
                                methlt_fl_rc, with_tsk, htmp_log, per_tissue, op, trsf,
                                cell_width, cell_height, fn, xlab_siz, ylab_siz,
                                xlab_rot, xcbar, ycbar, cbar_width, cbar_height, xheat,
                                yheat, heat_wf, heat_hf, rden_xf, cden_xf, rden_yf, cden_yf, rden_wf,
                                cden_wf, rden_hf, cden_hf, cbar_lab_siz, rclust, cclust,
                                gr_const, cl_spec):
        '''
            - creates clustermaps with mean protein usage of each subsystem (top 20 with highest variance) for tissues/cell lines
            - creates boxplot with protein usage of specific subsystems from all models when 'per_tissue' is False or from all tissues when 'per_tissue' is True
            - protein usage units are in mg/gDW
            :param mth_corr_fld_path: path to analysis folder
            :param algo: algorithm used during model reconstruction
            :param obj_id: objective used in the flux simulation (dict)
            :param constr_ex: either a list with exometabolites which fluxes are constraint or 'False' indicating no constraints are applied to exometabolite fluxes
            :param gen_md_pth: path to generic traditional GSMM (non-GECKO)
            :param gen_gecko_pth: path to generic gecko model
            :param methlt_fl_rc: path to .json file with added reactions related with DNA methylation
            :param with_tsk: bool indicating if cell-specific tasks were included during reconstruction
            :param htmp_log: bool indicating whether to logarithmize values for the heatmap
            :param per_tissue: bool indicating whether to aggregate values by tissue or not (keep individual cell lines)
            :param op: bool indicating whether to aggregate values of same subsystem by mean or sum
            :param trsf: which of 'min_max' normalization, 'std_scaling'  (i.e. standard scaling) and 'none' to apply to rows of heatmap
            :param gr_const: whether cell growth was constraint with experimental growth rates
            :param cl_spec: bool, whether cell line specific reaction of 'prodDNAtot' was used (True) or the generic was used instead (False).
            # other params used in applyclustermap
        '''
        if with_tsk: mth_corr_fld_path = os.path.join(mth_corr_fld_path, algo, 'including_tsks')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, algo, 'no_tsks')
        if constr_ex: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'w_flux_constr')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'no_flux_constr')
        if gr_const: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'biomass_constr')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'no_biomass_constr')
        if cl_spec: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'cl_spec_DNAtot')
        else: mth_corr_fld_path = os.path.join(mth_corr_fld_path, 'generic_DNAtot')
        file_lst = TestMds.create_file_list(path=mth_corr_fld_path, incl_wrd='fluxes')
        file_lst = [f for f in file_lst if list(obj_id.keys())[0] in f]
        if type(constr_ex) == list:
            file_lst = [f for f in file_lst if '-'.join(constr_ex) in f]  # use only the files corresponding to the constraints we are testing
        file_pth = file_lst[0]
        print(file_pth)
        df = pd.read_csv(file_pth, sep='\t', index_col=0)
        with open(methlt_fl_rc) as f:
            meth_rc = json.load(f)
        ms = load_matlab_model(gen_gecko_pth)  # generic gecko model
        # get correspondence between reaction and enzyme id for generic gecko model - gets rid of 'arm' and non-catalyzed reactions:
        fr = pd.DataFrame([(r.id, mt.id) for r in ms.reactions for mt in r.metabolites if 'prot' in mt.id])
        fr.index = fr[0]
        fr = fr.drop(0, axis=1)
        # get rid of 'draw_prot' and 'prot_pool_exchange' reactions:
        texc = [r for r in fr.index if r.startswith('draw_')]  # explanation below of why not use draw_protein reactions instead
        ofr = fr.drop(texc + ['prot_pool_exchange'])
        ofr = ofr.rename(columns={1: 'Enzyme'})
        # get reaction vs enzyme vs subsystem for generic traditional GSMM:
        m = load_yaml_model(gen_md_pth)  # generic traditional model
        sb_d = {r.id: r.subsystem[0] for r in m.reactions}
        excp = ['prodDNAtot', 'adaptbiomass']  # artificial reactions - pseudo-reactions
        mthrc = [r for r in (set(meth_rc) - set(excp)) if 'transp' not in r] + ['MAR08641']  # all DNA methylation reactions except transport and artificial (pseudo-reactions)
        sb_gek_d = {r: 'Transport reactions' if 'transp' in r else 'Artificial reactions' if r in excp else 'Dna (de)/methylation' if [r for mrc in mthrc if mrc in r] else sb_d[r.split('_')[0]] if '_REV' in r else sb_d[re.sub(r'No\d+', '', r)] if bool(re.search(r'No\d+', r)) else sb_d[r] for r in ofr.index}
        nd3 = ofr.join(pd.DataFrame(sb_gek_d, index=['Subsystem']).T)
        # - note: same enzyme may be used by different reactions of different subsystems, e.g. nd3.loc[nd3['Enzyme'] == 'prot_A0A087WXM9'],
        #         so, we can't use directly the flux of draw reaction of each enzyme.
        # get molecular weight (MW) of enzymes:
        mwd = pd.DataFrame({e: -(ms.reactions.get_by_id('draw_' + e).get_coefficient('prot_pool')) for e in nd3['Enzyme']}, index=['MMW']).T
        mwd['Enzyme'] = mwd.index
        # get Kcat of enzymes and merge with the MWs dataframe:
        rs = pd.DataFrame.from_dict({r: {m.id: -(1 / c) for m, c in ms.reactions.get_by_id(r).metabolites.items() if 'prot_' in m.id} for r in nd3.index})
        rs['Enzyme'] = rs.index
        rsf = pd.melt(rs, id_vars=['Enzyme'], value_vars=rs.columns[:-1], var_name='Reaction', value_name='Kcat').dropna()
        mrg = pd.merge(rsf, mwd, how='outer')
        # get dataframe with reaction, enzyme, subsystem, MW and Kcat of enzyme:
        nd3['Reaction'] = nd3.index
        mrg2 = pd.merge(mrg, nd3, how='outer', on=['Reaction', 'Enzyme'])
        # if True get mean value of each reaction per tissue, otherwise its per cell line:
        df2 = df.copy(deep=True)
        if per_tissue:
            df2.columns = ['_'.join(col.split('_')[1:]) for col in df2.columns]
            tis = list()
            for el in df2.columns:
                if el not in tis:
                    tis.append(el)
            lst = list()
            for ts in tis:
                r = df2[ts].mean(axis=1)
                r.name = ts
                lst.append(r)
            df2 = pd.concat(lst, axis=1)
        df2['Reaction'] = df2.index
        # create dataframe with simulated values and zero fluxes for all reactions
        # in the generic gecko that do not exist in any context-specific model:
        mrg3 = pd.merge(mrg2, df2, on='Reaction', how='left')  # how='left' to get rid of 'draw_prot'reactions, 'prot_pool_exchange' and 'arm_' reactions in df2, and to create space for reactions not in any reconstructed model
        mrg3.fillna(0, inplace=True)  # NAs here are reactions from generic model that do not exist in any tissue specific model
        # get protein usage in each individual reaction-enzyme combination in mg/gDW:
        divs = (mrg3.iloc[:, 5:]).div(mrg3['Kcat'], axis=0)
        wtfrc = divs.multiply(mrg3['MMW'], axis=0) * 1000 # note that MW is in KDa = 1000g/mol = 1g/mmol, and we want result in mg/gDW
        frc = pd.concat([mrg3.iloc[:, :5], wtfrc], axis=1)
        frc.drop(['Kcat', 'MMW'], axis=1, inplace=True)
        fld_allprot = os.path.join(mth_corr_fld_path, 'corr_pth_mth', )
        if not os.path.exists(fld_allprot): os.makedirs(fld_allprot)
        frc.to_csv(os.path.join(fld_allprot, 'allprot_rc_substm.tsv'), sep='\t')
        # use 'groupby' to get the 'mean flux' per each subsystem across different tissues or cell lines.
        # Note that we need to divide by the number of reaction-enzyme combinations,
        # to correct for the bias that subsystems with more reactions
        # and with reactions containing more enzymes have the tendency to use more protein:
        if op == 'mean':
            plt_f = frc.groupby('Subsystem').mean()  # 'mean()' divides by number of combinations (enzyme, reaction, subsystem)
            ylabel = 'Average Protein (mg/gDW)'
        elif op == 'sum':
            plt_f = frc.groupby('Subsystem').sum()
            ylabel = 'Total Protein (mg/gDW)'
        PathAnalysis.clust_box_plts(f=plt_f, nm_clst='prot_var_heatmap.svg', nm_bx='prot_boxplot.svg',
                       nm_bx_2='prot_boxplot_top5path.svg', ylabel=ylabel, htmp_log=htmp_log, trsf=trsf,
                       mth_corr_fld_path=mth_corr_fld_path, cell_width=cell_width, cell_height=cell_height, fn=fn, xlab_siz=xlab_siz, ylab_siz=ylab_siz,
                       xlab_rot=xlab_rot, xcbar=xcbar, ycbar=ycbar, cbar_width=cbar_width, cbar_height=cbar_height, xheat=xheat,
                       yheat=yheat, heat_wf=heat_wf, heat_hf=heat_hf, rden_xf=rden_xf, cden_xf=cden_xf, rden_yf=rden_yf, cden_yf=cden_yf, rden_wf=rden_wf,
                       cden_wf=cden_wf, rden_hf=rden_hf, cden_hf=cden_hf, cbar_lab_siz=cbar_lab_siz, rclust=rclust, cclust=cclust)
        # mass histogram distribution:
        Graphs.dist_plt(frm=plt_f, fld=mth_corr_fld_path, fnm='mass_distribution.svg')

    @staticmethod
    def r2p(r, N, two_tails=True):
        """
        Calculates the p value from a coefficient correlation.
        :param (float) r: the coefficient correlation value
        :param (int) N: Sample size
        :param (bool) two_tails: 1 or 2 tails. Default 2 tails.
        """
        r = np.abs(r)
        x = r * math.sqrt(N - 2) / math.sqrt(1 - r * r)
        p = t.sf(x, df=N - 2)
        if two_tails:
            return 2 * p
        else:
            return p
    @staticmethod
    def get_ensembl_mappings(fl_json):
        '''
        - get dictionary with mapping of ensembl gene ids to gene symbol
        - adapted from: 'https://autobencoder.com/2021-10-03-gene-conversion/'
        :param fl_json: path to .json file where to retrieve the mappings or where to save the mappings if the file does not exist yet
        '''
        if os.path.exists(fl_json): # if .json file already exists load it
            with open(fl_json, mode='r') as f:
                genesymbol_to_ensembl = json.loads(f.read())
        else:
            fld_sv = '/'.join(fl_json.split('/')[:-1])
            if not os.path.exists(fld_sv):
                os.makedirs(fld_sv)
            # set up connection to server:
            server = biomart.BiomartServer('http://useast.ensembl.org/biomart/martservice')
            mart = server.datasets['hsapiens_gene_ensembl']
            # list the types of data we want
            attributes = ['hgnc_symbol', 'ensembl_gene_id']
            # get the mapping between the attributes
            response = mart.search({'attributes': attributes})
            # to select directly the ids we need to convert we could do:
            # response = mart.search({'attributes': attributes, 'filters': {'hgnc_symbol': gene_symbol_lst}})
            # but if gene_symbol_lst is too big it complains it is too big of a request.
            # convert request results from a binary string to an easier-to-work-with text string:
            data = response.raw.data.decode('ascii')
            genesymbol_to_ensembl = list()
            for line in data.splitlines():
                line = line.split('\t')
                # the entries are in the same order as in the `attributes` variable:
                gene_symbol = line[0]
                ensembl_gene = line[1]
                if len(gene_symbol) > 0:  # if not empty
                    genesymbol_to_ensembl.append((gene_symbol, ensembl_gene))
            with open(fl_json, mode='w') as f:
                json.dump(genesymbol_to_ensembl, f)
        return genesymbol_to_ensembl
    @staticmethod
    def hypergemtest(sample_gene_sb, population_gene_sb, sample_sb):
        '''
        - do hypergeometric test and multiple test correction(FDR with Benjamin Hochberg)
        :param sample_gene_sb: dictionary {gene: [subsystem, ...]} where each gene is a gene of the sample
                               and each subsystem is subsystem of a reaction corresponding to taht gene in the generic model.
        :param population_gene_sb: dictionary {gene: [subsystem, ...]} where each gene is a gene of the population,
                               i.e. of the generic model, and each subsystem is subsystem of a reaction corresponding to taht gene in the generic model.
        :param sample_sb: set with subsystems of all reactions corresponding to the genes in the sample
        '''
        # get number of successes in sample for each pathway:
        ns_smp = dict()
        for v in sample_gene_sb.values():
            for el in v:
                if el in ns_smp:
                    ns_smp[el] += 1
                else:
                    ns_smp[el] = 1
        # get number of successes in population:
        ns_pop = dict()
        for v in population_gene_sb.values():
            for el in v:
                if el in ns_pop:
                    ns_pop[el] += 1
                else:
                    ns_pop[el] = 1
        # size of sample:
        N = len([el for v in sample_gene_sb.values() for el in v])
        # size of population:
        M = len([el for v in population_gene_sb.values() for el in v])
        # do the test:
        # "x-1" in the formula is explained at https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
        hyptest = dict()
        for sb in sample_sb:
            x = ns_smp[sb]  # number of successes in sample
            n = ns_pop[sb]  # number of successes in population
            pval = hypergeom.sf(x - 1, M, n, N)
            hyptest[sb] = [pval]
        hypfram = pd.DataFrame(hyptest).T
        hypfram.rename(columns={0: 'Pval'}, inplace=True)
        # apply multitest:
        hypfram['Pval'] = multipletests(pvals=hypfram['Pval'], method='fdr_bh')[1]
        # select significant:
        hypfram_sig = hypfram[hypfram['Pval'] < 0.05]
        # sort by pvalue:
        hypsig = hypfram_sig.sort_values('Pval', ascending=True)
        return hypsig

