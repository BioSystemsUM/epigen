import os

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats.mstats import spearmanr
import pandas as pd

class Graphs:
    '''class with methods to plot graphs'''

    @staticmethod
    def plots_biomass(df, grf, pth):
        '''
        - plots scatter plot with log10 values of simulated vs experimental growth rates of different cell lines
        - creates boxplot with relative errors between simulated and experimental growth rates of different cell lines
        :param df: dataframe with simulated biomass fluxes of several cell lines (cell lines on columns)
        :param grf: dataframe with experimental biomass fluxes of several cell lines (cell lines on columns)
        :param pth: path where to save the plots
        '''
        grf = grf['mean']
        grf = grf.loc[df.columns]
        df_plt = df.to_numpy().flatten()
        exp_plt = grf.to_numpy().flatten()
        # replace 0 by 2nd smallest value, because of log10(0) afterwards:
        rpl_v = min(np.sort(np.unique(df_plt))[1], np.sort(np.unique(exp_plt))[1])
        df_plt_f = np.where(df_plt==0.0, rpl_v, df_plt)
        exp_plt_f = np.where(exp_plt==0.0, rpl_v, exp_plt)
        plt_f = pd.DataFrame({'Predicted flux (log10 val.)': np.log10(df_plt_f), 'Measured flux (log10 val.)': np.log10(exp_plt_f)})
        def corrfunc(x, y, **kwargs):
            pearson, pp = pearsonr(x, y)
            spearman, ps = spearmanr(x, y)
            ax = plt.gca()
            ax.annotate(f'pearson = {round(pearson, 2)} \n p-value = {round(pp, 2)}',
                        xy=(.1, .9), xycoords=ax.transAxes, fontsize=13)
            ax.annotate(f'spearman = {round(spearman, 2)} \n p-value = {round(ps, 2)}',
                        xy=(.1, .7), xycoords=ax.transAxes, fontsize=13)
        ax = sns.scatterplot(x='Measured flux (log10 val.)', y='Predicted flux (log10 val.)', data=plt_f, s=10)
        corrfunc(x=plt_f['Measured flux (log10 val.)'], y=plt_f['Predicted flux (log10 val.)'])
        plt.xlabel(xlabel=r'Measured flux (log$\rm_{10}$ val.)', fontsize=13)
        plt.ylabel(ylabel=r'Predicted flux (log$\rm_{10}$ val.)', fontsize=14)
        mx = max(max(plt_f['Measured flux (log10 val.)']), max(plt_f['Predicted flux (log10 val.)']))
        mn = min(min(plt_f['Measured flux (log10 val.)']), min(plt_f['Predicted flux (log10 val.)']))
        plt.xlim(mn - 0.1, mx + 0.1)
        plt.ylim(mn - 0.1, mx + 0.1)
        plt.tight_layout()
        plt.savefig(pth, dpi=300)
        plt.close()
        print('####biomass####')
        dif = np.abs(plt_f['Predicted flux (log10 val.)'] - plt_f['Measured flux (log10 val.)'])
        print('biomass fluxes percentages:', np.sum(dif <= 1) / len(dif) * 100)
        # get relative error:
        plt_err = pd.DataFrame({'Predicted flux': df_plt, 'Measured flux': exp_plt})
        re = np.abs(plt_err['Predicted flux'] - plt_err['Measured flux'])/plt_err['Measured flux']
        sns.boxplot(y=re)
        plt.title('Relative error of growth rate', fontsize=18, pad=18)
        plt.yticks(fontsize=13)
        plt.tight_layout()
        pth_err = '/'.join(pth.split('/')[:-1])+'/rel_error_biomass_eleven.svg' if 'eleven' in pth.split('/')[-1] else '/'.join(pth.split('/')[:-1]) + '/rel_error_biomass.svg'
        plt.savefig(pth_err, dpi=300)
        plt.close()

    @staticmethod
    def plots_f(df, exp_df, pth):
        '''
        - plots scatter plot with log10 values of simulated vs experimental fluxes of some exchange reactions
        - plots flux values histograms with and without log10 of those simulated and experimental fluxes
        :param df: dataframe with simulated exchange fluxes of several cell lines (cell lines on columns)
        :param exp_df: dataframe with experimental exchange fluxes of several cell lines (cell lines on columns)
        :param pth: path where to save the plots
        '''
        df_i = df.to_numpy().flatten()
        exp_i = exp_df.to_numpy().flatten()
        df_plt_i = np.abs(df_i)
        exp_plt_i = np.abs(exp_i)
        # replace 0 by 2nd smallest value, to show log10(0) in the graph :
        rpl_v = min(np.sort(np.unique(df_plt_i))[1], np.sort(np.unique(exp_plt_i))[1])
        df_plt = np.where(df_plt_i == 0.0, rpl_v, df_plt_i)
        exp_plt = np.where(exp_plt_i == 0.0, rpl_v, exp_plt_i)
        plt_f = pd.DataFrame({'Predicted flux (log10 abs. val.)': np.log10(df_plt), 'Measured flux (log10 abs. val.)': np.log10(exp_plt)})
        def corrfunc(x, y, **kwargs):
            pearson, pp = pearsonr(x, y)
            spearman, ps = spearmanr(x, y)
            ax = plt.gca()
            ax.annotate(f'pearson = {round(pearson, 2)} \n p-value = {round(pp, 2)}',
                        xy=(.1, .9), xycoords=ax.transAxes, fontsize=13)
            ax.annotate(f'spearman = {round(spearman, 2)} \n p-value = {round(ps, 2)}',
                        xy=(.1, .7), xycoords=ax.transAxes, fontsize=13)
        ax = sns.scatterplot(x='Measured flux (log10 abs. val.)', y='Predicted flux (log10 abs. val.)', data=plt_f, s=4)
        corrfunc(x=plt_f['Measured flux (log10 abs. val.)'], y=plt_f['Predicted flux (log10 abs. val.)'])
        plt.xlabel(xlabel=r'Measured flux (log$\rm_{10}$ abs. val.)', fontsize=13)
        plt.ylabel(ylabel=r'Predicted flux (log$\rm_{10}$ abs. val.)', fontsize=14)
        mx = max(max(plt_f['Measured flux (log10 abs. val.)']), max(plt_f['Predicted flux (log10 abs. val.)']))
        mn = min(min(plt_f['Measured flux (log10 abs. val.)']), min(plt_f['Predicted flux (log10 abs. val.)']))
        plt.xlim(mn - 0.1, mx + 0.1)
        plt.ylim(mn - 0.1, mx + 0.1)
        plt.tight_layout()
        plt.savefig(pth, dpi=300)
        plt.close()
        print('####exchanges####')
        # sbs = plt_f.loc[(-4 < plt_f['Measured flux (log10 abs. val.)']) & (plt_f['Measured flux (log10 abs. val.)'] < 0), :]
        # subsbs = sbs.loc[(-4 < sbs['Predicted flux (log10 abs. val.)']) & (sbs['Predicted flux (log10 abs. val.)'] < 0), :]
        # print(subsbs.shape[0] / plt_f.shape[0] * 100) # 15.267695099818512 %
        dif = np.abs(plt_f['Predicted flux (log10 abs. val.)'] - plt_f['Measured flux (log10 abs. val.)'])
        print('exchange fluxes percentages:', np.sum(dif <=1)/len(dif) * 100)
        # histograms -log10:
        variables = list(plt_f.columns)
        # plt.figure(figsize=(20, 12))
        plt.figure(figsize=(17, 12))
        plt.subplots_adjust(hspace=0.5)
        for i, var in enumerate(variables):
            ax = plt.subplot(1, 2, i + 1)
            n, bins, edges = ax.hist(plt_f[var])
            varr = " ".join(var.split(' ')[:2]) + " " + r"(log$\rm_{10}$" + " " + "abs. val.)"
            ax.set_title(varr, fontsize=31, pad=20)
            # ax.set_xticks(np.round(bins[:len(bins):2], 1))
            ax.set_xticks(bins[:len(bins):2])
        plt.tight_layout(pad=20)
        # plt.show()
        pth_hist = '/'.join(pth.split('/')[:-1])+'/hist_log.svg'
        plt.savefig(pth_hist, dpi=300)
        plt.close()
        # histogram NOT log10 of abs val:
        plt_h = pd.DataFrame({'Predicted flux (abs. val.)': df_plt_i, 'Measured flux (abs. val.)': exp_plt_i})
        variables = list(plt_h.columns)
        # plt.figure(figsize=(20, 12))
        plt.figure(figsize=(10, 6))
        plt.subplots_adjust(hspace=0.5)
        for i, var in enumerate(variables):
            ax = plt.subplot(1, 2, i + 1)
            n, bins, edges = ax.hist(plt_h[var])
            ax.set_title(var, fontsize=27, pad=20)
            # ax.set_xticks(np.round(bins[:len(bins):2], 1))
            ax.set_xticks(bins[:len(bins):2])
        plt.tight_layout()
        pth_hist = '/'.join(pth.split('/')[:-1]) + '/hist.svg'
        plt.savefig(pth_hist, dpi=300)
        plt.close()

    @staticmethod
    def plot_exc_flx(exp_df, exc_flx, fld_fn):
        '''
        - prepare data for scatter plots of biomass reaction fluxes + boxplot with relative errors
        :param exp_df: dataframe with experimental values for exchange fluxes. rows are exchange reactions ids
        :param exc_flx: dict {model_id: {exchange_id: flux value}} for all models
        :param fld_fn: path to folder where to save results
        '''
        df = pd.DataFrame.from_dict(exc_flx)
        df.dropna(inplace=True)
        exp_df = exp_df.loc[df.index, df.columns]
        pth = os.path.join(fld_fn, 'exchanges_reconst_corr.svg')
        Graphs.plots_f(df=df, exp_df=exp_df, pth=pth)

    @staticmethod
    def prep_plot_biomass(mass_final_flx, fld_fn, grf):
        '''
        - prepare data for scatter plot of biomass reaction fluxes + boxplot with relative errors
        :param mass_final_flx: dict {model_id: {biomass_id: flux value}} for all models
        :param fld_fn: path to folder where to save results
        :param grf: dataframe with experimental biomass fluxes of several cell lines (cell lines on columns)
        '''
        df = pd.DataFrame.from_dict(mass_final_flx)
        pth = os.path.join(fld_fn, 'biomass_reconst_corr.svg')
        Graphs.plots_biomass(df=df, grf=grf, pth=pth)

    @staticmethod
    def applyclustermap(data, file_pth, cell_width, cell_height, fn, xlab_siz, ylab_siz, xlab_rot,
                        xcbar, ycbar, cbar_width, cbar_height, xheat, yheat, heat_wf, heat_hf,
                        rden_xf, cden_xf, rden_yf, cden_yf, rden_wf, cden_wf, rden_hf, cden_hf,
                        cbar_lab_siz, rclust, cclust):
        '''
        - does a clean clustermap graph
        :param data: data for heatmap
        :param file_pth: path to file where to save graph
        :param cell_width: width of each cell - controls figure size
        :param cell_height: height of each cell - controls figure size
        :param fn: param to increase overall font scale
        :param xlab_siz: fontsize of row labels
        :param ylab_siz: fontsize of column labels
        :param xlab_rot: degrees of row label rotation
        :param xcbar: x position of color bar
        :param ycbar: y position of color bar
        :param cbar_width: color bar width
        :param cbar_height: color bar height
        :param xheat: factor to multiply by current heatmap position to move it horizontally (> 1 increases, <1 decreases)
        :param yheat: factor to multiply by current heatmap position to move it vertically (> 1 increases, <1 decreases)
        :param heat_wf: factor to multiply by current heatmap width to to expand (>1) or shrink (<1) it horizontally
        :param heat_hf: factor to multiply by current heatmap height to to expand (>1) or shrink (<1) it vertically
        :param rden_xf: factor to multiply by current row dendogram position to move it horizontally (> 1 increases, <1 decreases)
        :param cden_xf: factor to multiply by current column dendogram position to move it horizontally (> 1 increases, <1 decreases)
        :param rden_yf: factor to multiply by current row dendogram position to move it vertically (> 1 increases, <1 decreases)
        :param cden_yf: factor to multiply by current column dendogram position to move it vertically (> 1 increases, <1 decreases)
        :param rden_wf: factor to multiply by current row dendogram width to to expand (>1) or shrink (<1) it horizontally
        :param rden_hf: factor to multiply by current row dendogram height to to expand (>1) or shrink (<1) it vertically
        :param cden_wf: factor to multiply by current column dendogram width to to expand (>1) or shrink (<1) it horizontally
        :param cden_hf: factor to multiply by current column dendogram height to to expand (>1) or shrink (<1) it vertically
        :param cbar_lab_siz: size of letters of colorbar label
        :param rclust: bool indicating whether to cluster rows or not
        :param cclust: bool indicating whether to cluster columns or not
        :return: saves the graph
        '''
        nc = data.shape[1]
        nr = data.shape[0]
        fcl = fn / nc
        sns.set(font_scale=fcl)
        colrs = sns.color_palette("dark:white_r", as_cmap=True)
        res = sns.clustermap(data, cmap=colrs, linecolor='black', xticklabels=True, yticklabels=True,
                             linewidths=0.1, figsize=(cell_width * nc, cell_height * nr),
                             cbar_pos=(xcbar, ycbar, cbar_width, cbar_height),
                             cbar_kws={'ticks': [data.to_numpy().flatten().min(), data.to_numpy().flatten().max()]},
                             row_cluster=rclust, col_cluster=cclust, rasterized=True)
        res.ax_cbar.tick_params(labelsize=cbar_lab_siz, pad=10)
        res.ax_cbar.set_yticklabels([str(np.round(data.to_numpy().flatten().min(), 1)), str(np.round(data.to_numpy().flatten().max(), 1))])
        res.ax_heatmap.set_ylabel('')
        xlabel = res.ax_heatmap.get_xmajorticklabels()
        ylabel = res.ax_heatmap.get_ymajorticklabels()
        hm = res.ax_heatmap.get_position()
        colden = res.ax_col_dendrogram.get_position()
        rowden = res.ax_row_dendrogram.get_position()

        res.ax_heatmap.set_xticklabels(xlabel, fontsize=xlab_siz, rotation=xlab_rot)
        res.ax_heatmap.set_yticklabels(ylabel, fontsize=ylab_siz)
        res.ax_heatmap.set_position([hm.x0 * xheat, hm.y0 * yheat, hm.width * heat_wf, hm.height * heat_hf])
        res.ax_row_dendrogram.set_position([rowden.x0 * rden_xf, rowden.y0 * rden_yf, rowden.width * rden_wf, rowden.height * rden_hf])
        res.ax_col_dendrogram.set_position([colden.x0 * cden_xf, colden.y0 * cden_yf, colden.width * cden_wf, colden.height * cden_hf])
        # plt.show()
        res.savefig(file_pth, dpi=300)
        plt.close('all')

    @staticmethod
    def plotbxplt(data, pth, ylabel):
        '''
        - plots boxplots in a row
        :param data: dataframe
        :param pth: path to file where to save the graph
        :param ylabel: Name in the vertical label
        '''
        sns.set_style('white')
        plt.figure(figsize=(21, 6.5))
        g = sns.boxplot(data=data, x='Subsystem', y=ylabel, palette='coolwarm')
        _, xlabels = plt.xticks()
        g.set_xticklabels(xlabels, size=16 * 1.45)
        _ = g.set_yticklabels(np.round(g.get_yticks(), 3), size=12 * 1.45)
        g.set_xlabel('')
        g.set_ylabel(g.get_ylabel(), fontsize=18 * 1.45)
        plt.tight_layout()
        plt.savefig(pth)
        plt.close()

    @staticmethod
    def dist_plt(frm, fld, fnm):
        '''
        - creates a histogram
        :param frm: dataframe
        :param fld: folder where to save
        :param fnm: name of file to save
        '''
        sns.set_style('white')
        plt.figure(figsize=(10, 6))
        plt.hist(frm.to_numpy().flatten())
        plt.tight_layout()
        pth = os.path.join(fld, fnm)
        plt.savefig(pth)
        plt.close()


