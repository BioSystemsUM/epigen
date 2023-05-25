home_path = '/home/tbarata'

import sys, os
content_root = 'epigen'
working_dir = os.path.join(home_path, content_root)
os.chdir(working_dir)
sys.path.extend([os.path.join(home_path, content_root, pkg) for pkg in ['code', 'code/src']])

from pathway_target_correlations import DataCorr

if '__name__' == '__main__':
    ### Pathways/individual genes/individual reactions which flux/protein usage correlates with overall methylation
    ### independently of growth rate:
    FLD_PTH='results/ecGEM_simul_human1/init/notsk_wdemethtsk/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    FLX_FL_NM='allflx_rc_substm.tsv'
    PRT_FL_NM='allprot_rc_substm.tsv'
    METHLT_FLD = 'data/methylation'
    ACHILLES_SMP_INFO = 'data/sample_info.csv'  # achilles dataset cell line info
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    ENSEMBL_SYMBOL_FILE = 'support/mappings/ensembl_symbol_human.json'
    div_bm=True # obtain results independently of growth rate
    # pathways which fluxes in average significantly correlates (directly/inversely) with DNA methylation level/growth rate:
    DataCorr.pthw_corr_meth_bmindp(fld_pth=FLD_PTH, fl_nm=FLX_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO, gen_gecko_pth=GEN_GECKO_PTH, div_bm=div_bm)
    # pathways which protein usage in average significantly correlates (directly/inversely) with DNA methylation level/growth rate:
    DataCorr.pthw_corr_meth_bmindp(fld_pth=FLD_PTH, fl_nm=PRT_FL_NM, methlt_dt_fld=METHLT_FLD, smp_info=ACHILLES_SMP_INFO, fl_nm_bm=FLX_FL_NM, gen_gecko_pth=GEN_GECKO_PTH, div_bm=div_bm, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE)

    ### Pathways over-represented in genes which gene expression correlates with overall DNA methylation:
    ### independently of growth rate:
    METHLT_FLD='data/methylation'
    TRSCP_PTH='data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.txt'
    ENVCOND='data/constraints/1218595databases1_corrected_further_cal.xls'
    GEN_GECKO_PTH='support/models/prodDNAtot/ecModel_batch.mat'
    ENSEMBL_SYMBOL_FILE='support/mappings/ensembl_symbol_human.json'
    ACHILLES_SMP_INFO='data/sample_info.csv'
    GEN_MD_PTH_SBML='support/models/prodDNAtot.xml'  # traditional generic model with DNA meth and demethylation reactions
    GEN_MD_PTH='data/models/Human-GEM.yml'
    FLD_PTH='results/ecGEM_simul_human1/init/notsk_wdemethtsk/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    DIV_GROWTH=True
    mth_trscp_corrvalue='direct'
    DataCorr.experimental_corr(methlt_dt_fld=METHLT_FLD, trscp_pth=TRSCP_PTH, mth_trscp_corrvalue=mth_trscp_corrvalue, envcond=ENVCOND,
                      gen_gecko_pth=GEN_GECKO_PTH, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, smp_info=ACHILLES_SMP_INFO,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, gen_md_pth=GEN_MD_PTH, fld_pth=FLD_PTH, div_growth=DIV_GROWTH)
    mth_trscp_corrvalue='inverse'
    DataCorr.experimental_corr(methlt_dt_fld=METHLT_FLD, trscp_pth=TRSCP_PTH, mth_trscp_corrvalue=mth_trscp_corrvalue, envcond=ENVCOND,
                      gen_gecko_pth=GEN_GECKO_PTH, ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, smp_info=ACHILLES_SMP_INFO,
                      gen_md_pth_sbml=GEN_MD_PTH_SBML, gen_md_pth=GEN_MD_PTH, fld_pth=FLD_PTH, div_growth=DIV_GROWTH)

    ### Pathways/individual genes/individual reactions which:
    ### - gene methylation correlates with gene expression
    ### - flux/protein usage correlates with growth rate
    ### - gene methylation correlates with growth rate
    ## with reaction flux:
    ENSEMBL_SYMBOL_FILE = 'support/mappings/ensembl_symbol_human.json'
    GEN_GECKO_PTH = 'support/models/prodDNAtot/ecModel_batch.mat'
    TRSCP_PTH = 'data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.txt'
    METH_FL = 'data/methylation/CCLE_RRBS_TSS1kb_20181022.txt'
    FLD_PTH = 'results/ecGEM_simul_human1/init/notsk_wdemethtsk/no_flux_constr/biomass_constr/cl_spec_DNAtot/corr_pth_mth'
    FLX_FL_NM = 'allflx_rc_substm.tsv'
    PRT_FL_NM = 'allprot_rc_substm.tsv'
    fl_nm = FLX_FL_NM
    DataCorr.mth_bm_trgt(ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, trscp_pth=TRSCP_PTH,
                meth_fl=METH_FL, fld_pth=FLD_PTH, fl_nm=fl_nm)
    # with protein usage:
    fl_nm = PRT_FL_NM
    fl_nm_bm=FLX_FL_NM
    DataCorr.mth_bm_trgt(ensembl_symbol_file=ENSEMBL_SYMBOL_FILE, gen_gecko_pth=GEN_GECKO_PTH, trscp_pth=TRSCP_PTH,
                meth_fl=METH_FL, fld_pth=FLD_PTH, fl_nm=fl_nm, fl_nm_bm=fl_nm_bm)

