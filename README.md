This repository contains the description of the steps and original python scripts used to do the analyses presented in the manuscript …...
All data (used in both original Python and adapted Matlab scripts) together with the adapted Matlab scripts are deposited at ... zenodo?

Python scripts were run on Linux. The MATLAB scripts were adapted from [Human-GEM GitHub repository](https://github.com/SysBioChalmers/Human-GEM) (version 1.3.0) and [here](https://doi.org/10.5281/zenodo.3577466), and were run on Windows, except the one run on a server (***generate_human_ecModels_NCI60_batch.m***).


### Required Python modules:
* [mewpy v0.1.12](https://github.com/BioSystemsUM/MEWpy)
* [troppo](https://github.com/BioSystemsUM/troppo)
* [cobamp v0.1.4](https://github.com/BioSystemsUM/cobamp)
* pandas v1.3.5
* cobrapy v0.25.0
* numpy v1.21.5
* matplotlib v3.2.2
* scipy v1.5.2
* seaborn v0.11.0
* json5 v0.9.5
* regex 2020.10.15
* faker v5.8.0
* scipy v1.5.2
* scikit-learn v0.23.2
* multiprocess v0.70.11.1
* pathos v0.2.7


### Software used to run code on MATLAB:
* MATLAB R2021b
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN), commit [f80bdbd](https://github.com/SysBioChalmers/RAVEN/commit/f80bdbd2f888731f9d45e8941ed26197605c4ce5)
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox) (version 3.0.6)
* The [Human-GEM GitHub repository](https://github.com/SysBioChalmers/Human-GEM) (version 1.3.0)
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for MATLAB (version 9.5.1)
* [GECKO repository](https://github.com/SysBioChalmers/GECKO) (version 1.3.5)
* [libSBML MATLAB API](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/) ](version 5.13.0)

### Steps to reproduce the analysis:
- replace ***home_path*** by your home directory in every python script.
#### A. Reconstruct models with FASTCORE using python pipeline based on [Richelle's article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006867)
1. run ***richelle_pipe.py***
    - creates a generic DNA methylation model: adds reactions, gene rules and metabolites related with DNA (de)methylation to ***Human1v12*** model
    - removes blocked reactions, and checks whether cell biomass is being produced and that (de)methylation reactions are not blocked
    - in generic model, checks which tasks should_fail=False can be done by generic model, and extract reactions necessary for those, both for "consensus" (tasks from Richelle et al. that may be done by some cell types but not others) and "essential" tasks (those done by all cell types). the "consensus" list includes extra DNA demethylation tasks.
    - converts transcriptomics data to gene scores using the best performing parameters discovered by [Richelle et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006867) for these cell lines: local threshold (lc_thr) = average of each gene, upper global threshold (ug_thr) = 75th percentile, lower global threshold (lg_thr) = 25th percentile. gene score = 5 * log(1 + (expression/threshold)), where threshold = ug_thr if lc_thr >= ug_thr, threshold = lg_thr if lc_thr <= lg_thr, otherwise threshold = lc_thr.
    - obtains reaction scores using GPR rules. ('AND', 'OR') = (min, max)
    - protects (gives the highest score) to reactions that are necessary for all "essential" tasks
    - gives results for with and without protection of the "consensus" tasks (cell line specific tasks).
    - runs FASTCORE, where core reactions have score > 5 * log(2)  
    - removes reactions FASTCORE indicates to remove but keeps exchange reactions and uncatalyzed DNA demethylation reactions that (because are neither essential for tasks, nor have an associated reaction score) would be always automatically excluded.
    - does parsimonious FBA with biomass as objective and models closed (uptake of only medium components)
    - checks whether models are feasible, produce biomass and methylate DNA
2. run scripts to create GECKO models in a cluster, adapted from [Human-GEM GitHub repository](https://github.com/SysBioChalmers/Human-GEM) (version 1.3.0) – [original data files here](https://doi.org/10.5281/zenodo.3577466).
* clone GECKO repository inside folder ***Human1_Publication_Data_Scripts\ec_GEMs\ComplementaryScripts***
* transfer files from ***epigen/support/models_richelle_pipe/fastcore*** or from ***epigen/support/models_richelle_pipe/fastcore/including_tsks*** to folder ***Human1_Publication_Data_Scripts/ec_GEMs/models/humanGEM_cellLines*** in a cluster,
to get GECKO models with or without cell-type specific tasks.
* remember to make folders executable
* adapt ***main.sh*** for number of jobs (models) to run and run ***main.sh*** from folder ***Human1_Publication_Data_Scripts/ec_GEMs/ComplementaryScripts***
* ***main.sh*** creates jobs (one for each model) by calling the script ***generate_human_ecModels_NCI60_batch.sh***, which in turn calls ***generate_human_ecModels_NCI60_batch.m***

#### B. Reconstruct models with tINIT using pipeline of [Robinson's article](https://www.science.org/doi/10.1126/scisignal.aaz1482)
1. run ***prepare_gen_md_for_matlab.py***
   * saves model as ***.mat*** object (***epigen/support/models/prodDNAtot.mat***)
   * saves transcriptomics data as ***.mat*** object (***epigen/data/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.mat***)
2. run scripts adapted from [Human-GEM GitHub repository](https://github.com/SysBioChalmers/Human-GEM) (version 1.3.0) – [original data files here](https://doi.org/10.5281/zenodo.3577466).
* clone [GECKO repository (v1.3.5)](https://github.com/SysBioChalmers/GECKO/archive/refs/tags/v1.3.5.zip) inside directory ***Human1_Publication_Data_Scripts\ec_GEMs\ComplementaryScripts***
* transfer files ***epigen/support/transcriptomics/CCLE_RNAseq_rsem_genes_tpm_20180929.mat***  and ***epigen/support/models/prodDNAtot.mat*** to ***Human1_Publication_Data_Scripts\tINIT_GEMs\data***.
* copy file ***data/tasks/metabolicTasks_Essential.xlsx*** to folder ***Human1_Publication_Data_Scripts\tINIT_GEMs\metabolic_tasks***, make sure the excel sheet name is ***TASKS***.
* replace ***x*** by ***b*** in the function ***Human-GEM\ComplementaryScripts\Functions\addBoundaryMets.m***
*  add ***gen_all_tINIT_models_two.m*** (adapted from the script: gen_all_tINIT_models.m) to folder ***Human1_Publication_Data_Scripts\tINIT_GEMs*** and run the script
    * final models did all generic metabolic tasks (done by all cell lines)
3. run ***inbetween.py***
* transfer files from folder ***Human1_Publication_Data_Scripts\tINIT_GEMs\run_tINIT_outputs*** to ***epigen/data/models_tINIT_human_pipe/init***
* run the script ***in_between.py***
    - gives results for with and without protection of the "consensus" tasks (cell line specific tasks).
    - always adds essential reactions for the DNA demethylation tasks if those tasks are suppose to be done by the cell line, even when testing without "consensus" tasks
    - adds non-catalyzed reactions involved in DNA demethylation that are not essential for demethylation tasks
    - does parsimonious FBA with biomass as objective and models closed (uptake of only medium components)
    - checks whether models are feasible, produce biomass and methylate DNA
4. run scripts to create GECKO models inside a cluster, adapted from [Human-GEM GitHub repository](https://github.com/SysBioChalmers/Human-GEM) (version 1.3.0) – [original data files here](https://doi.org/10.5281/zenodo.3577466).
* transfer files from ***epigen/support/models_tINIT_human_pipe/init/including_tsks*** or from ***epigen/support/models_tINIT_human_pipe/init*** or from ***epigen/support/models_tINIT_human_pipe/init/notsk_wdemethtsk*** to folder ***Human1_Publication_Data_Scripts/ec_GEMs/models/humanGEM_cellLines*** in a cluster, to get GECKO models with or without cell-type specific tasks or without cell-specific tasks except DNA demethylation ones.
* adapt ***main.sh*** for number of jobs (models) to run and run ***main.sh*** from folder ***Human1_Publication_Data_Scripts/ec_GEMs/ComplementaryScripts***
* ***main.sh*** creates jobs (one for each model) by calling the script ***generate_human_ecModels_NCI60_batch.sh***, which in turn calls ***generate_human_ecModels_NCI60_batch.m***

#### C. Obtain gecko model of generic traditional model
* move traditional generic GSMM file ***epigen/support/models/prodDNAtot.mat*** to folder ***Human1_Publication_Data_Scripts/ec_GEMs/models/humanGEM_cellLines*** in a cluster, and run script with one job

#### D. Simulations with traditional models created with Richelle's and Robinson's pipelines
* run the script ***GEM_simul.py***:
  - creates scatter plots with log10 of abs. val. of predicted fluxes vs measured fluxes of exchange reactions of 26 metabolites.
  - creates histograms with the distribution of absolute values of measured and simulated fluxes before and after logarithmization.
  - creates scatter plots with log10 val. of predicted vs measured growth rates
  - creates boxplots with relative errors of predicted growth rates.

#### E. Simulations with GECKO models created with Richelle's and Robinson's pipelines
* transfer folders with cell line names from ***Human1_Publication_Data_Scripts/ec_GEMs/models/*** to folder ***epigen/support/ecGEMs_richelle/fastcore/including_tsks*** or ***epigen/support/ecGEMs_richelle/fastcore/no_tsks*** or ***epigen/support/ecGEMs_human1/init/including_tsks*** or ***epigen/support/ecGEMs_human1/init/no_tsks*** or ***epigen/support/ecGEMs_human1/init/notsk_wdemethtsk***, depending on what cell-specific tasks models do or not and depending on the type of reconstruction pipeline applied
* run the script ***ecGEM_simul.py***:
   - if required, it allows the replacement of 'prodDNAtot' reaction by an equivalent one reflecting the cell line-specific ratio of DNA methylation
   - creates scatter plots with log10 of abs. val. of predicted fluxes vs measured fluxes of exchange reactions of 26 metabolites.
   - creates histograms with the distribution of abs. values of measured and simulated fluxes before and after logarithmization.
   - creates scatter plots with log10 of abs. val. of predicted vs measured growth rates
   - creates boxplots with relative errors of predicted growth rates.
   - creates scatter plots with correlation level between simulated flux of different reactions and DNA methylation levels in different genomic positions (enhancers, TSS, ...).
   - creates a scatter plot with correlation level between measured biomass and measured methylation level of region upstream of TSS.
* run script ***pathways.py***:
   - creates boxplots with average simulated flux and protein usage for top 5 pathways and central carbon + DNA (de)/methylation pathways across all cell lines/models
* run script ***pthw_target_corr.py***:
   - produces tables with pathways whose average flux and protein usage significantly correlated (p-value < 0.05) with overall DNA methylation levels independently of cell growth rate across the cell lines.
   - produces tables with individual reactions/enzymes whose flux and protein usage significantly correlated (p-value < 0.05) with overall DNA methylation levels independently of cell growth rate across the cell lines.
   - creates tables with pathways of reactions corresponding to genes whose transcription was significantly correlated with the ratio between overall DNA methylation level and experimental cell growth rate.
   - checks whether certain genes are in the list of genes whose transcription was significantly correlated with the ratio between overall DNA methylation level and experimental cell growth rate.
   - intersect:
     - genes whose promoter methylation significantly correlated with its transcription.
     - genes whose promoter methylation significantly correlated with the cell growth rate across the different cell lines.
     - genes associated with reactions whose flux, or the genes associated with enzymes whose protein usage, significantly correlated with the cell growth rate.
