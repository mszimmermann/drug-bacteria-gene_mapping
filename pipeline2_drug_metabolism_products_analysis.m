%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIPELINE 2: DRUG METABOLISM PRODUCTS ANALYSIS
% In this pipeline, the drug-bacteria screen untargeted metabolomics data 
% is analyzed, for each drug candidate metabolites are identified based on
% their occurance in the drug pools;
% candidate drug metabolites are filtered 
% mass deltas between drugs and candidate metabolites are calculated
% and chemical functional group enrichment is performed on the 
% most common mass deltas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline can be run from the main folder (Drug_bacteria_genes_mapping).
% Input and output file names and global variables are defined in the file
% Drug_bacteria_gene_mapping_variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change directory to Drug_bacteria_genes_mapping
% and add all subfolders
addpath(genpath(['.' filesep]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag whether to clear variables after calling each script
clear_var_flag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main drug candidate metabolites analysis workflow
pub_workflow_analyze_drug_metabolism_products
% calculate intensity fold changes of metabolites between drug-containing 
% pools and all other pools
% OPTIONAL TIME CONSUMING STEP:
% calculate automatic p(FDR) threshold for pooling based on positive
% controls (drugs found in pools) and negative controls (drug pooling
% scheme) - use flag calculate_pThresholdPoolingFlag = 0; to control
% whether re-calculate the p(FDR) threshold or use the default (10^(-6))
% merge and filter candidate metabolites for all drugs
% calculate fold changes for merged metabolites between 12h and 0h
% save metabolite filtering, mass deltas and fold change information to tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure 2a - volcano plot of metabolites in diltiazem pools in 
% B.thetaiotaomicron incubation
pub_fig2a_volcano_plot
% produces volcano plot provided metabolite data, targeted drug and species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure 2B histogram of metabolites per drug
pub_fig2B_drug_metabolite_histogram
% provided data in the filtered metabolite table S5
% created in pub_workflow_analyze_drug_metabolism_products
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure 2C histogram of mass deltas between metabolites and drugs
pub_fig2C_MZdelta_histogram
% provided data in the filtered metabolite table S5
% created in pub_workflow_analyze_drug_metabolism_products
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure 2d and table S4 of enrichment of chemical functional groups
% among drugs having the same mass differents to te corresponding metabolites
pub_fig2D_MZdelta_functional_enrichment
% provided the drug information and the filtered metabolite table S5
% created in pub_workflow_analyze_drug_metabolism_products
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Extended Data Figure 4b - barplot of dexamethasone metabolism by single species
pub_figED4b_dexamethasone_single_species_barplot
% provided metabolite and drug data created in pipeline 1 and pub_workflow_analyze_drug_metabolism_products
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Extended Data figure 6a - B.thetaiotaomicron specific mass deltas 
% between drugs and metabolites
pub_figED6A_btheta_specific_mzdeltas
% provided drug clustegram created in pipeline 1 and 
% metabolite data created in pub_workflow_analyze_drug_metabolism_products
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
