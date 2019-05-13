%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIPELINE 1: DRUG SCREEN ANALYSIS
% In this pipeline, the drug-bacteria screen data is analyzed
% and visualized as clustergram, barplots and heatmaps
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
% Main drug screen analysis workflow
pub_workflow_analyze_drug_screen
% calculate drug fold changes for each species 
% between t=12g and t=0h and to control
% save drug changes into a table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure 1B
pub_Fig1B_plot_PCA
% provided the PCA coordinates for the selected drugs and 
% DrugBANK drugs plot PCA plot characterizing the chemical space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure 1C and TableS3 in clustergram order
pub_Fig1C_plot_clustergram
% provided the output of pub_workflow_analyze_drug_screen (drug fold changes)
% make a clustergram and save drugs fold changes in clustergram order to a table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figures E1 B,C 
pub_FigED1BC_hist_MW_logP
% plot histograms of molecular weight and logP values for selected 
% drugs and drugs from DrugBank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure E1 D
pub_FigED1D_hist_intestinal_concentrations
% plot hostogram of estimated colon concentrations for selected drugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure E1 E
pub_FigED1E_drugs_number_per_threshold
% plot number of significantly changing drugs depending on the 
% consumption threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure E1 F,G
pub_FigED1FG_bar_drug_bacteria
% plot barplots of number of drugs metabolized by each species 
% and number of species metabolizing each drug
% at three thresholds: 20%, 50% and 80% drug consumption after 12h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Extended Data figure 2A
pub_FigED2A_bar_drugs_FC
% given the table of drug fold chanhges, 
% bar plot drug consumption % after 12h of incubation with each bacterium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce figure E2 B and Table S4
pub_FigED2B_drugs_functional_groups_and_enrichments
% provided drug functional group information from Table S2 
% metabolized drugs from Table S3 in clustergram order
% and manually defined sets of drugs
% (Cluster1 and Cluster2 from Figure 1C and Figure E3)
% NOTE: apart from clustergram-metabolized drugs, 
% drug sets from the CLuster I and II are manually defined inside the script 
% pub_FigED2B_drugs_functional_groups_and_enrichments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
