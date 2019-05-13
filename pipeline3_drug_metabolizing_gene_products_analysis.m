%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIPELINE 3: DRUG METABOLIZING GENE PRODUCTS ANALYSIS
% In this pipeline, the drug-bacteria screen data is analyzed
% together with identified bacterial gene products responsible for drug 
% metabolism for associations
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag whether to clear variables after calling each script
clear_var_flag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Figures 5b,d,e and Extended Data Figure 10e 
pub_Fig5bdeFigE10e_bar_drugmet_genesPident
% provided the drug-bacteria clustergram, drug-gene table and 
% BLAST results for genes and single species genomes, 
% plot barplots of single species drug metabolism and protein identity
% at given identity threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Figures 5b and Extended data Figure 10c -
% scatterplots between drug metabolism and protein identity
pub_Fig5bED10c_scatter_drugmetabolism_genePident
% provided the drug-bacteria clustergram, drug-gene table and 
% BLAST results for genes and single species genomes, 
% plot scatter plots of single species drug metabolism and protein identity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Figures 5d - drug-gene enrichment analysis
pub_Fig5c_drug_gene_enrichment_analysis
% provided the drug-bacteria clustergram, drug-gene table and 
% BLAST results for genes and single species genomes, 
% perform enrichment analysis of identity to specific genes
% among species metabolizing specific drugs
% and save results to file and table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if clear_var_flag
    clearvars -except clear_var_flag
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMUNITY DRUG METABOLISM ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% provided metabolic measurements of drugs and metabolites, 
% calculate drug and metabolite conversion clopes in human communities
pub_analyze_community_drug_metabolism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Figures 5g,h and Extended data figures 11 d,e,f -
% correlation bars between drug metabolism and
% community features (CFU, OTU, phylum, genus, protein abundances)
pub_Fig5ghED11def_correlate_community_metabolism_with_features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Extended data Figure 5f 
% - scatter plot of dexamethasone metabolism and community features
pub_FigED5f_community_correlation_scatter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if clear_var_flag
    clearvars -except clear_var_flag
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Figures 5f,g,h - drug and metabolite profiles in communities
pub_Fig5fgh_community_drug_profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if clear_var_flag
    clearvars -except clear_var_flag
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reproduce Figures 5f and Extended Data Figure 11a,b,e,f
% and calculate correlation between community
% drug metabolism and features
pub_Fig5f_FigE11abef_community_correlation_scatter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
