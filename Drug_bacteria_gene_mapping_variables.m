%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE INITIALIZATION: Drug_bacteria_gene_mapping_variables
% This file contains global variables and 
% input/output file names for pipelines 
% in Drug-gene-metabolite_mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if function "contains" exist, 
% If not, rename local function "mycontains" to "contains" to be used
% throughout the scipt
if ~exist('contains','file')
    mycontains = which('mycontains');
    mycontains_new = strrep(mycontains, [filesep 'mycontains.m'], [filesep 'contains.m']);
    movefile(mycontains, mycontains_new);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define data folders and check their existence
DataFolder = 'Data';
DataInput = [DataFolder filesep 'Input'];
DataOutput = [DataFolder filesep 'Output'];
DataPipeline1 = [DataFolder filesep 'Raw_DrugScreenAnalysis'];
DataPipeline2 = [DataFolder filesep 'Raw_DrugMetabolismProducts'];
DataPipeline3 = [DataFolder filesep 'Raw_DrugMetabolizingGeneProducts'];
DataPipeline4 = [DataFolder filesep 'Raw_GainOfFunction'];

testFolders = {DataFolder ...
   DataInput...
   DataOutput...
   DataPipeline1...
   DataPipeline2...
   DataPipeline3...
   DataPipeline4};

testFolderResults = zeros(length(testFolders),1);
for i=1:length(testFolders)
    testFolderResults(i) = ~exist(testFolders{i}, 'dir');
end
missingFolders = strjoin(testFolders(testFolderResults==1), '\n');

if nnz(testFolderResults)
    error_msg = ['Error: The following Data folders not found:\n'...
                 '%s\n\n'...
                 'Please download data files from FigShare:\n'...
                 'https://doi.org/10.6084/m9.figshare.8119058.v1\n'...
                 'and move the Data folder to the master directory.\n\n'...
                 'Expected Data directory tree:\n'...
                 '%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n'...
                 'If you have a different directory structure, please change the variables:\n'...
                 'DataFolder, DataInput, DataOutput, DataPipeline1, DataPipeline2, DataPipeline3, DataPipeline4, '...
                 'accordingly.\n\n'];
    error(error_msg, missingFolders,...
                    [pwd filesep DataFolder],...
                    [pwd filesep DataInput],...
                    [pwd filesep DataOutput],...
                    [pwd filesep DataPipeline1],...
                    [pwd filesep DataPipeline2],...
                    [pwd filesep DataPipeline3],...
                    [pwd filesep DataPipeline4])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL VARIABLES (THRESHOLDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pThreshold = 0.05; % for drug and metabolite fold change significance
intensityNoise = 5000; %background intensity (lowest threshold)
fcThresholdDrug = abs(log2(0.8)); % 20% reduction of drug over time
drug_outlierSTDcutoff = 5; % cutoffs to detect one-out-of-four outliers in the drug pools
drug_outlierMEANcutoff = 0.2; % cutoffs to detect one-out-of-four outliers in the drug pools
drug_fastMetabolizerThreshold = -5; % threshold for FC to control at t=0 for fast metabolizers
massThreshold = 0.002; % Mass threshold for identifying metabolites from untargeted metabolomics
RTthreshold = 0.15; % Retention time threshold for identifying metabolites from untargeted metabolomics
pThresholdPooling = 10^(-6); % significance threshold for metabolites detected in the drug pool vs other pools
fcThresholdPool = 1; % fold change threshold between metabolite in drug pool vs other pools
calculate_pThresholdPoolingFlag = 0; % whether to calculate significance threshold based on drug and random pooling scheme metabolite
                                     % distributions (time consuming step)
fcThresholdMetT12T0Pool = 1; % fold change threshold for metabolite between t12 and t0
smoothMZdeltas = 0.002; % smppthing window to smooth mass deltas between metabolites and drugs 
is_available_DrugBankFile = 0; % DrugBank file is not part of the pipeline and has to be downloaded and prepared separately
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT-SPECIFIC INPUT-OUTPUT FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables for PIPELINE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_workflow_analyze_drug_screen
infile_single_species = [DataInput filesep 'TableS1_AssayedSingleSpecies.csv'];
infile_drug_experiment = [DataInput filesep 'TableS2_DrugScreenInfo.csv'];
infolder_drug_data = [DataPipeline1 filesep];
outfile_Table3_drug_fold_changes = [DataOutput filesep 'TableS3_drug_fold_changes.csv'];
% save individual drug fold changes for bar plots
outfile_Table_individual_drug_fold_changes = [DataOutput filesep 'Table_individual_drug_fold_changes.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_Fig1B_plot_PCA
% takes as input infile_drug_experiment
infile_fig1B_drugbankPCA = 'drugbank_pca.txt';
infile_fig1B_pharmacomPCA = 'pharmacon_pca.txt';
outfile_fig1B_PCA = [DataOutput filesep 'Fig1B_plot_rdkit_drugbank_selected_drugs_PCA.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_Fig1C_plot_clustergram
% takes as input outfile_Table3_drug_fold_changes
outfile_fig1c_clustergram = [DataOutput filesep 'Fig1C_clustergram_drug_percentChange.pdf'];
outfile_table_FCclustergram = [DataOutput filesep 'TableS3_drug_fold_changes_in_clustergram_order.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_FigED1BC_hist_MW_logP
% takes as input infile_drug_experiment
% DrugBank file is not part of the pipeline and has to be downloaded and prepared separately
% infile_figED1_DrugBankDrugInfo = 'drugbank_approved_MW_150_1000_functional_groups_all.csv';
outfile_figED1B_MW = [DataOutput filesep 'FigE1B_hist_drug_mol_weight.pdf'];
outfile_figED1B_logP = [DataOutput filesep 'FigE1C_hist_drug_logP.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_FigED1C_hist_intestinal_concentrations
% takes as input infile_drug_experiment
outfile_figED1D_intconc = [DataOutput filesep 'FigE1D_hist_predicted_intestinal_concentrations.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_FigED1FG_bar_drug_bacteria
% takes as input outfile_Table3_drug_fold_changes
outfile_figED1F_bar_bacteria_per_drug = [DataOutput filesep 'FigE1F_bar_degrading_species_per_drug.pdf'];
outfile_figED1F_bar_bacteria_per_drug_zoomed = [DataOutput filesep 'FigE1F_bar_degrading_species_per_drug_only_degraded_drugs.pdf'];
outfile_figED1G_bar_drug_per_bacteria = [DataOutput filesep 'FigE1G_bar_degraded_drugs_per_species.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_FigEDE_drugs_number_per_threshold
% takes as input outfile_Table3_drug_fold_changes
outfile_figED1E_drug_number_per_threshold = [DataOutput filesep 'FigE1E_drug_changes_vs_PercentFCthreshold.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_FigED2A_bar_drugs_FC
% takes as input infile_drug_experiment
% takes as input outfile_Table3_drug_fold_changes
% takes as input outfile_table_FCclustergram
outfile_figED2A_bar_drugFC = [DataOutput filesep 'FigE2A_bar_drugFC_in_clustergram_order.ps'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_FigED2B_drugs_functional_groups_and_enrichments
% takes as input infile_drug_experiment
% takes as input outfile_table_FCclustergram
% DrugBank file is not part of the pipeline and has to be downloaded and prepared separately
% takes as input infile_figED1_DrugBankDrugInfo
outfile_tableS4_functional_group_enrichments = [DataOutput filesep 'TableS4_functional_group_enrichment_changing_drugs.csv'];
outfile_figED2B_drug_functional_groups = [DataOutput filesep 'FigE2B_bar_functional_groups_metabolized_drugs.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables for PIPELINE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_workflow_analyze_drug_metabolism_products
% takes as input infile_drug_experiment
% takes as input infile_single_species
% takes as input outfile_table_FCclustergram
infolder_drug_metabolite_data = [DataPipeline2 filesep];
outfile_tableS6_metaboliteIntensityFCData = [DataOutput filesep 'TableS6_metabolite_intensity_and_fold_changes.csv'];
outfile_tableS5_metabolite_filtering = [DataOutput filesep 'TableS5_drug_metabolite_candidates_filtering.csv'];
outfile_figure_pFDR_pooling_cutoff = [DataOutput filesep 'fig_automatic_pFDRthreshold_pooling_metabolites.pdf'];
% save single replicate intensities for dexamethasopne metabolite for
% barplot Extended Data Figure 4b
outfile_table_dexamethasone_metabolite_intensityT12 = [DataOutput filesep 'Table_dexamethasone_metabolite_int_t12t0.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_fig2a_volcano_plot
% takes as input infile_drug_experiment
% takes as input infile_single_species
% takes as input infolder_drug_metabolite_data
fig2a_volcano_poolingDrug = 'DILTIAZEM';
fig2a_volcano_poolingSpecies = 'Bacteroides thetaiotaomicron';
outfile_fig2a_volcano = [DataOutput filesep 'Fig2A_volcano_diltiazem.ps'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_fig2B_drug_metabolite_histogram
% % takes as input outfile_table_FCclustergram
% % takes as input outfile_tableS5_metabolite_filtering
outfile_fig2b_metabolites_per_drug = [DataOutput filesep 'Fig2B_histogram_metabolite_per_drugs.pdf'];
outfile_fig2b_metabolites_per_drug_zoomed = [DataOutput filesep 'Fig2B_histogram_metabolite_per_drugs_zoomed.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_fig2C_MZdelta_histogram
% % takes as input outfile_tableS5_metabolite_filtering
outfile_fig2c_metabolites_drug_mass_deltas = [DataOutput filesep 'Fig2C_histogram_metabolite_drug_MassDeltas.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_fig2D_MZdelta_functional_enrichment
% % takes as input infile_drug_experiment
% % takes as input outfile_tableS5_metabolite_filtering
outfile_fig2d_clustergram_MZdelta_FGenrichment = [DataOutput filesep 'Fig2D_clustergram_massdelta_FG_enrichment.pdf'];
outfile_TableS4_mzdelta_chemical_functional_group_enrichment = [DataOutput filesep 'TableS4_part2_massdelta_FG_enrichment.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_figED4b_dexamethasone_single_species_barplot
% takes as input outfile_table_FCclustergram
% takes as input outfile_tableS6_metaboliteIntensityFCData
outfile_figED4b_dexamethasone_barplot = [DataOutput filesep 'FigED4b_dexamethasone_barplot.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_figED6A_btheta_specific_mzdeltas
% takes as input outfile_table_FCclustergram
% takes as input outfile_tableS5_metabolite_filtering
% takes as input outfile_tableS6_metaboliteIntensityFCData
outfile_figED6A_btheta_drug_metabolites = [DataOutput filesep 'FigED6A_scatter_bteta_drug_metabolites.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables for PIPELINE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED6b_plot_screening_bars
inputFolderDataScreen1 = [DataPipeline4 filesep 'DataScreen1' filesep];
outfile_FigED6b_barplot_diltiazem_GoF_screen1 = [DataOutput filesep 'FigED6b_barplot_diltiazem_GoF_screen1.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED6c_plot_screening_heatmaps
inputFolderDataScreen2 = [DataPipeline4 filesep 'DataScreen2' filesep];
outfile_FigED6c_heatmap_diltiazem_GoF_screen2 = [DataOutput filesep 'FigED6c_heatmap_diltiazem_GoF_screen2.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables for PIPELINE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_Fig5deFigE10e_bar_drugmet_genesPident
% takes as input outfile_table_FCclustergram
% takes as input infile_single_species
infileTable14BLAST = [DataInput filesep 'TableS14_Gene_Genome_BLAST_results.csv'];
infile_gene_drug_matrix = [DataInput filesep 'TableS13_gene_drug_matrix.csv'];
pidentThreshold = 50; % threshold of protein identity for plotting
outfile_Fig5deFigE10e_bar_drugmetabolism_genes_pident = [DataOutput filesep 'Fig5deFigE10e_bar_drugmetabolism_genes_pident.ps'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_Fig5bFigED10c_scatter_drugmetabolism_genePident
% takes as input infileTable14BLAST % BLAST table of identified genes and genomes of single species
% takes as input infile_gene_drug_matrix % Table of curated gene-drug interactions
% takes as input outfile_table_FCclustergram % metabolized drug fol changes in clustergram order
% takes as input infile_single_species % single species information
outfile_Fig5bFigED10c_scatter_drugmetabolism_genes_pident = [DataOutput filesep 'Fig5bFigED10c_scatter_drugmetabolism_genes_pident.ps'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_Fig5c_drugmet_gene_enrichment_analysis
% takes as input infileTable14BLAST % BLAST table of identified genes and genomes of single species
% takes as input infile_gene_drug_matrix % Table of curated gene-drug interactions
% takes as input outfile_table_FCclustergram % metabolized drug fol changes in clustergram order
% takes as input infile_single_species % single species information
outfile_fig5c_drug_gene_enrichment = [DataOutput filesep 'Fig5c_clustergram_gene_drug_enrichment.pdf'];
outfile_fig5c_combined_enrcihment = [DataOutput filesep 'Fig5c_clustergram_combined_gene_drug_enrichment.pdf'];
outfile_table15_gene_enrcihment_table = [DataOutput filesep 'TableS15_gene_drug_enrichment_analysis.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_analyze_community_drug_metabolism.m
infolder_DataCommunitiesDrugMetabolism = [DataPipeline3 filesep];
outfile_table_community_drug_metabolism = [DataOutput filesep 'TableS16_community_drug_metabolism.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_Fig5ij_correlate_community_metabolism_with_features.m
infile_community_CFU_qPCR = [DataInput filesep 'TableS9S16_community_CFU_and_qPCR.csv'];
infile_community_OTU = [DataInput filesep 'TableS17_community_OTUs.csv'];
infile_community_ShortBRED_proteins = [DataInput filesep 'TableS18_community_ShortBRED_proteins.csv'];
outfile_community_correlation_bars = [DataOutput filesep 'Fig5gh_FigE11_community_correlation_bars.ps'];
outfile_community_correlation_table = [DataOutput filesep 'TableS19_community_drug_metabolism_feature_correlation.csv'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_FigED5f_community_correlation_scatter
% takes as input infile_gene_drug_matrix 
% takes as input infile_community_CFU_qPCR 
% takes as input infile_community_OTU 
% takes as input infile_community_ShortBRED_proteins 
% takes as input outfile_table_community_drug_metabolism
outfile_FigED5f_dexamethasone_scatter = [DataOutput filesep 'FigED5f_dexamethasone_community_correlation_scatter.pdf'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pub_Fig5fgh_community_drug_profiles.m
% takes as input outfile_table_community_drug_metabolism
% mycolormap_communities to reproduce colormap in the figures 
mycolormap_communities = jet(64);
%permjet = fliplr(randperm(size(mycolormap_communities,1)));
permjet = [6,54,59,32,50,8,24,2,55,18,47,41,44,14,49,23,9,33,4,37,38,11,36,13,45,43,7,10,34,56,52,30,51,39,57,22,12,17,62,21,46,16,15,42,31,3,20,28,61,5,53,64,58,40,19,26,35,48,25,63,29,1,27,60];
mycolormap_communities = [repmat([0 0 0],4,1);...
                            mycolormap_communities(permjet,:)];
outfile_community_drug_met_profiles = [DataOutput filesep 'Fig5fgh_community_drug_met_profiles.ps'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline pub_Fig5f_FigE11_community_correlation_scatter
% takes as input infile_gene_drug_matrix 
% takes as input infile_community_CFU_qPCR 
% takes as input infile_community_OTU 
% takes as input infile_community_ShortBRED_proteins 
% takes as input outfile_table_community_drug_metabolism
outfile_Fig5f_community_correlation_scatter = [DataOutput filesep 'Fig5fFigED11_community_correlation_scatters.ps'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

