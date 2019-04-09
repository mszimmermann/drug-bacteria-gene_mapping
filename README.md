# drug-bacteria-gene_mapping
Mapping human microbiome drug metabolism by gut bacteria and their genes

Data analysis and visualization pipelines for screening the capacity of common gut bacteria to metabolize oral drugs.
The project consists of four pipelines:

#####
Analysis schemes are provided in AnalysisSchemes folder.
Global variables and input/output file names for all pipelines are defined in Drug_bacteria_gene_mapping_variables.m.
#####

##### PIPELINE 1 ###########################################
1) Analyze drug screen (scheme is provided in pipeline_drug_screen_analysis.pdf)
Pipeline steps with comments are provided in pipeline_drug_screen_analysis.m)

Perform differential analysis and hierarchical clustering of the metabolism patterns of 271 oral drugs by 76 diverse human gut bacteria (and four pH controls). 
Perform selected drug characterization in terms of chemical properties and bacterial metabolism patterns. 
Perform chemical structure and functional group analysis of the drugs metabolized in a similar way across the tested bacteria. 
###### PIPELINE 1 FILES ####################################
Step-wise pipeline file: pipeline_drug_screen_analysis.m
Main analysis file: pub_workflow_analyze_drug_screen.m

Scripts producing the figures and the tables:
pub_Fig1B_plot_PCA.m
pub_Fig1C_plot_clustergram.m
pub_FigED1BC_hist_MW_logP.m
pub_FigED1D_hist_intestinal_concentrations.m
pub_FigEDE_drugs_number_per_threshold.m
pub_FigED1FG_bar_drug_bacteria.m
pub_FigED2A_bar_drugs_FC.m
pub_FigED2B_drugs_functional_groups_and_enrichments.m
############################################################


##### PIPELINE 2 ###########################################
2) Analyze drug metabolism products (scheme provided in pipeline_drug_metabolism_products_analysis.pdf)
Pipeline steps with comments are provided in pipeline_drug_metabolism_products_analysis.m)

Perform differential analysis of metabolites co-occuring with specific drugs measured by untargeted mass-spectrometry. 
Perform mass-delta analysis of the drug-metabolite pairs. 
Perform functional group enrichment analysis of the common mass-deltas between drug and metabolites.
############################################################
###### PIPELINE 2 FILES ####################################
Step-wise pipeline file: pipeline_drug_metabolism_products_analysis.m
Main analysis file: pub_workflow_analyze_drug_metabolism_products.m

Scripts producing the figures and the tables:
pub_fig2a_volcano_plot.m
pub_fig2B_drug_metabolite_histogram.m
pub_fig2C_MZdelta_histogram.m
pub_fig2D_MZdelta_functional_enrichment.m
pub_fig2F_dexamethasone_single_species_barplot.m
pub_fig3A_btheta_specific_mzdeltas.m
############################################################


##### PIPELINE 3 ###########################################
3) Analyze drug-metabolizing bacterial gene products (scheme provided in pipeline_drug_metabolizing_gene_products_analysis.pdf)
Pipeline steps with comments are provided in pipeline_drug_metabolizing_gene_products_analysis.m)

Perform protein identity analysis for the identified bacterial gene products. 
Perform gene product enrichment analysis for drug-metabolizing species.
Perform correlation analysis for single species and identified gene products identity and drug metabolizing capacity in human gut communities. 
############################################################
###### PIPELINE 3 FILES ####################################
Step-wise pipeline file: pipeline_drug_metabolizing_gene_products_analysis.m
Main analysis files for community data:
pub_analyze_community_drug_metabolism.m
pub_Fig5ij_correlate_community_metabolism_with_features.m

Scripts producing the figures and the tables:
pub_Fig2h_community_correlation_scatter.m
pub_Fig5b_scatter_drugmetabolism_genePident.m
pub_Fig5d_drug_gene_enrichment_analysis.m
pub_Fig5fgFigE8d_bar_drugmet_genesPident.m
pub_Fig5h_community_correlation_scatter.m
pub_Fig5hij_community_drug_profiles.m 


##### PIPELINE 4 ###########################################
4) Reproduce other manuscript figures
Pipeline steps with comments are provided in pipeline_other_manuscript_figures.m)

############################################################
###### PIPELINE 4 FILES ####################################
Step-wise pipeline file: pipeline_other_manuscript_figures.m

Scripts producing the figures:
pub_Fig3C_plot_screening_bars.m
pub_Fig3D_plot_screening_heatmaps.m



############################################################
#####
The code is developed with MatLab2017b and Python3.6 and is distributed under the terms of the GNU General Public License (please read copyright_and_license and LICENSE files for details.
############################################################
#####
This project is part of the work by Michael Zimmermann, Maria Zimmermann-Kogadeeva, Rebekka Wegmann and Andrew L. Goodman.
