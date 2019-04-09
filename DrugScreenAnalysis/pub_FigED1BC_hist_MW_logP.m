%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED1BC_hist_MW_logP
% read DrugBank and PharmaCon drug characteristics tables
% plot histograms of molecular weight (MW) and logP distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% takes as input infile_drug_experiment
% infile_figED1_DrugBankDrugInfo = 'drugbank_approved_MW_150_1000_functional_groups_all.csv';
% outfile_figED1B_MW = ['Output' filesep 'FigE1_B_hist_drug_mol_weight.pdf'];
% outfile_figED1B_logP = ['Output' filesep 'FigE1_C_hist_drug_logP.pdf'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import data
ExperimentParameters = import_experiment_parameters(infile_drug_experiment);
drugMolWeight_selected = ExperimentParameters.DrugMolWeight;
drugLogP_selected = ExperimentParameters.DrugLogP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if DrugBank file is available, plot ditribution of drug bank drugs 
if is_available_DrugBankFile
    drugBankInfoTable = readtable(infile_figED1_DrugBankDrugInfo);
    % columnDrugName = 'GENERIC_NAME';
    % columnDrugMolWeight = 'MOLECULAR_WEIGHT';
    % columnDrugLogP = 'JCHEM_LOGP';
    drugMolWeight_drugbank = drugBankInfoTable.MOLECULAR_WEIGHT;
    drugLogP_drugbank = drugBankInfoTable.JCHEM_LOGP;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hist molecular weight
    figure
    hold on

    histogram(drugMolWeight_drugbank)

    histogram(drugMolWeight_selected)
    xlabel('Molecular weight')
    ylabel('Number of drugs')
    legend({'DrugBank approved', 'Selected 271'})
else
    figure
    hold on

    histogram(drugMolWeight_selected)
    xlabel('Molecular weight')
    ylabel('Number of drugs')
    legend({'Selected 271'})
end
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            outfile_figED1B_MW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hist logP
% if DrugBank file is available, plot ditribution of drug bank drugs 
if is_available_DrugBankFile
    figure
    histogram(drugLogP_drugbank)
    hold on
    histogram(drugLogP_selected)
    xlabel('LogP')
    ylabel('Number of drugs')
    legend({'DrugBank approved', 'Selected 271 drugs'})
else
    figure
    histogram(drugLogP_selected)
    xlabel('LogP')
    ylabel('Number of drugs')
    legend({'Selected 271 drugs'})
end
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            outfile_figED1B_logP)

