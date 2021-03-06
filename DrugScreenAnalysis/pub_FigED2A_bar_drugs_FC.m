%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED2A_bar_drugs_FC
% read degraded drug fold change table and 
% plot drug FC per species in the clustergram order in a multipage ps file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% takes as input infile_drug_experiment
% takes as input outfile_Table3_drug_fold_changes
% takes as input outfile_table_FCclustergram
% takes as input outfile_Table_individual_drug_fold_changes
% outfile_figED2A_bar_drugFC = ['Output' filesep 'FigE2A_bar_drugFC_in_clustergram_order.ps'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import Experimental Parameters
ExperimentParameters = import_experiment_parameters(infile_drug_experiment);

% Make a list with abbreviated compound names
ExperimentParameters.DrugNamesAbbr = cellfun(@(x) x(1:6), ExperimentParameters.DrugNames, 'unif', 0);
%distingush betamethasone acetate and valerate
ExperimentParameters.DrugNamesAbbr{ismember(ExperimentParameters.DrugNames,...
                                            'BETAMETHASONE ACETATE')} = 'BETAMA';
ExperimentParameters.DrugNamesAbbr{ismember(ExperimentParameters.DrugNames,...
                                            'BETAMETHASONE VALERATE')} = 'BETAMV';
%distingush trimetazidine and trimethobenzamide
ExperimentParameters.DrugNamesAbbr{ismember(ExperimentParameters.DrugNames,...
                                            'TRIMETAZIDINE DIHYDROCHLORIDE')} = 'TRIMET';
ExperimentParameters.DrugNamesAbbr{ismember(ExperimentParameters.DrugNames,...
                                            'TRIMETHOBENZAMIDE')} = 'TRIMEH';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read drug FC information
drugAnalysisResults = readtable(outfile_Table3_drug_fold_changes,...
                                'HeaderLines',1); %skip first line
drugAnalysisColumnNames = drugAnalysisResults.Properties.VariableNames;                            
%convert table to cell array
drugAnalysisResults = table2cell(drugAnalysisResults);
FCcolumns = cellfun(@(x) contains(x, 'FC') &...
                         ~contains(x, 'FCSTD') &...
                         ~contains(x, 'DrugAdaptive'),...
                         drugAnalysisColumnNames);
FCSTDcolumns = cellfun(@(x) contains(x, 'FCSTD'),...
                         drugAnalysisColumnNames);
PFDRcolumns = cellfun(@(x) contains(x, 'p_FDR'),...
                         drugAnalysisColumnNames);                     
drugFC12to0_t0toCTRL_combined = cell2mat(drugAnalysisResults(:,FCcolumns))';
drugFC12to0_t0toCTRL_STDcombined = cell2mat(drugAnalysisResults(:,FCSTDcolumns))';
drugP12to0_t0toCTRL_combined = cell2mat(drugAnalysisResults(:,PFDRcolumns))';

DrugNames = drugAnalysisResults(:,ismember(drugAnalysisColumnNames, 'DrugName'));
drugFCadaptive = cell2mat(drugAnalysisResults(:,...
                          cellfun(@(x) contains(x, 'DrugAdaptiveFC'),...
                          drugAnalysisColumnNames)));
% get species names
totSpeciesNamesSingle = cellfun(@(x) x(strfind(x, 'FC')+2:end),...
                                     drugAnalysisColumnNames(FCcolumns), 'unif', 0)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read clustergram to plot in clustergram order
drugAnalysisClustergram = readtable(outfile_table_FCclustergram,...
                                'HeaderLines',1); %skip first line);
FCcolumns = cellfun(@(x) contains(x, 'FC') &...
                         ~contains(x, 'FCSTD') &...
                         ~contains(x, 'DrugAdaptive'),...
                         drugAnalysisClustergram.Properties.VariableNames);

clustergramSpecies = cellfun(@(x) x(strfind(x, 'FC')+2:end),...
                                  drugAnalysisClustergram.Properties.VariableNames(FCcolumns), 'unif', 0)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read individual data for dexamethasone metabolite in scindens
drugFoldChanges_individual = readtable(outfile_Table_individual_drug_fold_changes,...
                            'HeaderLines', 1);
drugAnalysisColumnNames_ind = drugFoldChanges_individual.Properties.VariableNames;
FCcolumns_ind = cellfun(@(x) contains(x, 'x_Consumed'),...
                         drugAnalysisColumnNames_ind);
drugFoldChanges_individual_data = drugFoldChanges_individual{:,FCcolumns_ind};
drugFC12to0_Species_ind = drugAnalysisColumnNames_ind(FCcolumns_ind);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedA', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedB', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedC', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = cellfun(@(x) strrep(x,'x_ConsumedD', ''),drugFC12to0_Species_ind, 'unif',0);
drugFC12to0_Species_ind = reshape(drugFC12to0_Species_ind,4,[])';

% intersect species names (individual drug FC)
[~, drug_orig_ind, drug_ind_idx] = intersect(lower(clustergramSpecies), lower(drugFC12to0_Species_ind(:,1)), 'stable');
% add control at the end
[~,idxdiff] = setdiff(lower(drugFC12to0_Species_ind(:,1)), lower(clustergramSpecies));
drug_ind_idx = [idxdiff; drug_ind_idx];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each drug, plot FC from 12 to 0 for the species
bwidth = 0.5;
[~, ~, idx] = intersect(clustergramSpecies, totSpeciesNamesSingle, 'stable');
% add control at the end
[~,idxdiff] = setdiff(totSpeciesNamesSingle, clustergramSpecies);
idx = [idxdiff; idx];
for j=1:length(ExperimentParameters.DrugNamesAbbr)
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    curFC = drugFC12to0_t0toCTRL_combined(:,j);
    curSTD = drugFC12to0_t0toCTRL_STDcombined(:,j);
    % translate in %drug
    curSTD = 100*curSTD.*(2.^curFC)*log(2);
    curFC = 100*(1-2.^curFC);
    
    % get individual dots
    % get the drugs and select current drug
    drugFoldChanges_individual_drugIDX = cellfun(@(x) contains(lower(x), lower(ExperimentParameters.DrugNames{j})),...
                                                drugFoldChanges_individual.DrugName);
    % select one drug and reshape accroding to species
    cur_individual_data = drugFoldChanges_individual_data(drugFoldChanges_individual_drugIDX,:);
    cur_individual_data = reshape(cur_individual_data,4,[])';
    cur_individual_data = cur_individual_data(drug_ind_idx,:);
    
    curP = drugP12to0_t0toCTRL_combined(:,j);
    curFC = curFC(idx);
    curSTD = curSTD(idx);
    curP = curP(idx);
    barFCns = curFC;
    barFCns(curP<=pThreshold & curFC>=drugFCadaptive(j)) = 0;
    bh = barh(barFCns, bwidth);
    bh(1).FaceColor = [.5 .5 .5];
    hold on
    barFCs = curFC;
    barFCs(curP>pThreshold | curFC<drugFCadaptive(j)) = 0;
    barh(barFCs, bwidth, 'r')
    % use horizontal errorbars if release > 2016a
    matlabrelease = version('-release');
    matlabyear = str2double(matlabrelease(1:4));
    if matlabyear > 2016 || isequal(matlabrelease, '2016b')
        errorbar(curFC,1:length(curFC), curSTD, 'horizontal', '.k');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add individual data points
    scatter(cur_individual_data(:),...
            repmat(1:size(cur_individual_data,1),1,4),...
            20,'b', 'filled')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    set(gca, 'YTick', 1:length(idx))
    set(gca, 'YTickLabel', strcat(totSpeciesNamesSingle(idx), '; ', num2str(curP,3)))
    set(gca, 'FontSize', 6)
    title(ExperimentParameters.DrugNames{j})
    xlim([min(-20, floor(min(curFC-abs(curSTD)))), max(100, ceil(max(curFC+curSTD)))])
    ylim([0 length(idx)+1])
    orient landscape

    print(fig, '-painters', '-dpsc', '-r600', '-bestfit', '-append',...
          outfile_figED2A_bar_drugFC)
    close(fig)     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


