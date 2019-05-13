%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_figED6A_btheta_specific_mzdeltas
% read data from metabolite table
% plot B.theta metabolites on top of other species drug metabolites
% for all drugs metabolized by B.theta
% save figure to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% takes as input outfile_table_FCclustergram
% takes as input outfile_tableS5_metabolite_filtering
% takes as input outfile_tableS6_metaboliteIntensityFCData
% outfile_figED6A_btheta_drug_metabolites
% massThreshold = 0.002;
% fcThresholdMetT12T0Pool = 1; (to select significnatly changing B.theta
% pThreshold = 0.05;            %metabolites between t=12h and t=0h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter only good and changing in at least one species
drugMetTable = readtable(outfile_tableS5_metabolite_filtering);
drugCustergram = readtable(outfile_table_FCclustergram,...
                            'HeaderLines', 1);
% trim spaces from drug names if there are any
drugMetTable.ParentDrug = cellfun(@(x) strtrim(x), drugMetTable.ParentDrug, 'unif', 0);
DrugNames = drugCustergram.DrugName;

drugAnalysisColumnNames = drugCustergram.Properties.VariableNames;
drugAnalysisResults = table2cell(drugCustergram);
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

drugFCadaptive = cell2mat(drugAnalysisResults(:,...
                          cellfun(@(x) contains(x, 'DrugAdaptiveFC'),...
                          drugAnalysisColumnNames)));
%convert to log2FC
drugFCadaptive = abs(log2((100-drugFCadaptive)/100))';

totSpeciesNamesSingle = cellfun(@(x) x(strfind(x, 'FC')+2:end),...
                                     drugAnalysisColumnNames(FCcolumns), 'unif', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metFoldChanges = readtable(outfile_tableS6_metaboliteIntensityFCData,...
                            'HeaderLines', 1);
% trim spaces from drug names if there are any
metFoldChanges.ParentDrug = cellfun(@(x) strtrim(x), metFoldChanges.ParentDrug, 'unif', 0);
FCcolumns = cellfun(@(x) contains(x, 'FCT_12hVsT_0h'),...
                         metFoldChanges.Properties.VariableNames);
pFDRcolumns = cellfun(@(x) contains(x, 'p_FDR_T_12hVsT_0h'),...
                         metFoldChanges.Properties.VariableNames);   
Int12columns = cellfun(@(x) contains(x, 'IntensityMeanT_12h'),...
                         metFoldChanges.Properties.VariableNames);
metFC12to0 = metFoldChanges{:,FCcolumns};
% check if metFC12to0 is numeric, might be text diet to formatting
if ~isnumeric(metFC12to0)
    metFC12to0 = cellfun(@(x) strrep(x,'#NAME?', 'NaN'), metFC12to0, 'unif',0);
    metFC12to0 = cellfun(@(x) str2double(x), metFC12to0);
end
metFDR12to0 = metFoldChanges{:,pFDRcolumns};
metInt12 = metFoldChanges{:,Int12columns};
metFCspecies = metFoldChanges.Properties.VariableNames(FCcolumns)';
% sort FC in the order of drugMetTable
idx = zeros(size(drugMetTable,1),1);
for i = 1:size(drugMetTable,1)
    idx(i) = find(ismember(metFoldChanges.ParentDrug, drugMetTable.ParentDrug{i}) &...
               (metFoldChanges.MZ == drugMetTable.MZ(i)) &...
               (metFoldChanges.RT == drugMetTable.RT(i)));
end
metFC12to0 = metFC12to0(idx,:);
metFDR12to0 = metFDR12to0(idx,:);
metInt12 = metInt12(idx,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load B.theta metabolizing drugs from targeted analysis
degradedDrugNamesBtheta_targeted = zeros(length(DrugNames),1);
species_idx = cellfun(@(x) contains(lower(x), 'theta'), totSpeciesNamesSingle);
for i=1:length(DrugNames)
    degradedDrugNamesBtheta_targeted(i) = ...
        sum(drugFC12to0_t0toCTRL_combined(species_idx,i)<=-drugFCadaptive(i) &...
            drugP12to0_t0toCTRL_combined(species_idx,i)<=0.05);
end

% remove artimisinine
selected_drugs_all = DrugNames(degradedDrugNamesBtheta_targeted==3);
selected_drugs_all(ismember(selected_drugs_all, 'ARTEMISININ')) = [];

selected_drugs_all = flipud(sort(selected_drugs_all));

metMZdifference_todrug = zeros(size(drugMetTable,1),1);
drug_coordinate_selected_drugs = zeros(size(drugMetTable,1),1);
for i=1:size(drugMetTable,1)
    curdrugidx = ismember(DrugNames,...
                          drugMetTable.ParentDrug(i));
    metMZdifference_todrug(i) = drugMetTable.DrugMassDelta(i);
    if nnz(ismember(selected_drugs_all, drugMetTable.ParentDrug(i)))
        drug_coordinate_selected_drugs(i) = find(ismember(selected_drugs_all, drugMetTable.ParentDrug(i)));
    end
end

drugMetTable.NumberOfIncreasedT12VsT0 = sum((metFC12to0>=fcThresholdMetT12T0Pool).*...
                                                (metFDR12to0<=pThreshold),2);


good_mets = drugMetTable.GoodFilter == 1 &...
            drugMetTable.DrugMassFlag==0 & ...
            drugMetTable.NumberOfDrugs==1 &...
            drugMetTable.NumberOfIncreasedT12VsT0>0 &...
            drugMetTable.DrugConsumedFlag==1;    


select_mets = drug_coordinate_selected_drugs~=0 &...
              good_mets &...
              abs(metMZdifference_todrug)>massThreshold;
metMZdifference_todrug_all = metMZdifference_todrug(select_mets);
drug_coordinate_all = drug_coordinate_selected_drugs(select_mets);


selectedSidx = cellfun(@(x) contains(lower(x), 'theta'), metFCspecies);

drugmet_change12_0_selectesSpecies = sum((metFC12to0(:,selectedSidx)>=fcThresholdMetT12T0Pool).*...
                                                (metFDR12to0(:,selectedSidx)<=pThreshold),2);
drugmet_change12_0_selectesSpecies_value = max(metInt12(:,selectedSidx).*...
                                              (metFDR12to0(:,selectedSidx)<=pThreshold),[],2);

selected_species_change = drugmet_change12_0_selectesSpecies(select_mets);
metFoundFC_t12t0_max = log10(drugmet_change12_0_selectesSpecies_value(select_mets));

fig = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(metMZdifference_todrug_all,drug_coordinate_all,...
        'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', [.5 .5 .5])
set(gca, 'YTick', 1:length(selected_drugs_all));
set(gca, 'YTickLabel', selected_drugs_all);
hold on
minX = min(metMZdifference_todrug_all)-100;
maxX = max(metMZdifference_todrug_all)+100;
for i=1:length(selected_drugs_all)
    plot([minX maxX], [i i], 'k--')
end
plot([0 0], [0 length(selected_drugs_all)], '--k')
xlim([minX maxX])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add btheta specific metabolites to the plot
scattersize=2.3;
species_number_threshold = 0;
scatter(metMZdifference_todrug_all(selected_species_change>species_number_threshold),...
        drug_coordinate_all(selected_species_change>species_number_threshold),...
        (scattersize).^metFoundFC_t12t0_max(selected_species_change>species_number_threshold),...
        'r', 'LineWidth',2)
% add FC legends
scatter(350*ones(4,1), 0.5+[1:4], (scattersize).^[4 5 6 7], 'r', 'LineWidth',2) %   *scatterSize
text(350*ones(4,1), 0.5+[1:4], {'4' '5' '6' '7'})
text(300, 6, {'Intensity at t=12h, log10'})
title({'B.theta drug metabolites in all species(grey) and B.theta',...
       [num2str(nnz(selected_species_change>species_number_threshold)) 'B.theta metabolites']})
orient landscape
set(gca, 'XTick', -300:50:400)
xlim([-300 400])
print(fig, '-painters', '-dpdf', '-r600', '-bestfit', ...
    outfile_figED6A_btheta_drug_metabolites)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
