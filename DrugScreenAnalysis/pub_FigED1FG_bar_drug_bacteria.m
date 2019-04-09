%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED1FG_bar_drug_bacteria
% read drug fold chaneg info during incubation with bacteria
% plot barplots of number of drugs degraded by bacteria
% and number of bacteria degrading each drug
% for three different degradation thresholds (20%, 50%, 80%)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % takes as input outfile_Table3_drug_fold_changes
% outfile_figED1F_bar_bacteria_per_drug = ['Output' filesep 'FigEV1F_bar_degrading_species_per_drug.pdf'];
% outfile_figED1F_bar_bacteria_per_drug_zoomed = ['Output' filesep 'FigEV1F_bar_degrading_species_per_drug_only_degraded_drugs.pdf'];
% outfile_figED1G_bar_drug_per_bacteria = ['Output' filesep 'FigEV1G_bar_degraded_drugs_per_species.pdf'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read drug FC table
drugAnalysisResults = readtable(outfile_Table3_drug_fold_changes,...
                                'HeaderLines',1); %skip first line
drugAnalysisColumnNames = drugAnalysisResults.Properties.VariableNames;                            
%convert table to cell array
drugAnalysisResults = table2cell(drugAnalysisResults);
FCcolumns = cellfun(@(x) ~isempty(strfind(x, 'FC')) &...
                         isempty(strfind(x, 'FCSTD')) &...
                         isempty(strfind(x, 'DrugAdaptive')),...
                         drugAnalysisColumnNames);
FCSTDcolumns = cellfun(@(x) ~isempty(strfind(x, 'FCSTD')),...
                         drugAnalysisColumnNames);
PFDRcolumns = cellfun(@(x) ~isempty(strfind(x, 'p_FDR')),...
                         drugAnalysisColumnNames);                     
drugFC12to0_t0toCTRL_combined = cell2mat(drugAnalysisResults(:,FCcolumns))';
drugFC12to0_t0toCTRL_STDcombined = cell2mat(drugAnalysisResults(:,FCSTDcolumns))';
drugP12to0_t0toCTRL_combined = cell2mat(drugAnalysisResults(:,PFDRcolumns))';

DrugNames = drugAnalysisResults(:,ismember(drugAnalysisColumnNames, 'DrugName'));
drugFCadaptive = cell2mat(drugAnalysisResults(:,...
                          cellfun(@(x) ~isempty(strfind(x, 'DrugAdaptiveFC')),...
                          drugAnalysisColumnNames)));
%convert to log2FC
drugFCadaptive = abs(log2((100-drugFCadaptive)/100));
if size(drugFCadaptive,1)>1
    drugFCadaptive = drugFCadaptive';
end

totSpeciesNamesSingle = cellfun(@(x) x(strfind(x, 'FC')+2:end),...
                                     drugAnalysisColumnNames(FCcolumns), 'unif', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot number of species degrading every drug
fcThresholds = [0.8 0.5 0.2];
degradingSpecies = zeros(length(fcThresholds), size(drugFC12to0_t0toCTRL_combined,2));
for i=1:length(fcThresholds)
    curFCadaptive = drugFCadaptive;
    curFCadaptive = max(curFCadaptive, -log2(fcThresholds(i)));
    curFCadaptive = repmat(curFCadaptive,size(drugFC12to0_t0toCTRL_combined,1),1);
    
    degradingSpecies(i,:) = sum( (drugFC12to0_t0toCTRL_combined <= -curFCadaptive) &...
                        (drugP12to0_t0toCTRL_combined <= pThreshold) );% &...
                        %(drugFDRtoCtrl_t12 <= pThreshold) );
end
[~, idx] = sort(degradingSpecies(1,:));
degradingSpecies = degradingSpecies(:,idx);
degradedDrugs = DrugNames(idx);

% plot figure
mycolors = [0.8 0.8 0.8;...
            0.5 0.5 0.5;...
            0.2 0.2 0.2;]
fig = figure('units','normalized','outerposition',[0 0 1 1]);
hold on 
for i=1:length(fcThresholds)
    bar(degradingSpecies(i,:), 'FaceColor', mycolors(i,:))
end
set(gca, 'XTick', 1:length(degradingSpecies))
set(gca, 'XTickLabel', lower(degradedDrugs))
set(gca, 'fontSize', 6)
xticklabel_rotate([], 90);
title('Number of degrading species per drug')
ylabel('Number of degrading species, adaptive FC threshold')                  
legend({'>20%', '>50%', '>80%'}, 'Location', 'NorthWest')
orient landscape
print(fig, '-painters', '-dpdf', '-r600', '-bestfit', ...
    outfile_figED1F_bar_bacteria_per_drug)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = degradingSpecies(1,:)>0;
fig = figure('units','normalized','outerposition',[0 0 1 1]);
hold on 
for i=1:length(fcThresholds)
    bar(degradingSpecies(i,idx), 'FaceColor', mycolors(i,:))
end
set(gca, 'XTick', 1:nnz(idx))
set(gca, 'XTickLabel', lower(degradedDrugs(idx)))
set(gca, 'fontSize', 6)
xticklabel_rotate([], 90);
title('Number of degrading species per drug')
ylabel('Number of degrading species, adaptive FC threshold')                  
legend({'>20%', '>50%', '>80%'}, 'Location', 'NorthWest')
orient landscape
print(fig, '-painters', '-dpdf', '-r600', '-bestfit', ...
    outfile_figED1F_bar_bacteria_per_drug_zoomed)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot number of drugs degraded by every species
degradingSpecies = zeros(size(drugFC12to0_t0toCTRL_combined,1),length(fcThresholds));
for i=1:length(fcThresholds)
    curFCadaptive = drugFCadaptive;
    curFCadaptive = max(curFCadaptive, -log2(fcThresholds(i)));
    curFCadaptive = repmat(curFCadaptive, size(drugFC12to0_t0toCTRL_combined,1),1);
    
    degradingSpecies(:,i) = sum( (drugFC12to0_t0toCTRL_combined <= -curFCadaptive) &...
                        (drugP12to0_t0toCTRL_combined <= pThreshold), 2);% &...
                        %(drugFDRtoCtrl_t12 <= pThreshold),2 );
end
[~, idx] = sort(degradingSpecies(:,1));
degradingSpecies = degradingSpecies(idx,:);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
for i=1:length(fcThresholds)
    bar(degradingSpecies(:,i), 'FaceColor', mycolors(i,:))
end
set(gca, 'XTick', 1:length(degradingSpecies))
set(gca, 'XTickLabel', lower(totSpeciesNamesSingle(idx)))
set(gca, 'fontSize', 6)
xticklabel_rotate([], 90);
legend({'>20%', '>50%', '>80%'}, 'Location', 'NorthWest')
title('Number of degraded drugs per species')
ylabel('Number of drugs degraded, adaptive FC threshold')                  

orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
                outfile_figED1G_bar_drug_per_bacteria)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

