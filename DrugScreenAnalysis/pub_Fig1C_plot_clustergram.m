%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_Fig1C_plot_clustergram
% read drug FC table
% plot clustergram and save to pdf file and to a table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables;
% takes as input outfile_Table3_drug_fold_changes
%outfile_fig1c_clustergram = ['Output' filesep 'Fig1C_clustergram_drug_percentChange.pdf'];
%outfile_table_FCclustergram = ['Output' filesep 'Table3_drug_fold_changes_in_clustergram_order.csv'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
drugFCadaptive = abs(log2((100-drugFCadaptive)/100))';

totSpeciesNamesSingle = cellfun(@(x) x(strfind(x, 'FC')+2:end),...
                                     drugAnalysisColumnNames(FCcolumns), 'unif', 0);

drugFC12to0_trunc = drugFC12to0_t0toCTRL_combined;
drugP_trunc = drugP12to0_t0toCTRL_combined;

clustDrugNames = DrugNames;
clustdrugFCadaptive = drugFCadaptive;

clustSpecies = totSpeciesNamesSingle;

%remove controls from the clustergram
drugFC12to0_trunc(cellfun(@(x) ~isempty(strfind(x, 'Control')), clustSpecies),:) = [];
drugP_trunc(cellfun(@(x) ~isempty(strfind(x, 'Control')), clustSpecies),:) = [];
clustSpecies(cellfun(@(x) ~isempty(strfind(x, 'Control')), clustSpecies)) = [];

removeDrugsFromClustergram = sum( (drugFC12to0_trunc<=-repmat(drugFCadaptive, size(drugFC12to0_trunc,1),1)) &...
                                  (drugP_trunc<=pThreshold));
drugFC12to0_trunc(:, removeDrugsFromClustergram==0) = [];
drugP_trunc(:,removeDrugsFromClustergram==0) = [];
clustDrugNames(removeDrugsFromClustergram==0) = [];
clustdrugFCadaptive(removeDrugsFromClustergram==0) = [];


drugP_trunc(isnan(drugP_trunc)) = 1;
drugFC12to0_trunc(isnan(drugFC12to0_trunc)) = 0;

clustdist = 'euclidean';

cmaplength = 8;
mycmap = repmat([0 0 linspace(0.4,1,cmaplength)]',1,3);

drugFC12to0_trunc = (drugFC12to0_trunc.*(drugP_trunc <= pThreshold)...
                    .*(drugFC12to0_trunc<=-repmat(clustdrugFCadaptive, size(drugFC12to0_trunc,1),1)));
% convert fold change to % consumed drug                
drugFC12to0_trunc = 1-2.^drugFC12to0_trunc;

cgo=clustergram(drugFC12to0_trunc,...
                             'DisplayRange', 1,...truncVal,...
                             'RowLabels', clustSpecies,...
                             'ColumnLabels', clustDrugNames,...
                             'Symmetric', 0,...
                             'ColumnPdist', clustdist, 'RowPdist', clustdist,...
                             'Colormap', mycmap);%, 'DisplayRatio', 0.6);

                         
                         
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 4)
% insert colorbar
cbButton = findall(gcf,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback');
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
cb  = findobj(gcf,'Tag','HeatMapColorbar');
cb.Label.String = 'Percent of drug consumed after 12h';

% save to file
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', outfile_fig1c_clustergram)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save drug fold change, std and p-values to file in the clustergram order
[plotSpecies, ~, species_idxdrug] = intersect(cgo.RowLabels, totSpeciesNamesSingle, 'stable');
[plotDrugs, ~, drugs_idxdrug] = intersect(cgo.ColumnLabels,DrugNames, 'stable');
drugs_idxdrug = flipud(drugs_idxdrug);

fid = fopen(outfile_table_FCclustergram, 'w');
fprintf(fid, ';;');
for i=1:length(species_idxdrug)
    fprintf(fid, ';%s;;;;', totSpeciesNamesSingle{species_idxdrug(i)});
end
fprintf(fid, '\n');
fprintf(fid, 'DrugName;Drug adaptive FC threshold %%'); 
for i=1:length(species_idxdrug)
    fprintf(fid, ';%% consumed %s;%% consumed STD %s;FC %s;FC STD %s; p(FDR) %s', ...
                totSpeciesNamesSingle{species_idxdrug(i)},totSpeciesNamesSingle{species_idxdrug(i)},...
                totSpeciesNamesSingle{species_idxdrug(i)}, totSpeciesNamesSingle{species_idxdrug(i)},...
                totSpeciesNamesSingle{species_idxdrug(i)});
end
fprintf(fid, '\n');
for j=1:length(drugs_idxdrug)
    fprintf(fid,'%s',DrugNames{drugs_idxdrug(j)}); 
    fprintf(fid, ';%.3f', drugFCadaptive(drugs_idxdrug(j)));
    for i=1:length(species_idxdrug)
        fprintf(fid, ';%.3f;%.3f;%.3f;%.3f;%.3f', ...
                    100*(1-2^drugFC12to0_t0toCTRL_combined(species_idxdrug(i),drugs_idxdrug(j))),...
                    100*((2.^drugFC12to0_t0toCTRL_combined(species_idxdrug(i),drugs_idxdrug(j)))*log(2).*drugFC12to0_t0toCTRL_STDcombined(species_idxdrug(i),drugs_idxdrug(j))),...
                    drugFC12to0_t0toCTRL_combined(species_idxdrug(i),drugs_idxdrug(j)),...
                    drugFC12to0_t0toCTRL_STDcombined(species_idxdrug(i),drugs_idxdrug(j)),...
                    drugP12to0_t0toCTRL_combined(species_idxdrug(i),drugs_idxdrug(j)));
    end
    fprintf(fid, '\n');
end
fclose(fid);