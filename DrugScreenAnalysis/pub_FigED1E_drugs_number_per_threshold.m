%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED1F_drugs_number_per_threshold
% read degraded drug fold change table and 
% plot how many drugs degraded depending on the degradation threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % takes as input outfile_Table3_drug_fold_changes
% outfile_figED1E_drug_number_per_threshold = ['Output' filesep 'FigE1E_drug_changes_vs_PercentFCthreshold.pdf'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read drug FC info
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
% calculate number of drugs changing at certain thresholds
pThreshold = 0.05;
changeTresholds = [0:0.1:0.9 0.95 0.99 0.999];
changeDrugs = zeros(size(changeTresholds));
changeDrugBacteriaPairs = zeros(size(changeTresholds));
for i=1:length(changeTresholds)
    curFCadaptive = drugFCadaptive;
    curFCadaptive(curFCadaptive==abs(log2(0.8))) = changeTresholds(i);
    curFCadaptive = max(1-2.^(-curFCadaptive), changeTresholds(i));
    curFCadaptive = repmat(curFCadaptive, size(drugFC12to0_t0toCTRL_combined,1),1);
    
    nocontrolFC = 1-2.^drugFC12to0_t0toCTRL_combined;
    nocontrolFC(ismember(totSpeciesNamesSingle(:,3), 'Control'),:) = 0;
    
    changeDrugs(i) = nnz(sum( (nocontrolFC>=curFCadaptive) &...
                              (drugP12to0_t0toCTRL_combined<pThreshold) ));
    changeDrugBacteriaPairs(i) = nnz((nocontrolFC>=curFCadaptive) &...
                                     (drugP12to0_t0toCTRL_combined<pThreshold) );
end
changeDrugs(1) = size(drugFC12to0_t0toCTRL_combined,2);

figure
plot(100*changeTresholds, changeDrugs, '-o')
xlabel('Drug consumption, %')
ylabel('Number of drugs')
axis square
title('Significant changes depending on FC threshold')

% text number of drugs for 20%, 50% and 80%
text(20, changeDrugs(changeTresholds==0.2), num2str(changeDrugs(changeTresholds==0.2)))
text(50, changeDrugs(changeTresholds==0.5), num2str(changeDrugs(changeTresholds==0.5)))
text(80, changeDrugs(changeTresholds==0.8), num2str(changeDrugs(changeTresholds==0.8)))

% save to file
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            outfile_figED1E_drug_number_per_threshold)
