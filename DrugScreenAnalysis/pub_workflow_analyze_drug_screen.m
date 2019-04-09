%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Analyse 21 drug pools of drug screen of the Single Species Experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_workflow_analyze_drug_screen.m
% read drug-bacteria screen data
% calculate fold changes between t=12h and t=0h
% identify fast metabolizers (drug comsumed at t=0h compared to control)
% save data to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% will use the following variables:
% infile_single_species = 'Table1_AssayedSingleSpecies.csv';
% infile_drug_experiment = 'Table2_DrugScreenInfo.csv';
% infolder_drug_data = ./Data/
% outfile_fig1B_PCA = 'OutputFiles\fig1B_plot_rdkit_drugbank_selected_drugs_PCA.pdf';
% intensityNoise = 5000;
% fcThresholdDrug = abs(log2(0.8)); % 20% reduction
% pThreshold = 0.05;


%load single spesies information 
AssaySpecies = table2cell(readtable(infile_single_species));

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

% import peak integrations from csv files, RT and AUC
PoolNo = 1:21;
TimePoints = {'T00','T12'};
for i = PoolNo % i = pools to be integrated
        TempStruct = ReadMixedTxt([infolder_drug_data sprintf('MZ002H_Pool%02.0f_Final.csv', i)], ',');
        PoolCount = sprintf('Pool%02.0f', i);
        RawData.(PoolCount) = struct;
        % set lowe intensity to defined noise level intensityNoise
        TempStruct.IntensitiesRaw( TempStruct.IntensitiesRaw < intensityNoise ) = intensityNoise; 
        RawData.(PoolCount) = TempStruct;
end

% internal standard with which to correct intensities
IS = {'IS_YOH', 'IS_CAF', 'IS_SUL', 'IS_IPR'}; 
for i = 1:length(PoolNo)
    IdxIS = ismember((cellfun(@(x) x(1:6), RawData.(sprintf('Pool%02.0f', i)).Compounds, 'UniformOutput', false)), IS);
    PoolCount = sprintf('Pool%02.0f', PoolNo(i));
    TempMatrix = RawData.(PoolCount).IntensitiesRaw;
    TempIS = TempMatrix(:,IdxIS); %intensities of internal std of this plate and time point
    MeanIS = mean(TempIS);
    CorrectionMatrix = ones(size(TempIS));
    for p = 1:length(MeanIS)
        CorrectionMatrix(:,p)= TempIS(:,p)./MeanIS(p); %calcuate fold change from the mean
    end
    Correction = mean(CorrectionMatrix,2);
    for p = 1: length(Correction)
        TempMatrix(p,:)= TempMatrix(p,:)./Correction(p);
    end
    % eliminate samples with bad injections
    TempInjectionSum = sum(TempMatrix,2);
    TempMean = mean(TempInjectionSum);
    Outliers = TempInjectionSum < (TempMean/2);
    TempMatrix(Outliers, :)=NaN;

    RawData.(PoolCount).IntensitiesCorrected = TempMatrix;
end
clear IS IdxIS PoolCount TempMatrix TempIS MeanIS CorrectionMatrix 
clear Correction TempInjectionSum TempMean TempSTD Outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make one large matrix with intensities (corrected)
totCompounds = {};
totMeasurements = 0;
for i = PoolNo % i = pools to be integrated
    PoolCount = sprintf('Pool%02.0f', i);
    totCompounds = union(totCompounds, RawData.(PoolCount).Compounds);
    totMeasurements = totMeasurements + size(RawData.(PoolCount).IntensitiesCorrected,1);
end
totDataMat = zeros(totMeasurements, length(totCompounds));
totDataRT = zeros(totMeasurements, length(totCompounds));
totDataPools = zeros(totMeasurements, 1);
totDataSampleNames = cell(totMeasurements, 1);
idx = 1;
for i = PoolNo % i = pools to be 
    PoolCount = sprintf('Pool%02.0f', i);
    % intersect current compounds with total compounds for index matching
    [~,totidx,compidx] = intersect(totCompounds, RawData.(PoolCount).Compounds, 'stable');
    % add current data to the total daza
    totDataMat(idx:idx+size(RawData.(PoolCount).IntensitiesCorrected,1)-1,totidx) = ...
        RawData.(PoolCount).IntensitiesCorrected(:, compidx);
    % combine RT matrices
    totDataRT(idx:idx+size(RawData.(PoolCount).RT,1)-1,totidx) = ...
        RawData.(PoolCount).RT(:, compidx);
    %combine sample names
    totDataSampleNames(idx:idx+size(RawData.(PoolCount).IntensitiesCorrected,1)-1) =...
        RawData.(PoolCount).SampleNames;
    %record the pools of the data
    totDataPools(idx:idx+size(RawData.(PoolCount).IntensitiesCorrected,1)-1) = i;
    %increment the index of the total matrix
    idx = idx+size(RawData.(PoolCount).IntensitiesCorrected,1);
end

% find outlier data due to bad injections, if sum of all compounds is
% lower than the mean of all 48 data sets minus two times their std 
% (=threshold value)
totMean = mean(totDataMat(:));
totSTD = std(totDataMat(:));
badInjections = zeros(size(totDataMat,1),1);
for i=1:size(totDataMat,1)
    badInjections(i) = sum(totDataMat(i,:))<(totMean-2*totSTD);
end
% remove bad injections
totDataMat(badInjections==1,:) = [];
totDataRT(badInjections==1,:) = [];
totDataPools(badInjections==1) = [];
totDataSampleNames(badInjections==1) = [];

% set 0 (not detected compounds) to min intensity
totDataMat(totDataMat==0) = intensityNoise;
% clear previos structures
clear RawData totidx totMean totMeasurements totSTD
clear badInjections compidx idx p PoolCount 
clear TempInjectionSum TempIS TempMatrix TempMean TempSTD TempStruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get species and time info from the sample names

totSpecies = cellfun(@(x) x(22:25), totDataSampleNames, 'unif', 0);
totTime = cellfun(@(x) str2double(x(28:29)), totDataSampleNames);
totSpecies_unique = unique(totSpecies);
% rename pericyazine to match the name used in the drug table
totCompounds = cellfun(@(x) strrep(x, 'PERICIAZINE', 'PERICYAZINE'), totCompounds, 'unif', 0);  
totCompounds_abbr = cellfun(@(x) x(1:6), totCompounds, 'unif', 0);
%distingush betamethasone acetate and valerate
totCompounds_abbr{ismember(totCompounds,...
                           'BETAMETHASONE ACETATE Results')} = 'BETAMA';
totCompounds_abbr{ismember(totCompounds,...
                           'BETAMETHASONE VALERATE Results')} = 'BETAMV';
%distingush trimetazidine and trimethobenzamide
totCompounds_abbr{ismember(totCompounds,...
                           'TRIMETAZIDINE DIHYDROCHLORIDE Results')} = 'TRIMET';
totCompounds_abbr{ismember(totCompounds,...
                           'TRIMETHOBENZAMIDE HYDROCHLORIDE Results')} = 'TRIMEH';


SpeciesNo = cellfun(@(x) strcat('S',num2str(x, '%03d')), AssaySpecies(:,1), 'unif', 0);

[~, ~, totSpeciesNames] = intersect(totSpecies_unique, SpeciesNo, 'stable');
totSpeciesNames = AssaySpecies(totSpeciesNames,3:5);
totSpeciesNamesSingle = strcat(totSpeciesNames(:,1), {' '},...
                               totSpeciesNames(:,2), {' '},...
                               totSpeciesNames(:,3));

% Get the control at pH 7
totCtrl = cellfun(@(x) ~isempty(strfind(x, 'Control pH 7')), totSpeciesNamesSingle);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create fold change table for each species and drug between time 12 and 0
drugFC12to0 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));
drugFCSTD12to0 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));

drugP12to0 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));
% define cutoffs to detect one-out-of-four outliers in the drug pools
% defined in Drug_bacteria_gene_mapping_variables
%drug_outlierSTDcutoff = 5;
%drug_outlierMEANcutoff = 0.2;

for i=1:length(totSpecies_unique)
     for j=1:length(ExperimentParameters.DrugNamesAbbr)
        IdxPool = find(ExperimentParameters.PoolingScheme(j,:)==1);
        curdrugidx = ismember(totCompounds_abbr, ExperimentParameters.DrugNamesAbbr{j});
        curData_t0 = totDataMat(ismember(totSpecies, totSpecies_unique{i}) &...
                                totTime==0 &...
                                ismember(totDataPools, IdxPool),:);
        curData_t12 = totDataMat(ismember(totSpecies, totSpecies_unique{i}) &...
                                 totTime==12 &...
                                 ismember(totDataPools, IdxPool),:);
        outliersData_t0 = detect_one_outlier(curData_t0(:,curdrugidx), drug_outlierSTDcutoff, drug_outlierMEANcutoff);
        outliersData_t12 = detect_one_outlier(curData_t12(:,curdrugidx), drug_outlierSTDcutoff, drug_outlierMEANcutoff);
       
        drugFC12to0(i,j) = nanmean(curData_t12(~outliersData_t12,curdrugidx)) / ...
                           nanmean(curData_t0(~outliersData_t0,curdrugidx));
        [~, drugP12to0(i,j)] = ttest2(curData_t12(~outliersData_t12,curdrugidx),...
                                 curData_t0(~outliersData_t0,curdrugidx),...
                                 'VarType', 'unequal');
        drugFCSTD12to0(i,j) = abs(drugFC12to0(i,j)).*...
                               sqrt( (nanstd(curData_t12(~outliersData_t12,curdrugidx))./...
                                        nanmean(curData_t12(~outliersData_t12,curdrugidx))).^2 + ...
                                     (nanstd(curData_t0(~outliersData_t0,curdrugidx))./nanmean(curData_t0(~outliersData_t0,curdrugidx))).^2 )';
     end
     fprintf('Calculated fold changes t=12h to t=0h for %s (%s)\n', totSpecies_unique{i}, totSpeciesNamesSingle{i});
end
clear IdxPool curData_t0 curData_t12 curdrugidx outliersData_t0 outliersData_t12
    
% adjust all drugs for multiple hypotheses testing 
% with Benjamini-Hochberg procedure
drugFDR12to0 = mafdr(drugP12to0(:),'BHFDR', 1);
drugFDR12to0 = reshape(drugFDR12to0, size(drugFC12to0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% to detect fast metabolizers (already at t=0), 
%%%% perform t-test to the changes of the no-bacteria control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create fold change table for each species and drug and time point

drugFCtoCtrl_t0 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));
drugFCtoCtrl_t12 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));
drugPtoCtrl_t0 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));
drugPtoCtrl_t12 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));
drugFCtoCtrlSTD_t0 = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));

for i=1:length(totSpecies_unique)
     for j=1:length(ExperimentParameters.DrugNamesAbbr)
        IdxPool = find(ExperimentParameters.PoolingScheme(j,:)==1);
        curdrugidx = ismember(totCompounds_abbr, ExperimentParameters.DrugNamesAbbr{j});
        curData_t0 = totDataMat(ismember(totSpecies, totSpecies_unique{i}) &...
                                totTime==0 &...
                                ismember(totDataPools, IdxPool),:);
        curData_t12 = totDataMat(ismember(totSpecies, totSpecies_unique{i}) &...
                                 totTime==12 &...
                                 ismember(totDataPools, IdxPool),:);
        curCtrl_t0 = totDataMat(ismember(totSpecies, totSpecies_unique(totCtrl)) &...
                                totTime==0 &...
                                ismember(totDataPools, IdxPool),:);
        curCtrl_t12 = totDataMat(ismember(totSpecies, totSpecies_unique(totCtrl)) &...
                                 totTime==12 &...
                                 ismember(totDataPools, IdxPool),:);
        outliersData_t0 = detect_one_outlier(curData_t0(:,curdrugidx), drug_outlierSTDcutoff, drug_outlierMEANcutoff);
        outliersCtrl_t0 = detect_one_outlier(curCtrl_t0(:,curdrugidx), drug_outlierSTDcutoff, drug_outlierMEANcutoff);
        
        drugFCtoCtrl_t0(i,j) = nanmean(curData_t0(~outliersData_t0,curdrugidx)) /...
                          nanmean(curCtrl_t0(~outliersCtrl_t0,curdrugidx));
        [~, drugPtoCtrl_t0(i,j)] = ttest2(curData_t0(~outliersData_t0,curdrugidx),...
                                          curCtrl_t0(~outliersCtrl_t0,curdrugidx),...
                                          'VarType', 'unequal');
        drugFCtoCtrlSTD_t0(i,j) = abs(drugFCtoCtrl_t0(i,j)).*...
                               sqrt( (nanstd(curData_t0(~outliersData_t0,curdrugidx))./...
                                        nanmean(curData_t0(~outliersData_t0,curdrugidx))).^2 + ...
                                     (nanstd(curCtrl_t0(~outliersCtrl_t0,curdrugidx))./nanmean(curCtrl_t0(~outliersCtrl_t0,curdrugidx))).^2 )';

                                 
        outliersData_t12 = detect_one_outlier(curData_t12(:,curdrugidx), drug_outlierSTDcutoff, drug_outlierMEANcutoff);
        outliersCtrl_t12 = detect_one_outlier(curCtrl_t12(:,curdrugidx), drug_outlierSTDcutoff, drug_outlierMEANcutoff);
      
        drugFCtoCtrl_t12(i,j) = nanmean(curData_t12(~outliersData_t12,curdrugidx)) /...
                          nanmean(curCtrl_t12(~outliersCtrl_t12,curdrugidx));
        [~, drugPtoCtrl_t12(i,j)] = ttest2(curData_t12(~outliersData_t12,curdrugidx),...
                                          curCtrl_t12(~outliersCtrl_t12,curdrugidx),...
                                          'VarType', 'unequal');
     end
     fprintf('Calculated fold changes to control for %s (%s)\n', totSpecies_unique{i}, totSpeciesNamesSingle{i});
end
clear IdxPool curData_t0 curData_t12 curdrugidx 
clear curCtrl_t0 curCtrl_t12 outliersData_t12  outliersData_t0 outliersCtrl_t0 outliersCtrl_t12

% adjust all drugs for multiple hypotheses testing 
drugFDRtoCtrl_t0 = mafdr(drugPtoCtrl_t0(:),'BHFDR', 1);
drugFDRtoCtrl_t0 = reshape(drugFDRtoCtrl_t0, size(drugFCtoCtrl_t0));

drugFDRtoCtrl_t12 = mafdr(drugPtoCtrl_t12(:),'BHFDR', 1);
drugFDRtoCtrl_t12 = reshape(drugFDRtoCtrl_t12, size(drugFCtoCtrl_t12));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add species that are fast metabolizers
% calculate the maximum fold change which is possible from t0 to lowest
% intensity

drugFC12to0_max = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));
drugFC12to0_t0toCTRL = zeros(length(totSpecies_unique), length(ExperimentParameters.DrugNamesAbbr));

for i=1:length(totSpecies_unique)
     for j=1:length(ExperimentParameters.DrugNamesAbbr)
        IdxPool = find(ExperimentParameters.PoolingScheme(j,:)==1);
        curdrugidx = ismember(totCompounds_abbr, ExperimentParameters.DrugNamesAbbr{j});
        curData_t0 = totDataMat(ismember(totSpecies, totSpecies_unique{i}) &...
                                totTime==0 &...
                                ismember(totDataPools, IdxPool),curdrugidx);
        curData_t12 = intensityNoise*ones(size(curData_t0));
        drugFC12to0_max(i,j) = nanmean(curData_t12) / nanmean(curData_t0);
        
        curCtrl_t0 = totDataMat(ismember(totSpecies, totSpecies_unique(totCtrl)) &...
                                totTime==0 &...
                                ismember(totDataPools, IdxPool),curdrugidx);
        drugFC12to0_t0toCTRL(i,j) = (abs(log2(drugFC12to0(i,j))) - ...
                                     abs(log2(drugFCtoCtrl_t0(i,j)))) .* ...
                                     (log2(drugFC12to0(i,j))<0 & log2(drugFCtoCtrl_t0(i,j))<0);
     end
end
clear IdxPool curData_t0 curData_t12 curdrugidx


drugFC12to0_max = log2(drugFC12to0_max);

% find drugs for which the difference between FC to ctrl at t0 and FC at
% t12 vs t0 is more than 2^5
[~, fastDrugs_idx] = ind2sub(size(drugFC12to0_t0toCTRL),...
                             find(drugFC12to0_t0toCTRL < drug_fastMetabolizerThreshold));
fastDrugs_idx = unique(fastDrugs_idx);
% for each of this drug, for each species choose max significant negative foldchange
% between FC to ctrl at t0 and FC at t12 vs t0
drugFC12to0_t0toCTRL_combined = drugFC12to0;
drugFC12to0_t0toCTRL_STDcombined = drugFCSTD12to0;

drugP12to0_t0toCTRL_combined = drugFDR12to0;

for idx=1:length(fastDrugs_idx)
    j = fastDrugs_idx(idx);
    change_fc = (log2(drugFCtoCtrl_t0(:,j))<=-fcThresholdDrug) &...
                (drugFDRtoCtrl_t0(:,j)<=pThreshold) &...
                (log2(drugFC12to0(:,j))>log2(drugFCtoCtrl_t0(:,j)));
    drugFC12to0_t0toCTRL_combined(change_fc,j) = drugFCtoCtrl_t0(change_fc,j);
    drugFC12to0_t0toCTRL_STDcombined(change_fc,j) = drugFCtoCtrlSTD_t0(change_fc,j);
    drugP12to0_t0toCTRL_combined(change_fc,j) = drugFDRtoCtrl_t0(change_fc,j);
end
% convert to log2 std
drugFC12to0_t0toCTRL_STDcombined = (drugFC12to0_t0toCTRL_STDcombined./...
                                    (drugFC12to0_t0toCTRL_combined*log(2)));
drugFC12to0_t0toCTRL_combined = log2(drugFC12to0_t0toCTRL_combined);                                
clear IdxPool curData_t0 curData_t12 curdrugidx change_fc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive threshold = mean+2*std of positive changes
drugFCadaptive = zeros(1,size(drugFC12to0,2));
for j=1:length(ExperimentParameters.DrugNamesAbbr)
    curFC = drugFC12to0_t0toCTRL_combined(:,j);
    maxThreshold = 0 + mean(curFC(curFC>0)) + 2*std(curFC(curFC>0));
    drugFCadaptive(j) = max(fcThresholdDrug, maxThreshold); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save drug fold change, std and p-values to file
fid = fopen(outfile_Table3_drug_fold_changes, 'w');
fprintf(fid, ';');
for i=1:length(totSpeciesNamesSingle)
    fprintf(fid, ';%s;;;;', totSpeciesNamesSingle{i});
end
fprintf(fid, '\n');
fprintf(fid, 'DrugName;Drug adaptive FC threshold %%'); 
for i=1:length(totSpeciesNamesSingle)
    fprintf(fid, ';%% consumed %s;%% consumed STD %s;FC %s;FC STD %s; p(FDR) %s', ...
                totSpeciesNamesSingle{i},totSpeciesNamesSingle{i},...
                totSpeciesNamesSingle{i}, totSpeciesNamesSingle{i},...
                totSpeciesNamesSingle{i});
end
fprintf(fid, '\n');
for j=1:size(drugFC12to0_t0toCTRL_combined,2)
    fprintf(fid,'%s',ExperimentParameters.DrugNames{j}); 
    fprintf(fid, ';%.3f', 100*(1-2.^(-drugFCadaptive(j))));
    for i=1:length(totSpeciesNamesSingle)
        fprintf(fid, ';%.3f;%.3f;%.3f;%.3f;%.3f', ...
                    max(0, 100*(1-2^drugFC12to0_t0toCTRL_combined(i,j))),...
                    100*((2.^drugFC12to0_t0toCTRL_combined(i,j))*log(2).*drugFC12to0_t0toCTRL_STDcombined(i,j)),...
                    drugFC12to0_t0toCTRL_combined(i,j),...
                    drugFC12to0_t0toCTRL_STDcombined(i,j),...
                    drugP12to0_t0toCTRL_combined(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


