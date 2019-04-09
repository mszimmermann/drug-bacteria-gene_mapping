%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_fig2a_volcano_plot.m
% read drug-bacteria untargeted metabolomics data
% for specified drug and species, 
% plot volcano plot of metabolites in the drug pools vs all other pools
% save volcano plot to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % takes as input infile_drug_experiment
% % takes as input infile_single_species
% infolder_drug_metabolite_data = ['DrugMetabolismProductsAnalysis' filesep 'DataFiles' filesep];
% pThresholdPooling = 10^(-6); % significance threshold for metabolite
% fcThresholdPool = 1; % fold change threshold between metabolite in drug pool vs other pools
% fig2a_volcano_poolingDrug = 'DILTIAZEM';
% fig2a_volcano_poolingSpecies = 'Bacteroides thetaiotaomicron';
% outfile_fig2a_volcano
% intensityNoise = 5000;

% load species information
AssaySpecies = table2cell(readtable(infile_single_species));
% Import general experimental parameter
ExperimentParameters = import_experiment_parameters(infile_drug_experiment);

ExperimentParameters.TimePoints = {'T00', 'T12'}; % define the timepoints to be extracted
% add control pools to the pooling scheme
addedPools = [97 98 99];
ExperimentParameters.PoolNumbers = [1:size(ExperimentParameters.PoolingScheme,2) addedPools];
ExperimentParameters.PoolingScheme = [ExperimentParameters.PoolingScheme ...
                                      zeros(size(ExperimentParameters.PoolingScheme,1), length(addedPools))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the specified species
AssayConditions = AssaySpecies(cellfun(@(x,y) contains(fig2a_volcano_poolingSpecies,x) &...
                                              contains(fig2a_volcano_poolingSpecies,y),...
                                              AssaySpecies(:,3), AssaySpecies(:,4)));
AssayConditions = cellfun(@(x) num2str(x), AssayConditions, 'unif',0);
AssayConditions = cellfun(@(x) strcat('S', repmat('0',1, 3-length(x)), x), AssayConditions, 'unif',0);

% get species names
AssayConditionsSpecies = cell(length(AssayConditions),1);
for i=1:length(AssayConditions)
    idx = cellfun(@(x) ismember(x, str2double(AssayConditions{i}(2:end))),...
                       AssaySpecies(:,1));
    AssayConditionsSpecies{i} = strjoin(AssaySpecies(idx,3:5),' ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create matrices of intensity fold change between drug-containing pools
% and other pools
fcMatrix_t0_cell = cell(length(AssayConditions),1);
fdrMatrix_t0_cell= cell(length(AssayConditions),1);
fcMatrix_t12_cell = cell(length(AssayConditions),1);
fdrMatrix_t12_cell = cell(length(AssayConditions),1);
curCompounds_cell = cell(length(AssayConditions),1);
curCompoundsSpectrum_cell = cell(length(AssayConditions),1);

for d = 1:length(AssayConditions)
    datafile = [infolder_drug_metabolite_data AssayConditions{d} '_RawDataListExportOptimized.csv'];
    % Import data of one experimental condition of the single species assays (24 pools, 2 timepoints)
    curIntensitiesRaw =  Import_Intensities(datafile, 6, inf);
    curCompounds =  Import_Compounds(datafile, 6, inf);
    curSamples = Import_Injections(datafile, 5, 5)';
    curCompoundsSpectrum =  Import_CompositeSpectrum(datafile, 6, inf);
    
    %check if intensities are likely on logarithmic scale
    if nnz(curIntensitiesRaw==0)>0 && nnz(curIntensitiesRaw>100000)<1000
        curIntensitiesRaw = 10.^(curIntensitiesRaw);
    end
    % find outliers due to bad injections, if sum of all compounds is
    % lower than the mean of all 48 data sets minus three times their std
    % threshold value
    Threshold = mean(sum(curIntensitiesRaw)-(3*(std(sum(curIntensitiesRaw)))));
    curExcludedDatasets = find(sum(curIntensitiesRaw) < Threshold);
    curIntensitiesNorm = curIntensitiesRaw;
    curIntensitiesNorm(:,curExcludedDatasets) = [];
    curSamples(curExcludedDatasets) = [];
    
    % exclude ions that were unique to the excluded dataset in all
    % parameters
    curIdxExcluded = find(sum(curIntensitiesNorm,2) ==0);
    curIntensitiesNorm(curIdxExcluded,:) =  [];
    curCompounds(curIdxExcluded,:) =  [];
    curCompoundsSpectrum(curIdxExcluded,:) =  [];
   
    curIntensitiesNorm(curIntensitiesNorm==1) = NaN;
    curIntensitiesNorm = quantilenorm(curIntensitiesNorm);
    curIntensitiesNorm(isnan(curIntensitiesNorm)) = intensityNoise; %noise level
 
    totTime = cellfun(@(x) str2double(x(28:29)), curSamples);
    totDataPools = cellfun(@(x) str2double(x(strfind(x,'Pool')+(4:5))), curSamples);

    totDataMat = curIntensitiesNorm';
    % Differential analysis of ions to identify drug metabolites
    tic
    [fcMatrix_t0, fdrMatrix_t0,...
     fcMatrix_t12, fdrMatrix_t12] = workflow_DifferentialAnalysisPool(totDataMat,...
                                                                       ExperimentParameters,...
                                                                       totTime,...
                                                                       totDataPools);
    toc
    fcMatrix_t0_cell{d} = fcMatrix_t0;
    fdrMatrix_t0_cell{d} = fdrMatrix_t0;
    fcMatrix_t12_cell{d} = fcMatrix_t12;
    fdrMatrix_t12_cell{d} = fdrMatrix_t12;
    curCompounds_cell{d} = curCompounds;
    curCompoundsSpectrum_cell{d} = curCompoundsSpectrum;
    fprintf('Calculated fold changes between pools for %s (%s)\n', AssayConditions{d}, AssayConditionsSpecies{d});
    
end
clear datafile curIntensitiesRaw curCompounds curSamples curCompoundsSpectrum
clear Threshold curExcludedDatasets curIntensitiesNorm curIdxExcluded 
clear totTime totDataPools totDataMat 
clear fcMatrix_t0 fdrMatrix_t0 fcMatrix_t12 fdrMatrix_t12 fcMatrix_t0_ind fcMatrix_t12_ind

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot volcano plots for t12 and t0 for each drug and species
drugidx = find(cellfun(@(x) contains(upper(x), fig2a_volcano_poolingDrug),...
               ExperimentParameters.DrugNames));

for i=1:length(fcMatrix_t12_cell)
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    % plot pool comparison an T=0
    sp = subplot(1,2,1);
    metFC = log2(fcMatrix_t0_cell{i}(:, drugidx));
    metPadj = fdrMatrix_t0_cell{i}(:, drugidx);
    metMass = cellfun(@(x) x(1:strfind(x,'@')-1), curCompounds_cell{i}, 'unif',0); 
    [upreg] = plotVolcano(metFC, metPadj, curCompounds_cell{i}, fcThresholdPool, pThresholdPooling, [], sp);
    % add mass labels to increased metabolites
    text(metFC(upreg), -log10(metPadj(upreg)), metMass(upreg))
    title({['T=0h ' ExperimentParameters.DrugNames{drugidx}],...
           AssayConditionsSpecies{i}});
    axis square
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot pool comparison an T=12
    sp = subplot(1,2,2);
    metFC = log2(fcMatrix_t12_cell{i}(:, drugidx));
    metPadj = fdrMatrix_t12_cell{i}(:, drugidx);
    metMass = cellfun(@(x) x(1:strfind(x,'@')-1), curCompounds_cell{i}, 'unif',0); 
    [upreg] = plotVolcano(metFC, metPadj, curCompounds_cell{i}, fcThresholdPool, pThresholdPooling, [], sp);
    % add mass labels to increased metabolites
    text(metFC(upreg), -log10(metPadj(upreg)), metMass(upreg))
    title({['T=12h ' ExperimentParameters.DrugNames{drugidx}],...
           AssayConditionsSpecies{i}});
    axis square
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save figure to file
    orient landscape
    print(fig, '-painters', '-dpsc', '-r600', '-bestfit', '-append',...
          outfile_fig2a_volcano)
    close(fig)
end    
