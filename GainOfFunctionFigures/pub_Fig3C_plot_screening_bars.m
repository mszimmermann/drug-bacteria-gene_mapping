%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_Fig3C_plot_screening_bars
% plot figure 3C: screening results from first Gain-of-Function screen
% in B.theta to search for gene products metabolizing diltiazem
% save figure to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% read following variables:
% inputFolderDataScreen1 data folder with raw data for plotting
% outfile_Fig3C_barplot_diltiazem_GoF_screen1 name of output figure file
% intensityNoise = 5000;

clear RawData

AssayedOD = {'OD333'};
TimePoints = {'T01' 'T02' 'T04' 'T06' 'T08' 'T12' 'T24'};
Time = [0 1 2 4 6 8 12 18 24];

fig = figure('units','normalized','outerposition',[0 0 1 1]);


for spi = 1:2
    subplot(2,1,spi)

    % import peak integrations from csv files, RT and AUC
    for i = 1:length(TimePoints)
           if spi==1
               TempStruct = ReadMixedTxt([inputFolderDataScreen1 ...
                                          sprintf([ 'MZ009F-', AssayedOD{1} '-ParentDrugs-', TimePoints{i} '.csv'])], ',');
           elseif spi==2
               TempStruct = ReadMixedTxt([inputFolderDataScreen1 ...
                                          sprintf(['MZ009F-', AssayedOD{1} '-Metabolites-', TimePoints{i} '.csv'])], ',');
           end
           RawData.(TimePoints{i}) = struct;
           TempStruct.IntensitiesRaw( TempStruct.IntensitiesRaw < intensityNoise)=intensityNoise; 
           RawData.(TimePoints{i}) = TempStruct;
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % correct intensities by internal standard
    IS = {'IS_YOH', 'IS_CAF', 'IS_SUL', 'IS_IPR'}; % internal standard with which to correctintensities
    for i = 1:length(TimePoints)
        IdxIS = ismember((cellfun(@(x) x(1:6), RawData.(TimePoints{i}).Compounds, 'UniformOutput', false)), IS);
        TempMatrix = RawData.(TimePoints{i}).IntensitiesRaw;
        TempIS = TempMatrix(:,IdxIS); %intensities of internal std of this time point
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
        TempSTD = std(TempInjectionSum);
        Outlayers = find(TempInjectionSum < (TempMean/2));
        TempMatrix(Outlayers, :)=NaN;

        RawData.(TimePoints{i}).IntensitiesCorrected = TempMatrix;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select diltiazem compounds 
    if spi==1
        compoundsIDX = find(cellfun(@(x) contains(x, 'DILTIAZEM'), RawData.T01.Compounds));
    elseif spi==2
        compoundsIDX = find(cellfun(@(x) contains(x, 'DILTIAZEM') &...
                                 contains(x, '372'), RawData.T01.Compounds));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create one matrix with samples versus metabolites
    RawDataCorrectedMat = zeros(length(fields(RawData))*size(RawData.T01.IntensitiesRaw,1),...
                                size(RawData.T01.IntensitiesRaw,2));
    RawDataSamples = cell(length(fields(RawData))*size(RawData.T01.IntensitiesRaw,1),1);
    RawDataTime = zeros(length(fields(RawData))*size(RawData.T01.IntensitiesRaw,1),1);
    idx = 1;
    for i=1:length(TimePoints)
        RawDataCorrectedMat(idx:idx+size(RawData.(TimePoints{i}).IntensitiesRaw,1)-1,:) = ...
            RawData.(TimePoints{i}).IntensitiesCorrected;
        RawDataSamples(idx:idx+size(RawData.(TimePoints{i}).IntensitiesRaw,1)-1) = ...
            cellfun(@(x) x(strfind(x,'Plate'):strfind(x,'OD')+4), RawData.(TimePoints{i}).SampleNames, 'unif',0);
        RawDataTime(idx:idx+size(RawData.(TimePoints{i}).IntensitiesRaw,1)-1) = str2double(TimePoints{i}(2:end));
        idx = idx+size(RawData.(TimePoints{i}).IntensitiesRaw,1);
    end
    RawDataSamples_unique = unique(RawDataSamples);
    RawDataTime_unique = unique(RawDataTime);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot stacked bar plots of compound change over time
    compoundNames = RawData.T01.Compounds;
    greyColor = [repmat(linspace(.9, 0.55, length(TimePoints)-1),3,1)';...
                    0.2 0.2 0.2];
    % plot timecourse for each plate
    plotDrugsFlag=1;

    for i=compoundsIDX
        legendarray = [];
        lidx = 1;
        legendplates = {};
        hold on
        curData = zeros(length(RawDataSamples_unique), length(RawDataTime_unique));
        for k = 1:length(RawDataTime_unique)
            curData(:,k) = (RawDataCorrectedMat(ismember(RawDataTime, RawDataTime_unique(k)),i));
        end
        plotData = zscore(curData);
        %outlier detection > or < than w*(q1-q3)
        w = 1.5;
        curData_outliers = zeros(size(plotData));
        for k=1:size(plotData,1)
            curData_outliers(k,:) = (plotData(k,:) > quantile(plotData(k,:),.75) + ...
                                                   w*(quantile(plotData(k,:),.75)-quantile(plotData(k,:),.25))) |...
                                    (plotData(k,:) < quantile(plotData(k,:),.25) - ...
                                                   w*(quantile(plotData(k,:),.75)-quantile(plotData(k,:),.25)));
        end
        curData_outliers(curData_outliers==1) = NaN;
        curData_outliers(curData_outliers==0) = curData(curData_outliers==0);
        %interpolate NaNs
        for i1 = 1:size(curData_outliers,1)
            for i2 = 1:size(curData_outliers,2)
                if isnan(curData_outliers(i1,i2))
                    if i2==1
                        nextnonan = find(~isnan(curData_outliers(i1,i2:end)));
                        curData_outliers(i1,i2) = curData_outliers(i1,i2+nextnonan(1)-1);
                    elseif i2==size(curData_outliers,2)
                        prevnonan = find(~isnan(curData_outliers(i1,1:i2)));
                        curData_outliers(i1,i2) = curData_outliers(i1,prevnonan(end));
                    else
                        nextnonan = find(~isnan(curData_outliers(i1,i2:end)));
                        prevnonan = find(~isnan(curData_outliers(i1,1:i2)));
                        if ~isempty(nextnonan) %there is a measurement on the right
                            curData_outliers(i1,i2) = (curData_outliers(i1,prevnonan(end))+...
                                                   curData_outliers(i1,i2+nextnonan(1)-1))/2;
                        else %take closest previous measurements
                            curData_outliers(i1,i2) = curData_outliers(i1,prevnonan(end));
                        end

                    end
                end
            end
        end

        hold on
        for k = 1:length(RawDataTime_unique)
            bar(1:length(curData_outliers(:,k)), curData_outliers(:,k), 'FaceColor', greyColor(k,:))
        end
        ylabel('Normalized intensity', 'FontSize', 12)

        %get the "center" of curData plots
        centerCurData = mean(nanmean(curData_outliers));
        % get the y limits
        curData_ylim = ylim();
        % get the ratio between "upper" and "lower" parts
        curData_ylim = (curData_ylim(2)-centerCurData)/(centerCurData-curData_ylim(1));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot z-scores at the right y axis
        yyaxis right

        plotData = mean(plotData,2);
        plot(plotData, 'Color', [84 134 186]/256, 'LineWidth', 2)
        ylabel('Mean Z-score across time points', 'FontSize', 12)

        xlim([0 length(curData_outliers(:,k))+1])
        legend(legendarray,TimePoints, 'FontSize', 12, 'Location', 'EastOutside')
        xlabel('Plate No', 'FontSize', 12 )

        % get y limits and make 0 at mean of curData
        curylim = ylim();
        if curylim(2)/curylim(1) < curData_ylim
            curylim(2) = round(abs(curylim(1))*curData_ylim);
        else
            curylim(2) = round(-abs(curylim(1))*(1-curData_ylim));
        end
        ylim(curylim)

        title(compoundNames{i},'FontSize', 12)

    end  
end
% print figure to file
orient landscape
print(fig, '-painters', '-dpdf', '-r600', '-bestfit', ...
        outfile_Fig3C_barplot_diltiazem_GoF_screen1)