%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_Fig3D_plot_screening_heatmaps
% plot figure 3C: screening results from second Gain-of-Function screen
% in B.theta to search for gene products metabolizing diltiazem
% save figure to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% read following variables:
% inputFolderDataScreen2 data folder with raw data for plotting
% outfile_Fig3D_heatmap_diltiazem_GoF_screen2 name of output figure file
% intensityNoise = 5000;

clear RawData
AssayedOD = {'Plate6'};
TimePoints = {'T01' 'T04' 'T02' 'T12' 'T24'};

greyColor = [repmat(linspace(.9, 0.55, length(TimePoints)-1),3,1)';...
                    0.2 0.2 0.2];
                
for spi=1:2
    
    % import peak integrations from csv files, RT and AUC
    for i = 1:length(TimePoints)
        if spi==1
           TempStruct = ReadMixedTxt([inputFolderDataScreen2...
                                      sprintf(['4Pool-', AssayedOD{1} ,'-' TimePoints{i} '-Drugs' '.csv'])], ',');
        elseif spi==2
           TempStruct = ReadMixedTxt([inputFolderDataScreen2...
                                      sprintf(['4Pool-', AssayedOD{1} ,'-' TimePoints{i} '-Metabolites' '.csv'])], ',');
        end
       RawData.(TimePoints{i}) = struct;
       TempStruct.IntensitiesRaw( TempStruct.IntensitiesRaw < intensityNoise)=intensityNoise; % detection limit
       RawData.(TimePoints{i}) = TempStruct;
    end
  
    % correct by internal standard
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
    
    % select compounds
    compoundNames = RawData.T01.Compounds';
    if spi==1
        compoundsIDX = find(cellfun(@(x) contains(x, 'DILTIAZEM'), compoundNames));
    elseif spi==2
        compoundsIDX = find(cellfun(@(x) contains(x, 'DILTIAZEM') &...
                                 contains(x, '372'), compoundNames));
    end
    
    i=compoundsIDX;

    plotDrugsFlag=1;

    fig = figure('units','normalized','outerposition',[0 0 1 1]);
   
    % plot timecourse for each plate
    for mylabels = 1:2
        if mylabels == 1
            myxlabel = 'Columns';
            RawDataSamples_unique = RawData.T01.SampleNames;
            currentSamplesIDX = find(cellfun(@(x) contains(x, '-C'), RawDataSamples_unique));
            subplot(2,1,1);
        elseif mylabels == 2
            myxlabel = 'Rows';
            RawDataSamples_unique = RawData.T01.SampleNames;
            currentSamplesIDX = find(cellfun(@(x) contains(x, '-R'), RawDataSamples_unique));
            subplot(2,1,2);
        end
        [currentSamples, sortidx] = sort(RawDataSamples_unique(currentSamplesIDX));
        currentSamplesIDX = currentSamplesIDX(sortidx);

        for i=compoundsIDX %16 diltiazem compoundsIDX
            legendarray = [];
            lidx = 1;
            legendplates = {};
            hold on
            curData = zeros(length(currentSamplesIDX), length(TimePoints));
            for k = 1:length(TimePoints)
                curData(:,k) = RawData.(TimePoints{k}).IntensitiesCorrected(currentSamplesIDX,i);
            end
            plotData = zscore(curData);
     
            % detect outliers in intensity and not z-score
            w = 1.5;
            curData_outliers = zeros(size(plotData));
            for k=1:size(plotData,1)
                curData_outliers(k,:) = (curData(k,:) > quantile(curData(k,:),.75) + ...
                                                       w*(quantile(curData(k,:),.75)-quantile(curData(k,:),.25))) |...
                                        (curData(k,:) < quantile(curData(k,:),.25) - ...
                                                       w*(quantile(curData(k,:),.75)-quantile(curData(k,:),.25)));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

            % insert a break in the middle
            curData_outliers = [curData_outliers(1:length(curData_outliers)/2,:);...
                       nan(2,size(curData_outliers,2));...
                       curData_outliers((length(curData_outliers)/2+1):end,:)];


            hold on
            for k = 1:length(TimePoints)
                %plot(1:length(curData_outliers(:,k)), curData_outliers(:,k), 'Color', greyColor(k,:), 'LineWidth', 2)
                bar(1:length(curData_outliers(:,k)), curData_outliers(:,k), 'FaceColor', greyColor(k,:))
            end
            ylabel('Drug Concentration', 'FontSize', 12)

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
            % insert a break in the middle
            plotData = [plotData(1:length(plotData)/2,:);...
                       nan(2,size(plotData,2));...
                       plotData((length(plotData)/2+1):end,:)];



            plot(plotData, 'Color', [84 134 186]/256, 'LineWidth', 2)
            ylabel('Mean Z-score across time points', 'FontSize', 12)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot outlier line and detect candidate pools
            % plot outliers as q1-w*(q3-q1) for drugs
            % plot outliers as q3+w*(q3-q1) for drug compounds
            if plotDrugsFlag
                quantLine = quantile(plotData,.25)-w*(quantile(plotData,.75)-quantile(plotData,.25));
                plot([1 length(curData_outliers(:,k))], [quantLine quantLine], 'k--')
                candidatePools = strjoin(arrayfun(@(x) num2str(x),find(plotData<quantLine), 'unif',0),',');
            else
                w=3;
                quantLine = quantile(plotData,.75)+w*(quantile(plotData,.75)-quantile(plotData,.25));
                plot([1 length(curData_outliers(:,k))], [quantLine quantLine], 'k--')
                candidatePools = strjoin(arrayfun(@(x) num2str(x),find(plotData>quantLine), 'unif',0),',');
            end

            xlim([0 length(curData_outliers(:,k))+1])
            legend(legendarray,TimePoints, 'FontSize', 12)
            xlabel(myxlabel, 'FontSize', 12 )

            % get y limits and make 0 at mean of curData
            curylim = ylim();
            if abs(curylim(2)/curylim(1)) < curData_ylim
                curylim(2) = (abs(curylim(1))*curData_ylim);
            else
                curylim(1) = (-abs(curylim(2)/curData_ylim));
            end
            ylim(curylim)

            title([{i}, compoundNames{i}, {candidatePools}],'FontSize', 12)

            %     set(gca,'FontSize',6)      
        end  
    end
    filename = strrep(outfilename_Fig3D_heatmap_GoF, '.pdf', ['_part',num2str(spi),'_1.pdf']);
    print(fig, '-painters', '-dpdf', '-r600', '-bestfit', filename)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make heat map with the layout of the original plates => save all to PDF
    cmaplength = 40;
    fsize = 6;
    mycmap = ([linspace(8, 256,cmaplength);...
              linspace(81, 256,cmaplength);...
              linspace(156, 256,cmaplength)]/256)';
    mycmap = flipud(mycmap);

    i = compoundsIDX;
    Time = {'T24'};%{'T12'};

    curData = zeros(size(RawData.(TimePoints{k}).IntensitiesCorrected,1), length(TimePoints));
    for k = 1:length(TimePoints)
        curData(:,k) = RawData.(TimePoints{k}).IntensitiesCorrected(:,i);
    end
    plotData = zscore(curData);
    plotData = mean(plotData,2);


    myfig = figure('units','normalized','outerposition',[0 0 1 1]);

    mytitle = RawData.(Time{1}).Compounds(:,i);
    mytitle = strrep(mytitle, 'Results', '');

    legendarray = [];
    reshapedMatrix = reshapePlateIntensities(plotData,...RawData.(Time{1}).IntensitiesCorrected(:,i),...
                                     RawData.(Time{1}).SampleNames);
    h(1) = subplot(2,2,1);
    imagesc(reshapedMatrix(1:16, 1:24))
    set(gca, 'XTick', 1:24)
    set(gca, 'YTick', 1:16)
    set(gca, 'YTickLabel', arrayfun(@(x) char(x), 'A'+[0:15]'))
    caxis([ceil(min(min(reshapedMatrix))), floor(max(max(reshapedMatrix)))]) %color axis
    %caxis([floor(min(min(reshapedMatrix))), round(max(max(reshapedMatrix)))]) %color axis
    colorbar
    set(gca, 'fontSize', fsize)

    h(2) = subplot(2,2,2);
    imagesc(reshapedMatrix(1:16, 25:end))
    set(gca, 'XTick', 1:24)
    set(gca, 'YTick', 1:16)
    set(gca, 'YTickLabel', arrayfun(@(x) char(x), 'A'+[0:15]'))
    caxis([ceil(min(min(reshapedMatrix))), floor(max(max(reshapedMatrix)))]) %color axis
    %caxis([floor(min(min(reshapedMatrix))), round(max(max(reshapedMatrix)))]) %color axis
    colorbar
    set(gca, 'fontSize', fsize)

    h(3) = subplot(2,2,3);
    imagesc(reshapedMatrix(17:end, 1:24))
    set(gca, 'XTick', 1:24)
    set(gca, 'YTick', 1:16)
    set(gca, 'YTickLabel', arrayfun(@(x) char(x), 'A'+[0:15]'))
    caxis([ceil(min(min(reshapedMatrix))), floor(max(max(reshapedMatrix)))]) %color axis
    %caxis([floor(min(min(reshapedMatrix))), round(max(max(reshapedMatrix)))]) %color axis
    colorbar
    set(gca, 'fontSize', fsize)

    h(4) = subplot(2,2,4);
    imagesc(reshapedMatrix(17:end, 25:end))
    set(gca, 'XTick', 1:24)
    set(gca, 'YTick', 1:16)
    set(gca, 'YTickLabel', arrayfun(@(x) char(x), 'A'+[0:15]'))
    colormap(mycmap)
    caxis([ceil(min(min(reshapedMatrix))), floor(max(max(reshapedMatrix)))]) %color axis
    %caxis([floor(min(min(reshapedMatrix))), round(max(max(reshapedMatrix)))]) %color axis
    colorbar
    set(gca, 'fontSize', fsize)

    suptitle([[AssayedOD{1},', Time ' Time{1}], mytitle])

    p1 = get(h(1),'Position');
    p2 = get(h(2),'Position');
    p3 = get(h(3),'Position');
    p4 = get(h(4),'Position');
    % set right limit to match first subplot
    p1(3) = p1(3)*1.3;
    p2(3) = p1(3);
    p3(3) = p1(3);
    p4(3) = p1(3);
    % update positions of the subplots
    set(h(1),'Position',p1);
    set(h(2),'Position',p2);
    set(h(3),'Position',p3);
    set(h(4),'Position',p4);

    filename = strrep(outfilename_Fig3D_heatmap_GoF, '.pdf', ['_part',num2str(spi),'_2.pdf']);
    print(myfig, '-painters', '-dpdf', '-r600', '-bestfit', filename)
end

    