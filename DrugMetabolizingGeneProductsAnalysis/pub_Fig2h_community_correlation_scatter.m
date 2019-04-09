%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Pipeline pub_Fig2h_community_correlation_scatter.m
% plot community drug and metabolite correlation scatterplots
% load data from file prepared in the script
% pub_analyze_community_drug_metabolism.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables;
% infile_gene_drug_matrix = 'tables13_gene_drug_matrix.csv' ;
% infile_community_CFU_qPCR = 'table_community_CFU_and_qPCR.csv';
% infile_community_OTU = 'table_community_OTU.csv';
% infile_community_ShortBRED_proteins = 'table_community_ShortBRED_proteins.csv';
% outfile_table_community_drug_metabolism = 'table_community_drug_metabolism.csv'
% intensityNoise = 5000; % noise intensity level 
% outfile_Fig2h_FigE4f_dexamethasone_scatter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read community CFU and qPCR info
community_cfu_table = readtable(infile_community_CFU_qPCR);
% read community drug metabolism table
community_drug_metabolism_table = readtable(outfile_table_community_drug_metabolism);
% get the slope of the drug and metabolite sfrom the table
column_slope = cellfun(@(x) contains(lower(x),'slope') &...
                            contains(lower(x), 'dexamethasone'), community_drug_metabolism_table.Properties.VariableNames);
column_slope = find(column_slope);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intersect community names across OTU, protein, cfu and drug metabolism
% tables
[commonCommunities, idxMET, idxCFU] = intersect(community_drug_metabolism_table.Community, community_cfu_table.Community, 'stable');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correlation table
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spX = 2;
spY = 2;
spidx = 1;
for i = 1:length(column_slope)
        curDrugName = community_drug_metabolism_table.Properties.VariableNames{column_slope(i)};
        % remove Slope_ and leave only drug name 
        curDrugName = strsplit(curDrugName,'_');
        curDrugName = curDrugName{2};
        
        curFeaturesMatrix = [community_cfu_table.CFUMean(idxCFU)...
                             community_cfu_table.Cscindens16S_qPCR_(idxCFU)];
        curFeatures = {'Community CFU'...
                       'S.scindens 16S (qPCR)'};

        curSlope = community_drug_metabolism_table{idxMET, column_slope(i)};

        [corrRvectorSlope, corrPvectorSlope] = corr(curSlope, curFeaturesMatrix);
        % round to 2 decimals
        corrRvectorSlope = round(corrRvectorSlope*100)/100;
        corrPvectorSlope = round(corrPvectorSlope*100)/100;

            
        for corr_i = 1:length(corrRvectorSlope)
            if spidx>spX*spY
                suptitle(strrep(community_drug_metabolism_table.Properties.VariableNames{column_slope(i)}, 'Slope_', ''))
                orient landscape
                print(fig, '-painters', '-dpsc', '-r600', '-append', '-bestfit',...
                         outfile_Fig5h_community_correlation_scatter);
                close(fig);  
                fig = figure('units','normalized','outerposition',[0 0 1 1]);
                spidx=1;
            end
            subplot(spX, spY, spidx)
            hold on
            scatter(curFeaturesMatrix(:,corr_i), curSlope, 20,...
                    mycolormap_communities(5:4+length(curSlope),:),...
                    'filled')
            lsline;
            axis square
            ylabel('Slope, a.u.', 'fontSize', 6)
            xlabel(curFeatures{corr_i},'interpreter', 'none', 'fontSize', 6)
            set(gca, 'fontSize', 6)
            title({ strrep(community_drug_metabolism_table.Properties.VariableNames{column_slope(i)}, 'Slope_', ''),...
                    sprintf('PCC=%.2f, p=%.2f', corrRvectorSlope(corr_i),...
                                                corrPvectorSlope(corr_i))},...
                                               'fontSize', 6,'interpreter', 'none')
            spidx = spidx+1;
        end
end   
orient landscape
print(fig, '-painters', '-dpdf', '-r600', '-append', '-bestfit',...
     outfile_Fig2h_FigE4f_dexamethasone_scatter);
close(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

