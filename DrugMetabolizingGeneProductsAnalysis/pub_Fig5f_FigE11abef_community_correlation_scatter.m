%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline pub_Fig5f_FigE11abef_community_correlation_scatter
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
% outfile_Fig5f_community_correlation_scatter
% intensityNoise = 5000; % noise intensity level 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read gene_drug_metabolite table
gene_drug_table = readtable(infile_gene_drug_matrix);
% extract information on gene_drugs (find parent columns in the matrix)
column_types = varfun(@class,gene_drug_table,'OutputFormat','cell');
column_parent = zeros(size(column_types));
for i=1:length(column_types)
    if isequal(column_types{i}, 'cell')
        column_parent(i) = contains(lower(gene_drug_table{1,i}),'parent');
    end
end
% create list of gene-drug pairs
gene_drug_list = cell(size(gene_drug_table,1)*size(gene_drug_table,2),2);
idx = 1;
for i=2:size(gene_drug_table,1)
    curdrugmat = cellfun(@(x) str2double(x), gene_drug_table{i,column_parent==1});
    curdrugs = gene_drug_table.Properties.VariableNames(column_parent==1);
    curdrugs = curdrugs(curdrugmat==1);
    gene_drug_list(idx:idx+length(curdrugs)-1,:) = [repmat(gene_drug_table.Gene(i),size(curdrugs)); curdrugs]';
    idx = idx+length(curdrugs);
end
gene_drug_list(idx:end,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read community CFU and qPCR info
community_cfu_table = readtable(infile_community_CFU_qPCR);
% read community OTU info
community_OTU_table = readtable(infile_community_OTU);
% get communty OTU names and select which are species levels
community_OTU_names = community_OTU_table.ID';
community_OTU_names_species = cellfun(@(x) contains(x, '|s') & ~contains(x, '|t'),  community_OTU_names);
% read community ShortBRED quantified proteins table
community_protein_table = readtable(infile_community_ShortBRED_proteins);
% read community drug metabolism table
community_drug_metabolism_table = readtable(outfile_table_community_drug_metabolism);
% get the slope of the drug and metabolite sfrom the table
column_slope = cellfun(@(x) contains(lower(x),'slope'), community_drug_metabolism_table.Properties.VariableNames);
column_slope = find(column_slope);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read clustergram table and get drug metabolizing species
drugAnalysisClustergram = readtable(outfile_table_FCclustergram,...
                                    'HeaderLines', 1);
FCcolumns = cellfun(@(x) contains(x, 'FC') &...
                         ~contains(x, 'FCSTD') &...
                         ~contains(x, 'DrugAdaptive'),...
                         drugAnalysisClustergram.Properties.VariableNames);

clustergramSpecies = cellfun(@(x) x(strfind(x, 'FC')+2:end),...
                                  drugAnalysisClustergram.Properties.VariableNames(FCcolumns), 'unif', 0)';
clustDrugNames = drugAnalysisClustergram.DrugName;
% remove spaces to match other data drug names
clustDrugNames = cellfun(@(x) strrep(x, ' ', ''), clustDrugNames, 'unif', 0);
drugFC12to0_t0toCTRL_combined = drugAnalysisClustergram{:, FCcolumns};
drugFCadaptive = repmat(drugAnalysisClustergram.DrugAdaptiveFCThreshold_, 1,size(drugFC12to0_t0toCTRL_combined,2));
drugFC12to0_metabolism = (drugFC12to0_t0toCTRL_combined <= pThreshold) &...
                         (drugFC12to0_t0toCTRL_combined <= -drugFCadaptive);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intersect community names across OTU, protein, cfu and drug metabolism
% tables
[commonCommunities, idxMET, idxCFU] = intersect(community_drug_metabolism_table.Community, community_cfu_table.Community, 'stable');
[~, ~, idxOTU] = intersect(commonCommunities, community_OTU_table.Properties.VariableNames);
[~, ~, idxProt] = intersect(commonCommunities, community_protein_table.Properties.VariableNames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correlation table
for i = 1:length(column_slope)
        curDrugName = community_drug_metabolism_table.Properties.VariableNames{column_slope(i)};
        % remove Slope_ and leave only drug name 
        curDrugName = strsplit(curDrugName,'_');
        curDrugName = curDrugName{2};
        
        metSpecies = clustergramSpecies(drugFC12to0_metabolism(ismember(clustDrugNames,...
                                  upper(curDrugName)),:)~=0);
        % remove strain IDs
        for si = 1:length(metSpecies)
            curspace = find(isstrprop(metSpecies{si},'upper') | isstrprop(metSpecies{si},'digit'));
            if length(curspace)>2
                %curspace = strfind(metSpecies{si}, ' ');
                metSpecies{si} = metSpecies{si}(1:curspace(3)-1);
            end
            if length(curspace)>1
                metSpecies{si} = [metSpecies{si}(1:curspace(2)-1)...
                                  ' '...
                                  lower(metSpecies{si}(curspace(2):end))];
            end
        end
        % remove duplicate species
        metSpecies = unique(metSpecies);
        %add _ to match OTU names
        metSpecies = cellfun(@(x) strrep(x, ' ', '_'), metSpecies, 'unif', 0);
        MVOTUmetSpecies = zeros(length(community_OTU_names),1);
        for si = 1:length(metSpecies)
            MVOTUmetSpecies = MVOTUmetSpecies+...
                              (cellfun(@(x) contains(x,metSpecies{si}), community_OTU_names) &...
                                           community_OTU_names_species)';
        end
        MVOTUmetSpecies = community_OTU_names(MVOTUmetSpecies~=0)';
        MVOTUmetSpecies_genus = unique(cellfun(@(x) x(1:strfind(x, '|s')-1), MVOTUmetSpecies, 'unif', 0));
        MVOTUmetSpecies_phylum = unique(cellfun(@(x) x(1:strfind(x, '|c')-1), MVOTUmetSpecies, 'unif', 0));
        % find indeces of the relevant OTUS
        [~, MVOTUmetSpecies] = intersect(community_OTU_names, MVOTUmetSpecies);
        [~, MVOTUmetSpecies_genus] = intersect(community_OTU_names,MVOTUmetSpecies_genus);
        [~, MVOTUmetSpecies_phylum] = intersect(community_OTU_names,MVOTUmetSpecies_phylum);
        
        % get proteins that metabolize the drug
        curDrugName = strrep(curDrugName, 'HYDROCHLORIDE','');
        curGeneNames = gene_drug_list(ismember(upper(gene_drug_list(:,2)), curDrugName));
        curGeneNames = unique(curGeneNames);
        [~, curProteinIDX] = intersect(community_protein_table.QuantifiedProtein, curGeneNames); 
        
        % compile matrix of features relevant for the drug: 
        % CFU, OTU (Phylum, genus, species), Proteins
        if ~isempty(curProteinIDX)
            curFeaturesMatrix = [community_cfu_table.CFUMean(idxCFU)...
                                 community_OTU_table{MVOTUmetSpecies_phylum,idxOTU}'...
                                 community_OTU_table{MVOTUmetSpecies_genus,idxOTU}'...
                                 community_OTU_table{MVOTUmetSpecies,idxOTU}'...
                                 community_protein_table{curProteinIDX, idxProt}'];
            curFeatures = ['Community CFU'...
                           cellfun(@(x) x(strfind(x, '|p'):end), community_OTU_names(MVOTUmetSpecies_phylum), 'unif',0)...
                           cellfun(@(x) x(strfind(x, '|g'):end), community_OTU_names(MVOTUmetSpecies_genus), 'unif',0)...
                           cellfun(@(x) x(strfind(x, '|s'):end), community_OTU_names(MVOTUmetSpecies), 'unif',0)...
                           community_protein_table.QuantifiedProtein{curProteinIDX}];


            curSlope = community_drug_metabolism_table{idxMET, column_slope(i)};

            [corrRvectorSlope, corrPvectorSlope] = corr(curSlope, curFeaturesMatrix);
            % round to 2 decimals
            corrRvectorSlope = round(corrRvectorSlope*100)/100;
            corrPvectorSlope = round(corrPvectorSlope*100)/100;

            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            set(groot, 'DefaultTextInterpreter', 'none')
            spX = 4;
            spY = 5;
            spidx = 1;
            if ~isempty(corrRvectorSlope)
                for corr_i = 1:length(corrRvectorSlope)
                    if spidx>spX*spY
                        suptitle(strrep(community_drug_metabolism_table.Properties.VariableNames{column_slope(i)}, 'Slope_', ''))
                        orient landscape
                        print(fig, '-painters', '-dpsc', '-r600', '-append', '-bestfit',...
                                 outfile_Fig5f_community_correlation_scatter);
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
                    title(sprintf('PCC=%.2f, p=%.2f', corrRvectorSlope(corr_i),...
                                                      corrPvectorSlope(corr_i)),...
                                                      'fontSize', 6)
                    spidx = spidx+1;
                end
                suptitle(strrep(community_drug_metabolism_table.Properties.VariableNames{column_slope(i)}, 'Slope_', ''))
            end
       
           orient landscape
           print(fig, '-painters', '-dpsc', '-r600', '-append', '-bestfit',...
                 outfile_Fig5f_community_correlation_scatter);
           close(fig);
       end  

end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

