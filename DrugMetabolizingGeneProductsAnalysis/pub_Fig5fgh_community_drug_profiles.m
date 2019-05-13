%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pipeline pub_Fig5fgh_community_drug_profiles
% plot community drug and metabolite profiles 
% load data from file prepared in the script
% pub_analyze_community_drug_metabolism.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables;
% outfile_table_community_drug_metabolism = 'table_community_drug_metabolism.csv'
% outfile_community_drug_met_profiles 
% intensityNoise = 5000; % noise intensity level 
% mycolormap_communities 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read community drug metabolism table
community_drug_metabolism_table = readtable(outfile_table_community_drug_metabolism);
% get the slope of the drug and metabolite sfrom the table
column_Mean = cellfun(@(x) contains(lower(x),'mean_'), community_drug_metabolism_table.Properties.VariableNames);
column_Mean = find(column_Mean);
community_time = community_drug_metabolism_table{1,column_Mean};
community_Mean = community_drug_metabolism_table{2:end,column_Mean};
community_compounds = community_drug_metabolism_table.Properties.VariableNames(column_Mean);
for i=1:length(community_compounds)
    cursplit = strsplit(community_compounds{i},'_');
    if mod(length(cursplit),2)
        community_compounds{i} = strjoin(cursplit(2:end-1),'_');
    else
        community_compounds{i} = strjoin(cursplit(2:end),'_');
    end
end
community_compounds_unique = unique(community_compounds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(community_compounds_unique)
    curidx = ismember(community_compounds, community_compounds_unique{i});
    fig = figure;
    hold on
    for j=1:size(community_Mean,1)
        plot(community_time(curidx), community_Mean(j,curidx),'LineWidth', 2,...
             'Color', mycolormap_communities(j,:));
    end
    title(community_compounds_unique{i}, 'interpreter', 'none')
    orient landscape
    print(fig, '-painters', '-dpsc', '-r600', '-append', '-bestfit',...
        outfile_community_drug_met_profiles);
    close(fig);  
end

