%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED2B_drugs_functional_groups_and_enrichments
% read degraded drug fold change table and 
% perform chemical functional group enrichment analysis of degraded drugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables;
% takes as input infile_drug_experiment
% takes as input infile_figED1_DrugBankDrugInfo
% takes as input outfile_table_FCclustergram
%outfile_tableS4_functional_group_enrichments = ['Output' filesep 'TableS4_functional_group_enrichment_changing_drugs_all.csv'];
%outfile_figED2B_drug_functional_groups = ['Output' filesep 'FigE2_B_bar_functional_groups_metabolized_drugs.pdf'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExperimentParameters = import_experiment_parameters(infile_drug_experiment);
ExperimentParameters.DrugNames = cellfun(@(x) strtrim(x), ExperimentParameters.DrugNames, 'unif', 0);
if is_available_DrugBankFile 
    % get drugbank FG information
    drugBankTable = readtable(infile_figED1_DrugBankDrugInfo);
    % get column names
    columnNames = drugBankTable.Properties.VariableNames;
    column_fg = find(cellfun(@(x) contains(x, 'Number_of_'), columnNames));
    drugbankFGnames = columnNames(column_fg);
    drugbankFG = drugBankTable{:, column_fg};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read functional frops of the selected drugs in experiment
drugTable = readtable(infile_drug_experiment);
% get column names
columnNames = drugTable.Properties.VariableNames;
drugNames = drugTable.MOLENAME;
drugNames = cellfun(@(x) strtrim(x), drugNames, 'unif',0);
column_fg = find(cellfun(@(x) contains(x, 'Number_of_'), columnNames));
drugFGnames = columnNames(column_fg);
drugFG = drugTable{:, column_fg};

[drugNames, idxSelected, idxAll] = intersect(ExperimentParameters.DrugNames, drugNames, 'stable');
drugFG = drugFG(idxAll,:);

% read clustergram info
drugAnalysisClustergram = readtable(outfile_table_FCclustergram,...
                                'HeaderLines',1); %skip first line);
% get the names of consumed drugs from the clustergram
clustDrugNames = drugAnalysisClustergram.DrugName;
clustDrugNames = cellfun(@(x) strtrim(x), clustDrugNames, 'unif',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform functional group enrichment analysis
changingDrugsEnrichments = {'Metabolized', 'Cluster1', 'Cluster2'};

clustDrugNames1 = {'NORETHINDRONE ACETATE'
    'VILAZODONE'
    'FAMCICLOVIR'
    'ROXATIDINE ACETATE'
    'RACECADOTRIL'};
clustDrugNames2 = {'SULFASALAZINE'
    'PHENAZOPYRIDINE'
    'NITRENDIPINE'
    'TINIDAZOLE'
    'ENTACAPONE'};

changingDrugs{1} = ismember(drugNames, clustDrugNames);
changingDrugs{2} = ismember(drugNames, clustDrugNames1);
changingDrugs{3} = ismember(drugNames, clustDrugNames2);
groupEnrichments_cell = cell(size(changingDrugsEnrichments));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(changingDrugsEnrichments)
    groupEnrichments = zeros(length(drugFGnames),6);
    for idx = 1:length(drugFGnames)
        nGroupCountainingMolecules = nnz(drugFG(:,idx));% this is K - number of molecules containing group
        % this is k - number of molecules containing the group which are also in the changing group
        nChangingGroup = nnz(changingDrugs{i}>0 & drugFG(:,idx)>0);
        % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
        % where N - total number of genes, n - size of the group
        nMolecules = length(drugNames);
        nChanging = nnz(changingDrugs{i});
        score = sum( hygepdf(nChangingGroup:nGroupCountainingMolecules,...
                             nMolecules,...
                             nGroupCountainingMolecules,...
                             nChanging));
        groupEnrichments(idx,1) = score;
        groupEnrichments(idx,2) = score;
        groupEnrichments(idx,3:6) = [nChangingGroup,nMolecules, nGroupCountainingMolecules,nChanging];
    end
    groupEnrichments(:,2) = mafdr(groupEnrichments(:,2),'BHFDR', 1); 
    groupEnrichments_cell{i} = groupEnrichments;
end

fid = fopen(outfile_tableS4_functional_group_enrichments,'w');
fprintf(fid,'Functional group');
for i=1:length(changingDrugsEnrichments)
    fprintf(fid,';%s p-value;%s FDR;%s stat',changingDrugsEnrichments{i},...
        changingDrugsEnrichments{i}, changingDrugsEnrichments{i});
end
fprintf(fid,'\n');
for i=1:length(drugFGnames)
    fprintf(fid,'%s',drugFGnames{i});
    for j=1:length(groupEnrichments_cell)
        fprintf(fid,';%.4f;%.4f;(%d,%d,%d,%d)',groupEnrichments_cell{j}(i,:));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each functional group, plot how many are in cluster and how many
% are not
if is_available_DrugBankFile 
    drugbankFGsum = sum(drugbankFG>0);
    drugFGsum = sum(drugFG>0);
    curFGnames = drugFGnames;
    for i=1:length(curFGnames)
        drugbankidx = ismember(drugbankFGnames, curFGnames{i});
        curFGnames{i} = sprintf('%s (selected drugs %d/ DrugBank %d)', curFGnames{i}, drugFGsum(i), drugbankFGsum(drugbankidx));
    end
else
    drugFGsum = sum(drugFG>0);
    curFGnames = drugFGnames;
    for i=1:length(curFGnames)
        curFGnames{i} = sprintf('%s (selected drugs %d)', curFGnames{i}, drugFGsum(i));
    end
end

i = 1; % cluster "metabolized"
groupEnrichments = groupEnrichments_cell{i};
percentIn = groupEnrichments(:,3)./groupEnrichments(:,5);
[~, sortidx] = sort(percentIn);
sortidx(groupEnrichments(sortidx,5)==0)=[];
curFGnames = curFGnames(sortidx);
curFGnames = cellfun(@(x) strrep(x,'Number_of_', ''), curFGnames, 'unif', 0);
curFGnames = cellfun(@(x) strrep(x,'_', ' '), curFGnames, 'unif', 0);

fig = figure('units','normalized','outerposition',[0 0 1 1]);
barh(1:length(sortidx), groupEnrichments(sortidx,3), 'b')
hold on
barh(1:length(sortidx), groupEnrichments(sortidx,3)-groupEnrichments(sortidx,5),'r')
set(gca, 'Ytick', 1:length(sortidx))
set(gca, 'YTickLabel', curFGnames)
%set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(gca, 'Xtick', [-100 -50 0 75 150])
set(gca, 'XTickLabel', {-100 'Not metabolized', 0 'Metabolized' 150})
orient landscape
print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
            outfile_figED2B_drug_functional_groups)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
