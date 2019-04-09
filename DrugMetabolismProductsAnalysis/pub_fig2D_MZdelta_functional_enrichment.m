%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_fig2D_MZdelta_functional_enrichment.m
% read drug-metabolite candidate table and perform functional group enrichment
% of drugs with the frequent MZ deltas to metabolites
% save chemical functional group enrichment analysis to table
% save clustergram of mzdeltas and functional groups to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % takes as input infile_drug_experiment
% % takes as input outfile_tableS5_metabolite_filtering
% outfile_fig2d_clustergram_MZdelta_FGenrichment
% outfile_TableS4_mzdelta_chemical_functional_group_enrichment

% read functional frops of the selected drugs in experiment
drugTable = readtable(infile_drug_experiment);
% get column names
columnNames = drugTable.Properties.VariableNames;
drugNames = drugTable.MOLENAME;
drugNames = cellfun(@(x) strtrim(x), drugNames, 'unif',0);
column_fg = find(cellfun(@(x) contains(x, 'Number_of_'), columnNames));
drugFGnames = columnNames(column_fg);
drugFG = drugTable{:, column_fg};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load mass deltas from file
drugMetTable = readtable(outfile_tableS5_metabolite_filtering);

drugmets_drugmass_delta = drugMetTable.DrugMassDeltaSmoothed;

drugmets_drugmass_delta(isnan(drugmets_drugmass_delta)) = [];
% calculate number of occurances per smooth mass
drugmets_drugmass_delta_unique = unique(drugmets_drugmass_delta);
drugmets_drugmass_delta_unique_counts = zeros(size(drugmets_drugmass_delta_unique));
for i=1:length(drugmets_drugmass_delta_unique)
    drugmets_drugmass_delta_unique_counts(i) = nnz(drugmets_drugmass_delta==...
                                    drugmets_drugmass_delta_unique(i) &...
                                  (drugMetTable.GoodFilter==1) &...
                                  (drugMetTable.DrugMassFlag==0) &...
                                  (drugMetTable.NumberOfIncreasedT12vsT0 > 0) &...
                                  (drugMetTable.DrugConsumedFlag > 0) &...
                                  (drugMetTable.NumberOfDrugs == 1));
end

% remove 0 mz delta
drugmets_drugmass_delta_unique_counts(drugmets_drugmass_delta_unique==0) = 0;

% sort deltaMZ according to number of drugs
[drugmets_drugmass_delta_unique_counts, sortidx] = ...
    sort(drugmets_drugmass_delta_unique_counts, 'descend');
drugmets_drugmass_delta_unique = drugmets_drugmass_delta_unique(sortidx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform functional group enrichments for deltas with drugs >=3
% create table with FG percentage
delta_clusters = find(drugmets_drugmass_delta_unique_counts>=2 &...
                     abs(drugmets_drugmass_delta_unique)>1);
groupEnrichments_cell = cell(size(delta_clusters));
groupEnrichment_fraction = zeros(length(drugFGnames), length(delta_clusters));
groupEnrichment_fdr = zeros(length(drugFGnames), length(delta_clusters));
changingDrugsEnrichments = cell(size(delta_clusters));
changingDrugsEnrichments_num = zeros(size(delta_clusters));

for i = 1:length(changingDrugsEnrichments)
    changingDrugsEnrichments{i} = num2str(drugmets_drugmass_delta_unique(delta_clusters(i)));
    clustDrugNames = unique(drugMetTable.ParentDrug(drugmets_drugmass_delta == drugmets_drugmass_delta_unique(delta_clusters(i)) &...
                                  (drugMetTable.GoodFilter==1) &...
                                  (drugMetTable.DrugMassFlag==0) &...
                                  (drugMetTable.NumberOfIncreasedT12vsT0 > 0) &...
                                  (drugMetTable.DrugConsumedFlag > 0) &...
                                  (drugMetTable.NumberOfDrugs == 1)));
    clustDrugNames = cellfun(@(x) strtrim(x), clustDrugNames, 'unif', 0);
    changingDrugs = ismember(drugNames, clustDrugNames);
    
    groupEnrichments = zeros(length(drugFGnames),6);
    for idx = 1:length(drugFGnames)
        nGroupCountainingMolecules = nnz(drugFG(:,idx));% this is K - number of molecules containing group
        % this is k - number of molecules containing the group which are also in the changing group
        nChangingGroup = nnz(changingDrugs>0 & drugFG(:,idx)>0);
        % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
        % where N - total number of genes, n - size of the group
        nMolecules = length(drugNames);
        nChanging = nnz(changingDrugs);
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
    groupEnrichment_fraction(:,i) = groupEnrichments(:,3)./groupEnrichments(:,6);
    groupEnrichment_fdr(:,i) = groupEnrichments(:,1);
    changingDrugsEnrichments_num(i) = nChanging;
end

fid = fopen(outfile_TableS4_mzdelta_chemical_functional_group_enrichment,'w');
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
selectFG = ( sum(groupEnrichment_fraction>=0.5,2)>0) &...
           ( sum(groupEnrichment_fdr<=0.05,2)>0 ) ;
selectMZdelta = drugmets_drugmass_delta_unique_counts(delta_clusters)>3 ;
selectMZdelta = changingDrugsEnrichments_num>3;

curFGnames = drugFGnames;
curFGnames = cellfun(@(x) strrep(x,'Number_of_', ''), curFGnames, 'unif', 0);
curFGnames = cellfun(@(x) strrep(x,'_', ' '), curFGnames, 'unif', 0);

mycmap = [0 0 0;...
          0 0 0;...
          0 0 0;...
          0 0 0;...
          0 0 0;...
          0.5 0.5 0.5;...
          0.5 0.5 0.5;...
          0.75 0.75 0.75;...
          0.75 0.75 0.75;...
          1 1 1];

clustMatrix = (groupEnrichment_fraction(selectFG,selectMZdelta).*...
              (groupEnrichment_fraction(selectFG,selectMZdelta)>=0.5))';
% clustMatrix = (groupEnrichment_fraction(selectFG,selectMZdelta).*...
%               (groupEnrichment_fdr(selectFG,selectMZdelta)<=0.05))';
clustRows = changingDrugsEnrichments(selectMZdelta);
clustCols = curFGnames(selectFG);

removeCol = sum(clustMatrix>=0.5)==0;
clustCols(removeCol) = [];
clustMatrix(:,removeCol)=[];

clustdist = 'euclidean';%'cityblock';%'correlation';%'spearman';%'hamming';%
cgo=clustergram(clustMatrix,...
                 'DisplayRange', 1,...
                 'RowLabels', clustRows,...
                 'ColumnLabels', clustCols,...
                 'Symmetric', 0,...
                 'ColumnPdist', clustdist, 'RowPdist', clustdist,...
                 'Colormap', mycmap);%, 'DisplayRatio', 0.6);

                         
                         
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 6)
% insert colorbar
cbButton = findall(gcf,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback');
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
cb  = findobj(gcf,'Tag','HeatMapColorbar');
cb.Label.String = 'Fraction of drugs having the functional group';
 
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            outfile_fig2d_clustergram_MZdelta_FGenrichment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

