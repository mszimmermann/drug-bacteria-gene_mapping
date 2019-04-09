%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_fig2B_drug_metabolite_histogram.m
% read drug-metabolite candidate table and plot a histogram with 
% the number of metabolites per drug before and after filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % takes as input outfile_table_FCclustergram
% % takes as input outfile_tableS5_metabolite_filtering
% outfile_fig2b_metabolites_per_drug
% outfile_fig2b_metabolites_per_drug_zoomed

drugMetTable = readtable(outfile_tableS5_metabolite_filtering);
drugClustergram = readtable(outfile_table_FCclustergram, 'HeaderLines',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each drug, calculate how many metabolites occured at least in
% one species in its pool at t=0 and t=12 (1st and 2nd column)
drugmets_drugNames = drugMetTable.ParentDrug;
drugmets_drugNames_unique = unique(drugmets_drugNames);
% select only drugs that are metabolized
drugmets_drugNames_unique = intersect(drugmets_drugNames_unique, drugClustergram.DrugName);
changingMets_per_drug_filtered =  zeros(length(drugmets_drugNames_unique),1);
changingMets_per_drug_unfiltered =  zeros(length(drugmets_drugNames_unique),1);

% select only metabolites that pass filtering criteria
good_mets = ones(size(drugMetTable,1),1);        
for i=1:length(drugmets_drugNames_unique)
   changingMets_per_drug_unfiltered(i) = nnz(ismember(drugmets_drugNames(good_mets==1), drugmets_drugNames_unique{i}));
end

good_mets = drugMetTable.GoodFilter==1 &...
            drugMetTable.DrugMassFlag==0 & ...
            drugMetTable.NumberOfDrugs==1 &...
            drugMetTable.NumberOfIncreasedT12vsT0>0 &...
            drugMetTable.DrugConsumedFlag==1;        
for i=1:length(drugmets_drugNames_unique)
   changingMets_per_drug_filtered(i) = nnz(ismember(drugmets_drugNames(good_mets==1), drugmets_drugNames_unique{i}));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filtered and unfiltered in one plot
origEdges = unique([changingMets_per_drug_filtered;...
                    changingMets_per_drug_unfiltered]);
origHeight_unfiltered = histcounts(changingMets_per_drug_unfiltered,origEdges);
origHeight_filtered = histcounts(changingMets_per_drug_filtered,origEdges);

origMiddles = origEdges(1:end-1)+(origEdges(2:end)-origEdges(1:end-1))/2;
figure;
bar(origMiddles, origHeight_unfiltered, 1,...
    'FaceColor', [.5 .5 .5], 'EdgeColor', [.5 .5 .5]);
hold on
bar(origMiddles, origHeight_filtered, 1,...
    'FaceColor', [1 0 0], 'EdgeColor', [1 0 0]);

xlabel('Number of metabolites detected per drugs pool')
title('Number of potential drug metabolites per drug')
ylabel('Number of drugs')
orient landscape
legend('Unfiltered', 'Filtered')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
     outfile_fig2b_metabolites_per_drug)
% zoom to 0..100
xlim([0 100])
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    outfile_fig2b_metabolites_per_drug_zoomed);

