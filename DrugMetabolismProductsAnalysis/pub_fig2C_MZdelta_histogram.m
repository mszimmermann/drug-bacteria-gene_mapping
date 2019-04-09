%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_fig2C_MZdelta_histogram.m
% read drug-metabolite candidate table and plot a histogram with 
% mass difference between metabolites and drugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % takes as input outfile_tableS5_metabolite_filtering
% outfile_fig2c_metabolites_drug_mass_deltas

drugMetTable = readtable(outfile_tableS5_metabolite_filtering);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drugmets_drugmass_delta_trunc = drugMetTable.DrugMassDeltaSmoothed;

drugmets_drugmass_delta_trunc(isnan(drugmets_drugmass_delta_trunc)) = [];
% calculate number of occurances per smooth mass
drugmets_drugmass_delta_trunc_unique = unique(drugmets_drugmass_delta_trunc);
drugmets_drugmass_delta_trunc_unique_counts = zeros(size(drugmets_drugmass_delta_trunc_unique));
for i=1:length(drugmets_drugmass_delta_trunc_unique)
    drugmets_drugmass_delta_trunc_unique_counts(i) = nnz(drugmets_drugmass_delta_trunc==...
        drugmets_drugmass_delta_trunc_unique(i) &...
                                  (drugMetTable.GoodFilter==1) &...
                                  (drugMetTable.DrugMassFlag==0) &...
                                  (drugMetTable.NumberOfIncreasedT12vsT0 > 0) &...
                                  (drugMetTable.DrugConsumedFlag > 0) &...
                                  (drugMetTable.NumberOfDrugs == 1));
end

% remove 0 mz delta
drugmets_drugmass_delta_trunc_unique_counts(drugmets_drugmass_delta_trunc_unique==0) = 0;
delta_frequencies = drugmets_drugmass_delta_trunc_unique_counts;

plot_delta = drugmets_drugmass_delta_trunc_unique;
plot_counts = drugmets_drugmass_delta_trunc_unique_counts;
plot_delta = [plot_delta;...
              drugmets_drugmass_delta_trunc_unique(drugmets_drugmass_delta_trunc_unique_counts<2)+1; ...
              drugmets_drugmass_delta_trunc_unique(drugmets_drugmass_delta_trunc_unique_counts<2)-1];
plot_counts = [plot_counts;...
               zeros(2*nnz(drugmets_drugmass_delta_trunc_unique_counts<2),1)];
[plot_delta, sortidx] = sort(plot_delta);
plot_counts = plot_counts(sortidx);

figure;
plot(plot_delta,plot_counts)
xlabel('Mass delta between metabolite and the drug')
ylabel('Number of drug-metabolite pairs')
% text labels of the most frequent mass shifts
delta_frequencies_max = find(plot_counts>4);
delta_frequencies_ycoord = 5:((max(delta_frequencies)-5)/(nnz(delta_frequencies_max)-1)):max(delta_frequencies);
delta_frequencies_xcoord = zeros(length(delta_frequencies_max),1);
delta_frequencies_text = cell(length(delta_frequencies_max),1);

for i=1:length(delta_frequencies_max)
    curx = plot_delta(delta_frequencies_max(i));
    curxtext = num2str(curx);
    curx(curx>0) = curx(curx>0)+5;
    curx(curx<0) = curx(curx<0)-50;
    
    delta_frequencies_xcoord(i) = curx;
    delta_frequencies_text{i} = curxtext;
end
[~, idx] = sort( plot_counts(delta_frequencies_max) );
delta_frequencies_xcoord = delta_frequencies_xcoord(idx);
delta_frequencies_text = delta_frequencies_text(idx);
% add text with most frequent deltas
text(delta_frequencies_xcoord,delta_frequencies_ycoord,delta_frequencies_text)
ylim([0 max(delta_frequencies)+1])
xlim([-600 600])
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    outfile_fig2c_metabolites_drug_mass_deltas)
