%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_fig2F_dexamethasone_single_species_barplot
% plot dexamethasone drug and metabolite fold changes sorted by 
% species order in the clustergram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% takes as input outfile_table_FCclustergram
% takes as input outfile_tableS6_metaboliteIntensityFCData
% outfile_fig2f_dexamethasone_barplot
% massThreshold = 0.002;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load clustergram of changing drugs
drugClustergram = readtable(outfile_table_FCclustergram,...
                            'HeaderLines', 1);
drugAnalysisColumnNames = drugClustergram.Properties.VariableNames;
FCcolumns = cellfun(@(x) contains(x, 'FC') &...
                         ~contains(x, 'FCSTD') &...
                         ~contains(x, 'DrugAdaptive'),...
                         drugAnalysisColumnNames);
FCSTDcolumns = cellfun(@(x) contains(x, 'FCSTD'),...
                         drugAnalysisColumnNames);
PFDRcolumns = cellfun(@(x) contains(x, 'p_FDR'),...
                         drugAnalysisColumnNames);                     
drugFC12to0_t0toCTRL_combined = drugClustergram{:,FCcolumns}';
drugFC12to0_t0toCTRL_STDcombined = drugClustergram{:,FCSTDcolumns}';
drugP12to0_t0toCTRL_combined = drugClustergram{:,PFDRcolumns}';
drugFC12to0_Species = drugAnalysisColumnNames(FCcolumns);
% replace FC with empty
drugFC12to0_Species = cellfun(@(x) strrep(x, 'FC', ''), drugFC12to0_Species, 'unif',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metFoldChanges = readtable(outfile_tableS6_metaboliteIntensityFCData,...
                            'HeaderLines', 1);
FCcolumns = cellfun(@(x) contains(x, 'FCT_12hVsT_0h'),...
                         metFoldChanges.Properties.VariableNames);
pFDRcolumns = cellfun(@(x) contains(x, 'p_FDR_T_12hVsT_0h'),...
                         metFoldChanges.Properties.VariableNames);   
Int12columns = cellfun(@(x) contains(x, 'IntensityMeanT_12h'),...
                         metFoldChanges.Properties.VariableNames);
Int12STDcolumns = cellfun(@(x) contains(x, 'IntensitySTDT_12h'),...
                         metFoldChanges.Properties.VariableNames);

metFC12to0 = metFoldChanges{:,FCcolumns};
metFDR12to0 = metFoldChanges{:,pFDRcolumns};
metInt12 = metFoldChanges{:,Int12columns};
metInt12_STD = metFoldChanges{:,Int12STDcolumns};

metFCspecies = metFoldChanges.Properties.VariableNames(FCcolumns)';
% replace FC with empty
metFCspecies = cellfun(@(x) strrep(x, 'FCT_12hVsT_0h_log2_', ''), metFCspecies, 'unif',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intersect species names
[~, drugnameidx, metnameidx] = intersect(drugFC12to0_Species, metFCspecies);

% find drug and metabolite indeces
drug_idx = find(cellfun(@(x) contains(x, 'DEXAME'),drugClustergram.DrugName));
met_idx = find(cellfun(@(x) contains(x, 'DEXAME'), metFoldChanges.ParentDrug) &...
               arrayfun(@(x) abs(x-332.179)<massThreshold, metFoldChanges.MZ));

drugFC = drugFC12to0_t0toCTRL_combined(:,drug_idx);
drugSTD = drugFC12to0_t0toCTRL_STDcombined(:,drug_idx);

drugSTD = 100*drugSTD.*(2.^drugFC)*log(2);
drugFC = 100*(1-2.^drugFC);
    
figure
subplot(2,1,1)
bar(drugFC)
hold on
errorbar(1:length(drugFC12to0_Species),...
         drugFC,...
         drugSTD, 'k.')
ylim([-20 120])
xlim([0 77])
ylabel('Percent of drug consumed')
title(drugClustergram.DrugName(drug_idx))

subplot(2,1,2)
bar(drugnameidx, metInt12(met_idx,metnameidx))
hold on
errorbar(drugnameidx,...
         metInt12(met_idx,metnameidx),...
         metInt12_STD(met_idx,metnameidx), 'k.')
ylim([0 3*10^6])
xlim([0 77])
ylabel('Intensity at t=12h')
title(strcat(metFoldChanges.ParentDrug(met_idx),{' '},num2str(metFoldChanges.MZ(met_idx))))

orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            outfile_fig2f_dexamethasone_barplot);
        