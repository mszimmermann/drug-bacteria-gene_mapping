%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_FigED1D_hist_intestinal_concentrations
% read experiment drug info with predicted intestinal concantrations from
% Maier, L., Pruteanu, M., Kuhn, M., Zeller, G., Telzerow, A., Anderson, E.E., 
% Brochado, A.R., Fernandez, K.C., Dose, H., Mori, H. and Patil, K.R., 2018.
% Extensive impact of non-antibiotic drugs on human gut bacteria.
% Nature, 555(7698), p.623.
% plot histogram of intestinal concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% % takes as input infile_drug_experiment
% outfile_figED1C_intconc = 'AnalysisAfterCleanup\OutputFiles\FigEV1C_hist_predicted_intestinal_concentrations.pdf';

%Import Experimental Parameters
ExperimentParameters = import_experiment_parameters(infile_drug_experiment);
displayConc = ExperimentParameters.DrugEstimatedColonConcentration;
displayConc(isnan(displayConc)) = [];

figure
hist(log10(displayConc))
xlim([-2 4])
ylim([0 20])
axis square
ylabel('Number of drugs')
xlabel('Estimated intestinal concentration, uM, log10')
text(-1,18,sprintf('Mean conc = %d', round(mean(displayConc))));
text(-1,16,sprintf('Median conc = %d', round(median(displayConc))));

title(sprintf('Predicted intestinal concentrations of selected drugs (%d)', length(displayConc)))

% save figure to file
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
     outfile_figED1D_intconc)
  


