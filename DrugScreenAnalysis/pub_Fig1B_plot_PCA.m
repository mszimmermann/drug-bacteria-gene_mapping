%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script pub_Fig1B_plot_PCA
% read drug principle component coordinates for drugBANK and PharmaCon libraries
% plot 3D PCA plot and save to pdf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Maria Zimmermann-Kogadeeva and Michael Zimmermann (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read global variables defining thresholds and file names
Drug_bacteria_gene_mapping_variables
% will use the following variables:
% infile_drug_experiment = 'Table2_DrugScreenInfo.csv';
% infile_fig1B_drugbankPCA = 'drugbank_pca.txt';
% infile_fig1B_pharmacomPCA = 'pharmacon_pca.txt';
% outfile_fig1B_PCA = 'Output\fig1B_plot_rdkit_drugbank_selected_drugs_PCA.pdf';

%Import Experimental Parameters
ExperimentParameters = import_experiment_parameters(infile_drug_experiment);

drugbank_pca = readtable(infile_fig1B_drugbankPCA);
drugbank_pca = [drugbank_pca.Var1 drugbank_pca.Var2 drugbank_pca.Var3];

pharmacon_pca = readtable(infile_fig1B_pharmacomPCA);
pharmacon_names = pharmacon_pca.Var1;
pharmacon_pca = [pharmacon_pca.Var2 pharmacon_pca.Var3 pharmacon_pca.Var4];

[~, idxpharm] = intersect(pharmacon_names, ExperimentParameters.DrugNames);

figure;
basecolor = [.5 .5 .5];
maincolor = [222 45 38]/256;
scatter1 = scatter3(drugbank_pca(:,1),drugbank_pca(:,2),drugbank_pca(:,3),...
                    'MarkerFaceColor',basecolor,'MarkerEdgeColor',basecolor); 
% Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
scatter1.MarkerFaceAlpha = .3;
scatter1.MarkerEdgeAlpha = .3;

hold on
scatter2 = scatter3(pharmacon_pca(idxpharm,1),pharmacon_pca(idxpharm,2),...
      pharmacon_pca(idxpharm,3), 's', 'MarkerEdgeColor', maincolor,...
      'MarkerFaceColor', maincolor);
grid on
xlim([min(drugbank_pca(:,1)) max(drugbank_pca(:,1))])
ylim([min(drugbank_pca(:,2)) max(drugbank_pca(:,2))])
zlim([min(drugbank_pca(:,3)) max(drugbank_pca(:,3))])
xlabel('Principal component 1')
ylabel('Principal component 2')
zlabel('Principal component 3')
legend({'Drug bank', 'Selected drugs'})
orient landscape
print(gcf, '-painters','-dpdf','-r600','-bestfit', ...
       outfile_fig1B_PCA);
   