function ExperimentParameters = import_experiment_parameters(filename)  
% import experiment parameters from filename:
% drug names, 
% drug pooling scheme, 
% drug MZ and RT
ExperimentTable = readtable(filename);
ExperimentParameters.DrugNames = ExperimentTable.MOLENAME;
drugPoolScheme = zeros(length(ExperimentParameters.DrugNames), 21);
for i=1:size(ExperimentTable,1)
    curpools = ExperimentTable.DrugPoolNumbers{i};
    curpools = strsplit(curpools);
    curpools = cellfun(@(x) str2double(x), curpools);
    drugPoolScheme(i, curpools) = 1;
end
ExperimentParameters.PoolingScheme = drugPoolScheme;
ExperimentParameters.DrugMolWeight = ExperimentTable.molWeight;
ExperimentParameters.DrugLogP = ExperimentTable.molLogP;
ExperimentParameters.DrugEstimatedColonConcentration = ExperimentTable.EstimatedColonConcentrationMaier2018uM;
ExperimentParameters.DrugMasses = ExperimentTable.TargetedMZ;
ExperimentParameters.TargetedRT = ExperimentTable.TargetedRT;