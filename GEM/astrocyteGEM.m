gtex_data = readtable('D:/GEM/mydata/mean.txt');
gtex_data(1:5,1:2)

% extract the tissue and gene names
data_struct.tissues = gtex_data.Properties.VariableNames(2:end)';  % sample (tissue) names
data_struct.genes = gtex_data.genes;  % gene names
data_struct.levels = table2array(gtex_data(:, 2:end));  % gene TPM values
data_struct.threshold = 1;
data_struct

%load GEM
load('Mouse-GEM.mat');  % loads model as a structure named "ihuman"
mouseGEM = addBoundaryMets(mouseGEM);
essentialTasks = parseTaskList('D:/GEM/metabolicTasks_Essential.txt')

% see what the other inputs mean by typing "help checkTasks"
checkTasks(mouseGEM, [], true, false, false, essentialTasks);

refModel = mouseGEM;          % the reference model from which the GEM will be extracted
tissue = 'astrocyte';           % must match the tissue name in data_struct.tissues
celltype = [];              % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];               % data structure containing protein abundance information (not used here)
arrayData = data_struct;    % data structure with gene (RNA) abundance information
metabolomicsData = [];      % list of metabolite names if metabolomic data is available
removeGenes = true;         % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];              % we already loaded the task file, so this input is not required
useScoresForTasks = true;   % (default) use expression data to decide which reactions to keep
printReport = true;         % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params = [];                % additional optimization parameters for the INIT algorithm
paramsFT = [];              % additional optimization parameters for the fit-tasks algorithm


astrocyteGEM = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

checkTasks(astrocyteGEM, [], true, false, false, essentialTasks);

save('D:/GEM/mydata/astrocyteGEM.mat', 'astrocyteGEM');

