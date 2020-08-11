% Model validation: SNM 1.1rev4
% Created 14 Apr 2020 JR

%% Import simulation predictions (pred)
homedir = "data\";
valdir = "1_ModelValidation\";
plotdir = "plots\";
addpath(homedir)
addpath(strcat(homedir,valdir));
fprintf("Loading Data...\n")
filedir = "snm_1_1_rev4_validation_inputs.mat";
load(filedir)
filedir = "snm_1_1_rev4_validation_predictions.mat";
load(filedir)
filedir = "snm_1_1_rev4_validation_speciesNames.mat";
load(filedir)

modelpath = 'snm_1_1_rev4.xlsx';
sheetname = "reactions";
opts = detectImportOptions(modelpath,"Sheet",sheetname);
snm = readtable(modelpath, opts);
idx_rxns = find(contains(snm.module,"output"));
snm_rxns = split(snm.Rule(idx_rxns)," => ");
snm_outputs = unique(snm_rxns(:,end));
outputs_idx = zeros(length(snm_outputs),1);
% ID model outputs
for output = 1:length(snm_outputs)
    outputs_idx(output) = find(strcmp(speciesNames,snm_outputs(output)));
end

% create prediction tables
pred_full = array2table(predictions,'VariableNames',speciesNames.');
inputNames = cell2table(inputNames.','VariableNames',{'input'});
pred_full = [inputNames pred_full];
% subset for outputs only
outputs_pred = pred_full.Properties.VariableNames(outputs_idx+1);
pred = pred_full(:,outputs_idx+1);
% subset for intermediates
pred_med = pred_full;
pred_med(:,outputs_idx+1) = [];

%% Import SNM 1.0 predictions (pred_1_0)
% filedir = "snm_1_0_validation_predictions.mat";
% load(filedir)
% filedir = "snm_1_0_validation_speciesnames.mat";
% load(filedir)

% % add 'input' variable
% pred_1_0_full = array2table(predictions_SNM_1_0,'VariableNames',speciesNames_1_0.');
% pred_1_0_full = [inputNames pred_1_0_full];
% % subset for outputs only
% outputs_1_0_idx = [91,92,70,8,40,80,79,81,82,78,76,77,74,87,25];
% outputs_1_0_pred = pred_1_0_full.Properties.VariableNames(outputs_1_0_idx);
% pred_1_0 = pred_1_0_full(:,outputs_1_0_idx);
% % subset for intermediates
% pred_1_0_med = pred_1_0_full;
% pred_1_0_med(:,outputs_1_0_idx) = [];

%% Import experimental validation (expt)
filedir = "snm_1_1_validation.xlsx";
opts = detectImportOptions(filedir);
expt = readtable(filedir,opts);

% convert 'Measurement' strings to match prediction table (ie. integers)
expt.Measurement(strcmp(expt.Measurement,'Increase')) = {1};
expt.Measurement(strcmp(expt.Measurement,'Decrease')) = {-1};
expt.Measurement(strcmp(expt.Measurement,'No Change')) = {0};
expt.Measurement = cell2mat(expt.Measurement);
% convert 'Prediction' strings to match prediction table
expt.Prediction(strcmp(expt.Prediction,'Increase')) = {1};
expt.Prediction(strcmp(expt.Prediction,'Decrease')) = {-1};
expt.Prediction(strcmp(expt.Prediction,'No Change')) = {0};
expt.Prediction(strcmp(expt.Prediction,'NA')) = {0};
expt.Prediction = cell2mat(expt.Prediction);
expt.Properties.VariableNames(5) = {'SNM_1_0_1pct'};

% subset into in-out/in-med
expt_med = expt(strcmp(expt.in_out,"in-med"),:);
expt(strcmp(expt.in_out,"in-med"),:) = [];

%% Combine expt/pred tables (in-out):
% Convert pred table to long-form
pred.Input = pred_full.input;
pred_long = stack(pred,1:width(pred)-1);    % stack all vars except input
pred_long.Properties.VariableNames = {'Input','Output','SNM_1_1'};
pred_long.ID_pred = [1:height(pred_long)].';
pred_long.Output = cellstr(string(pred_long.Output));    % match class with expt
% join pred with expt based on matching fields: Input, Output
joined = innerjoin(expt,pred_long,'Keys',{'Input','Output'});

% % Repeat for pred_1_0
% pred_1_0.Input = pred_1_0_full.input;
% pred_1_0_long = stack(pred_1_0,1:width(pred_1_0)-1);    % stack all vars except input
% pred_1_0_long.Properties.VariableNames = {'Input','Output','SNM_1_0_5pct'};
% pred_1_0_long.ID_pred_1_0 = [1:height(pred_1_0_long)].';
% pred_1_0_long.Output = cellstr(string(pred_1_0_long.Output));    % match class with expt
% % Convert 'Output' strings to match expt table (CI, CIII)
% pred_1_0_long.Output(strcmp(pred_1_0_long.Output,'CI')) = {'proCI'};
% pred_1_0_long.Output(strcmp(pred_1_0_long.Output,'CIII')) = {'proCIII'};
% join pred with expt based on matching fields: Input, Output
% joined = innerjoin(joined,pred_1_0_long,'Keys',{'Input','Output'});

% Test for equivalence between measurements/predictions
% joined.Correct_1_0_1pct = (joined.Measurement == joined.SNM_1_0_1pct);
% joined.Correct_1_0_5pct = (joined.Measurement == joined.SNM_1_0_5pct);
joined.Correct_1_1 = (joined.Measurement == joined.SNM_1_1);

% Get percentage of correct predictions
s = summary(joined);
% pct_1_0_1pct = s.Correct_1_0_1pct.True / s.Correct_1_0_1pct.Size(1);
% pct_1_0_5pct = s.Correct_1_0_5pct.True / s.Correct_1_0_5pct.Size(1);
pct_1_1 = s.Correct_1_1.True / s.Correct_1_1.Size(1);
% fmt = 'SNM 1.0 Output Accuracy (0.01 Threshold): %2.3f \n';
% fprintf(fmt,pct_1_0_1pct)
% fmt = 'SNM 1.0 Output Accuracy (0.05 Threshold): %2.3f \n';
% fprintf(fmt,pct_1_0_5pct)
fmt = 'SNM 1.1 Output Accuracy: %2.3f \n';
fprintf(fmt,pct_1_1)

%% Combine expt/pred tables (in-med):
% Convert pred table to long-form
pred_med_long = stack(pred_med,2:width(pred_med));    % stack all vars except input
pred_med_long.Properties.VariableNames = {'Input','Output','SNM_1_1'};
pred_med_long.ID_pred = [1:height(pred_med_long)].';
pred_med_long.Output = cellstr(string(pred_med_long.Output));    % match class with expt
% join pred with expt based on matching fields: Input, Output
joined_med = innerjoin(expt_med,pred_med_long,'Keys',{'Input','Output'});

% % Repeat for pred_1_0
% pred_1_0_med_long = stack(pred_1_0_med,2:width(pred_1_0_med));    % stack all vars except input
% pred_1_0_med_long.Properties.VariableNames = {'Input','Output','SNM_1_0_5pct'};
% pred_1_0_med_long.ID_pred_1_0 = [1:height(pred_1_0_med_long)].';
% pred_1_0_med_long.Output = cellstr(string(pred_1_0_med_long.Output));    % match class with expt
% % Convert 'Output' strings to match expt table (CI, CIII)
% pred_1_0_med_long.Output(strcmp(pred_1_0_med_long.Output,'CI')) = {'proCI'};
% pred_1_0_med_long.Output(strcmp(pred_1_0_med_long.Output,'CIII')) = {'proCIII'};
% % join pred with expt based on matching fields: Input, Output
% joined_med = innerjoin(joined_med,pred_1_0_med_long,'Keys',{'Input','Output'});
% 
% % Test for equivalence between measurements/predictions
% joined_med.Correct_1_0_5pct = (joined_med.Measurement == joined_med.SNM_1_0_5pct);
joined_med.Correct_1_1 = (joined_med.Measurement == joined_med.SNM_1_1);

% Get percentage of correct predictions
s_med = summary(joined_med);
% pct_1_0_5pct = s_med.Correct_1_0_5pct.True / s_med.Correct_1_0_5pct.Size(1);
pct_1_1 = s_med.Correct_1_1.True / s_med.Correct_1_1.Size(1);
% fmt = 'SNM 1.0 Intermediate Accuracy (0.05 Threshold): %2.3f \n';
% fprintf(fmt,pct_1_0_5pct)
fmt = 'SNM 1.1 Intermediate Accuracy: %2.3f \n';
fprintf(fmt,pct_1_1)

%% Visualize predictions:
joined = sortrows(joined,{'Input','Output'});   %arrange rows for uniformity
joined_med = sortrows(joined_med,{'Input','Output'});
cmap = [0 0.502 0.752; 0.75 0.75 0.75; 1 0 0];
nmax = 11;  % maximum number of outputs  for all subsets
xdata = ["Experiment","Model 1.0 1pct","Simulation"];
pos2 = 0:9;
pos2 = 0.082 + 0.096*pos2;

%%%%%%% Fig. 2A: Qualitative heatmap, in-out relationships %%%%%%%
overall = figure('position',[100 40 1000 350]);
joined_inputs = unique(joined.Input);
for input = 1:length(joined_inputs)
    % subset_input = string(inputNames.input(input));
    subset = joined(string(joined.Input) == joined_inputs{input},[4 6 5 16]);    % filter for individual input and prediction cols
    subset.Properties.VariableNames(2:end) = xdata(:);
    subset = stack(subset,{'Experiment','Simulation'},'NewDataVariableNames','Value');
    n = length(unique(subset.Output));  % number of outputs for current subset
    % set subplot positions based on length of each subset
    offset = n/nmax;
    % plot heatmap in each subplot
    subplot('Position', [pos2(input) 0.21 0.0253 0.4950*offset]);
    fig = heatmap(subset,'Value_Indicator','Output','ColorVariable','Value');
    set(gca, ...
        "Colormap", cmap, ...
        "ColorbarVisible", "off", ...
        "ColorLimits", [-1 1], ...
        "CellLabelColor", 'none', ...
        "Title", joined_inputs{input}, ...
        "XLabel", "", ...
        "YLabel", "", ...
        "FontName","Arial")
    if input == 1
        set(gca,'YLabel',"Measured Output");
    end
end

%%%%%%% Fig. 2B: Qualitative heatmap, in-med relationships %%%%%%%
overall_med = figure('position',[100 400 1000 350]);
joined_inputs_med = unique(joined_med.Input);
for input = 1:length(joined_inputs_med)
    % subset_input = string(inputNames.input(input));
    subset = joined_med(string(joined_med.Input) == joined_inputs_med{input},[4 6 5 16]);    % filter for individual input and prediction cols
    subset.Properties.VariableNames(2:end) = xdata(:);
    subset = stack(subset,{'Experiment','Simulation'},'NewDataVariableNames','Value');
    n = length(unique(subset.Output));  % number of outputs for current subset
    % set subplot positions based on length of each subset
    offset = n/nmax;
    % plot heatmap in each subplot
    subplot('Position', [pos2(input) 0.21 0.0253 0.4950*offset]);
    fig = heatmap(subset,'Value_Indicator','Output','ColorVariable','Value');
    set(gca, ...
        "Colormap", cmap, ...
        "ColorbarVisible", "off", ...
        "ColorLimits", [-1 1], ...
        "CellLabelColor", 'none', ...
        "Title", joined_inputs_med{input}, ...
        "XLabel", "", ...
        "YLabel", "")
    if input == 1
        set(gca,'YLabel',"Measured Intermediate");
    end
end

saveas(overall,strcat(plotdir,"2A_ModelValidation_outputs.fig"));
saveas(overall_med,strcat(plotdir,"2B_ModelValidation_intermediates.fig"));
