% Model validation: SNM 1.1rev4
% Created 14 Apr 2020 JR

homedir = "data\";
valdir = "1_ModelValidation\";
plotdir = "plots\";
addpath(homedir)
addpath(strcat(homedir,valdir));
% Import simulation data
fprintf("Loading Data...\n")
filedir = "snm_1_1_rev4_val_paramsweep_inputs.mat";
load(filedir)
filedir = "snm_1_1_rev4_val_paramsweep_predictions.mat";
load(filedir)
filedir = "snm_1_1_rev4_val_paramsweep_speciesNames.mat";
load(filedir)

%% Import/process experimental validation (expt)
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

modelpath = 'snm_1_1_rev4.xlsx';
sheetname = "reactions";
opts = detectImportOptions(modelpath,"Sheet",sheetname);
snm = readtable(modelpath, opts);
idx_rxns = find(strcmpi(snm.module,"output"));
snm_rxns = split(snm.Rule(idx_rxns)," => ");
snm_outputs = unique(snm_rxns(:,end));
outputs_idx = zeros(length(snm_outputs),1);
% ID model outputs
for output = 1:length(snm_outputs)
    outputs_idx(output) = find(strcmp(speciesNames,snm_outputs(output)));
end

%% Process simulation predictions (pred)
% 05.11.2020 JR: loop for parameter sweep values
n_sweep = 1.05:0.05:1.75;
EC50_sweep = 0.2:0.05:0.8;
correct_sweep = [];
inputNames = cell2table(inputNames.','VariableNames',{'input'});

for n_i = 1:size(predictions,3)
    for EC50_i = 1:size(predictions,4)
        % create prediction tables
        pred_full = array2table(predictions(:,:,n_i,EC50_i),'VariableNames',speciesNames.');
        pred_full = [inputNames pred_full];
        % subset for outputs only
        outputs_pred = pred_full.Properties.VariableNames(outputs_idx+1);
        pred = pred_full(:,outputs_idx+1);
        % subset for intermediates
        pred_med = pred_full;
        pred_med(:,outputs_idx+1) = [];
        
        %% Combine expt/pred tables (in-out):
        % Convert pred table to long-form
        pred.Input = pred_full.input;
        pred_long = stack(pred,1:width(pred)-1);    % stack all vars except input
        pred_long.Properties.VariableNames = {'Input','Output','SNM_1_1'};
        pred_long.ID_pred = (1:height(pred_long)).';
        pred_long.Output = cellstr(string(pred_long.Output));    % match class with expt
        % join pred with expt based on matching fields: Input, Output
        joined = innerjoin(expt,pred_long,'Keys',{'Input','Output'});
        % Test for equivalence between measurements/predictions
        joined.Correct_1_1 = (joined.Measurement == joined.SNM_1_1);
        % Get percentage of correct predictions
        s = summary(joined);
        pct_1_1 = s.Correct_1_1.True / s.Correct_1_1.Size(1);
        
        %% Combine expt/pred tables (in-med):
        % Convert pred table to long-form
        pred_med_long = stack(pred_med,2:width(pred_med));    % stack all vars except input
        pred_med_long.Properties.VariableNames = {'Input','Output','SNM_1_1'};
        pred_med_long.ID_pred = [1:height(pred_med_long)].';
        pred_med_long.Output = cellstr(string(pred_med_long.Output));    % match class with expt
        % join pred with expt based on matching fields: Input, Output
        joined_med = innerjoin(expt_med,pred_med_long,'Keys',{'Input','Output'});
        
        % Test for equivalence between measurements/predictions
        joined_med.Correct_1_1 = (joined_med.Measurement == joined_med.SNM_1_1);
        
        % Get percentage of correct predictions
        s_med = summary(joined_med);
        pct_1_1_med = s_med.Correct_1_1.True / s_med.Correct_1_1.Size(1);
        
        % Save results
        correct_sweep = [correct_sweep;n_sweep(n_i) EC50_sweep(EC50_i) pct_1_1 pct_1_1_med];
        fmt = 'n=%.2f | EC50=%.2f | In-Out: %2.3f | In-Med: %2.3f\n';
        % Keep highest-scoring conditions (based on total accuracy)
        accuracy_cur_total = pct_1_1+pct_1_1_med;
        if n_i==1 && EC50_i==1
            accuracy_prev_total = 0;
        else
            accuracy_prev = max(correct_sweep(1:end-1,:),[],1);
            accuracy_prev_total = accuracy_prev(end-1)+accuracy_prev(end);
        end
        if accuracy_cur_total > accuracy_prev_total
            joined_keep = joined;
            joined_keep_med = joined_med;
            params_keep = [n_sweep(n_i) EC50_sweep(EC50_i)];
            fprintf(fmt,n_sweep(n_i),EC50_sweep(EC50_i),pct_1_1,pct_1_1_med)
        end
    end
end

%% Visualize predictions:
%%%%%%% Fig. S2: parameter sweep results %%%%%%%
correct_sweep = array2table(correct_sweep,'VariableNames',["n" "EC50" "in-out" "in-med"]);
fig_sweep = figure('position',[100 200 750 250]);
subplot(1,2,1);
fig = heatmap(correct_sweep,"n","EC50","ColorVariable","in-out");
title("Input-Output Accuracy");
subplot(1,2,2);
fig = heatmap(correct_sweep,"n","EC50","ColorVariable","in-med");
title("Input-Intermediate Accuracy");
saveas(fig_sweep,strcat(plotdir,"S2_ModelValidation_paramsweep.fig"));

%%%%%%% Fig. 2A: Qualitative heatmap, in-out relationships %%%%%%%
joined = sortrows(joined_keep,{'Input','Output'});   %arrange rows for uniformity
joined_inputs = unique(joined_keep.Input);
joined_med = sortrows(joined_keep_med,{'Input','Output'});
joined_inputs_med = unique(joined_keep_med.Input);
cmap = [0 0.502 0.752; 0.75 0.75 0.75; 1 0 0];
nmax = 14;  % maximum number of outputs  for all subsets
xdata = ["Experiment","Model 1.0 1pct","Current Model"];
pos2 = 0:9;
pos2 = 0.94 - 0.096*pos2;

overall = figure('position',[100 40 425 750]);
for input = 1:length(joined_inputs)
    subset_input = joined_inputs{input};
    subset = joined_keep(string(joined_keep.Input) == subset_input,[4 6 5 16]);    % filter for individual input and prediction cols
    subset.Properties.VariableNames(2:end) = xdata(:);
    subset = stack(subset,{'Experiment','Current Model'},'NewDataVariableNames','Value');
    n = length(unique(subset.Output));  % number of outputs for current subset
    % set subplot positions based on length of each subset
    offset = n/nmax;
    % plot heatmap in each subplot
    col = input;
    subplot(length(joined_inputs),1,col);
    fig = heatmap(subset,'Output','Value_Indicator','ColorVariable','Value');
    pos = get(fig,"Position");
    pos(1) = pos(1)+0.05;
    pos(2) = pos2(col);
    pos(3) = (pos(3)-0.08)*offset;
    pos(4) = pos(4)-0.035;
    set(gca,"Position", pos, ...
        "Colormap", cmap, ...
        "ColorbarVisible", "off", ...
        "ColorLimits", [-1 1], ...
        "CellLabelColor", 'none', ...
        "Title", subset_input, ...
        "XLabel", "", ...
        "YLabel", "", ...
        "FontSize", 7)
    if input == length(joined_inputs)
        set(gca,'XLabel',"Measured Output");
    end
end

%%%%%%% Fig. 2B: Qualitative heatmap, in-med relationships %%%%%%%
overall_med = figure('position',[200 40 400 750]);
for input = 1:length(joined_inputs_med)
    subset_input = joined_inputs_med{input};
    subset = joined_keep_med(string(joined_keep_med.Input) == subset_input,[4 6 5 16]);    % filter for individual input and prediction cols
    subset.Properties.VariableNames(2:end) = xdata(:);
    subset = stack(subset,{'Experiment','Current Model'},'NewDataVariableNames','Value');
    n = length(unique(subset.Output));  % number of outputs for current subset
    % set subplot positions based on length of each subset
    offset = n/nmax;
    
    % plot heatmap in each subplot
    col = input;
    subplot(length(joined_inputs),1,col);
    fig = heatmap(subset,'Output','Value_Indicator','ColorVariable','Value');
    pos = get(fig,"Position");
    pos(1) = pos(1)+0.1;
    pos(2) = pos2(col);
    pos(3) = (pos(3)-0.08)*offset;
    pos(4) = pos(4)-0.035;
    set(gca,"Position", pos, ...
        "Colormap", cmap, ...
        "ColorbarVisible", "off", ...
        "ColorLimits", [-1 1], ...
        "CellLabelColor", 'none', ...
        "Title", subset_input, ...
        "XLabel", "", ...
        "YLabel", "", ...
        "FontSize", 7)
    if input == length(joined_inputs_med)
        set(gca,'XLabel',"Measured Output");
    end
end

