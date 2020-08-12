% Drug Validation Analysis: Angiotensin Receptor Blocker
% 05.04.2020 JR
% 
% Purpose: simulate effect of ARB in both infarct and remote areas. Input
% levels are derived from idealized fits from rat infarcts (Zeigler et al.
% Matrix Biology 2020), and tension levels are varied manually within same
% range.

homedir = "data\";
valdir = "4_DrugCaseStudies\";
cmapdir = "utils\BrewerMap-master\";
plotdir = "plots\";
addpath(strcat(homedir,valdir));
addpath(cmapdir);

%% Calculate input levels
% infarct levels: using InputCurve script for t=2wks post-infarct
% remote levels: apply fraction to infarct levels, then add to baseline

peak = 0.4;                     % peak input level for model
weeks = 4;
t = weeks*7*24;                 % time in h
t0_analysis = 168;              % 168 h added as baseline
t_analysis = t+t0_analysis;     % matches InputCurve script
[InputCsim,~,inputNode_ramirez,~,~] = InputCurve_12_19(peak,peak);
inputs_baseline = InputCsim(:,t0_analysis);
inputs_ramirez_infarct = InputCsim(:,t_analysis);

inputs_fraction = 1;
inputs_ramirez_remote = (inputs_fraction*(inputs_ramirez_infarct-inputs_baseline))+inputs_baseline;

% manually construct input levels for von Leuder + Burke studies
inputs_vonleuder = 0.5;
inputNode_vonleuder = 1;

inputs_burke = [0.4 0.4 0.6];
inputNode_burke = [1 2 10];

%% Run simulations
% Ramirez study: using min/max for input levels
tension_remote = 0.1;
tension_infarct = 0.6;
perturblvls_ramirez = -0.5;
perturbnodes_ramirez = 2;
[sim_ramirez_remote,~] = DrugSimulation(inputs_ramirez_remote,inputNode_ramirez,perturblvls_ramirez,perturbnodes_ramirez,tension_remote);
[sim_ramirez_infarct,~] = DrugSimulation(inputs_ramirez_infarct,inputNode_ramirez,perturblvls_ramirez,perturbnodes_ramirez,tension_infarct);

% von Leuder study: perturbs estimated manually
tension_cc = 0.3;
perturblvls_vonleuder_npi = [0 0.05; 0 0.5; 0 1.5; 0 4];
perturblvls_vonleuder_val = [-0.01 0; -0.05 0; -0.15 0; -0.5 0];
perturblvls_vonleuder_both = [0 4; -0.01 0; -0.01 4; -0.05 0; -0.05 4; -0.15 0; -0.15 4; -0.5 0; -0.5 4];
perturbnodes_vonleuder = [2 27];
[sim_vonleuder_npi,~] = DrugSimulation(inputs_vonleuder,inputNode_vonleuder,perturblvls_vonleuder_npi,perturbnodes_vonleuder,tension_cc);
[sim_vonleuder_val,~] = DrugSimulation(inputs_vonleuder,inputNode_vonleuder,perturblvls_vonleuder_val,perturbnodes_vonleuder,tension_cc);
[sim_vonleuder_both,~] = DrugSimulation(inputs_vonleuder,inputNode_vonleuder,perturblvls_vonleuder_both,perturbnodes_vonleuder,tension_cc);

% Burke study: perturbs estimated manually
perturblvls_burke = [-0.5 0; 0 4; -0.5 4];
perturbnodes_burke = [2 27];
[sim_burke,~] = DrugSimulation(inputs_burke,inputNode_burke,perturblvls_burke,perturbnodes_burke,tension_cc);

%% Plot top-ranking changes in activity
load(strcat(homedir,valdir,"speciesNames.mat"), "speciesNames");

% top 10 downregulated nodes
% top = 10;
% [ranked_delta_infarct,ranked_idx_infarct] = sort(act_delta_infarct,'ascend');
% ranked_delta_remote = act_delta_remote(ranked_idx_infarct);
% figure("Position",[100 100 550 250]);
% bar(cat(1,ranked_delta_infarct(1:top),ranked_delta_remote(1:top)).');
% xticks(1:top);
% xticklabels(speciesNames(ranked_idx_infarct));
% xtickangle(45);
% ylabel("\DeltaActivity (KO - Baseline)");
% legend("Infarct Area","Remote Area","Location","best");
% title(strcat("Top ",string(top)," Downregulated Nodes"));
% saveas(gcf,strcat(homedir,valdir,"topnodes_downregulation.fig"));

% top 10 upregulated nodes
% [ranked_delta_infarct,ranked_idx_infarct] = sort(act_delta_infarct,'descend');
% ranked_delta_remote = act_delta_remote(ranked_idx_infarct);
% figure("Position",[200 100 550 250]);
% bar(cat(1,ranked_delta_infarct(1:top),ranked_delta_remote(1:top)).');
% xticks(1:top);
% xticklabels(speciesNames(ranked_idx_infarct));
% xtickangle(45);
% ylabel("\DeltaActivity (KO - Baseline)");
% legend("Infarct Area","Remote Area","Location","best");
% title(strcat("Top ",string(top)," Upregulated Nodes"));
% saveas(gcf,strcat(homedir,valdir,"topnodes_upregulation.fig"));
% close all

%% Import experimental data: von Leuder, Burke, Ramirez studies
filename = "snm_1_1_rev4.xlsx";
snmdir = strcat(homedir,filename);

filename = "experimentaldata_Ramirez_outputs.xlsx";
dbdir = strcat(homedir,valdir,filename);
[expt_ramirez,sheetnames_ramirez] = ImportSupplementalData(dbdir,snmdir);
expt_ramirez_all = [expt_ramirez{1};expt_ramirez{3}];
expt_ramirez_all(4,:) = [];     % remove duplicate condition (day0)
expt_ramirez_sem = [expt_ramirez{2};expt_ramirez{4}];
expt_ramirez_sem(4,:) = [];

filename = "experimentaldata_vonLeuder_outputs.xlsx";
dbdir = strcat(homedir,valdir,filename);
[expt_vonleuder,sheetnames_vonleuder] = ImportSupplementalData(dbdir,snmdir);
% expt_vonleuder_all = [expt_vonleuder{1};expt_vonleuder{3}];
% expt_vonleuder_sem = [expt_vonleuder{2};expt_vonleuder{4}];

filename = "experimentaldata_Burke_outputs.xlsx";
dbdir = strcat(homedir,valdir,filename);
[expt_burke,sheetnames_burke] = ImportSupplementalData(dbdir,snmdir);


%% Process simulation + experimental data
% find output indices from experimental data variables
outputs_names = expt_ramirez_all.Properties.VariableNames;
outputs_idx_ramirez = [];
for name = outputs_names(2:end)
    outputs_idx_ramirez = [outputs_idx_ramirez find(strcmp(speciesNames,name))];
end
outputs_idx_ramirez2 = outputs_idx_ramirez([1:3 5:6]);

outputs_names = expt_vonleuder{1}.Properties.VariableNames;
outputs_idx_vonleuder = [];
for name = outputs_names(2:end)
    outputs_idx_vonleuder = [outputs_idx_vonleuder find(strcmp(speciesNames,name))];
end

outputs_names = expt_burke{1}.Properties.VariableNames;
outputs_idx_burke = [];
for name = outputs_names(2:end)
    outputs_idx_burke = [outputs_idx_burke find(strcmp(speciesNames,name))];
end

% subset simulation data for outputs + reshape + normalize
normalize = @(x) (x ./ x(1,:));

sim_ramirez_all = [sim_ramirez_remote(:,outputs_idx_ramirez2);sim_ramirez_infarct(2:end,outputs_idx_ramirez2)];
sim_ramirez_all = normalize(sim_ramirez_all);

sim_vonleuder_npi = sim_vonleuder_npi(:,outputs_idx_vonleuder);
sim_vonleuder_npi = normalize(sim_vonleuder_npi);
sim_vonleuder_val = sim_vonleuder_val(:,outputs_idx_vonleuder);
sim_vonleuder_val = normalize(sim_vonleuder_val);
sim_vonleuder_both_outputs = sim_vonleuder_both(:,outputs_idx_vonleuder);
sim_vonleuder_both_outputs = normalize(sim_vonleuder_both_outputs);
sim_vonleuder_both_reordered = [
    [sim_vonleuder_both_outputs(1) 0]; sim_vonleuder_both_outputs([2 3]).'; 
    sim_vonleuder_both_outputs([4 5]).'; sim_vonleuder_both_outputs([6 7]).'; 
    sim_vonleuder_both_outputs([8 9]).'; sim_vonleuder_both_outputs([10 11]).'];
expt_vonleuder_both_reordered = [
    [expt_vonleuder{1,5}{1,2} 0]; expt_vonleuder{1,5}{[2 3],2}.';
    expt_vonleuder{1,5}{[4 5],2}.'; expt_vonleuder{1,5}{[6 7],2}.';
    expt_vonleuder{1,5}{[8 9],2}.'; expt_vonleuder{1,5}{[8 9],2}.'];
expt_vonleuder_both_error = [
    [expt_vonleuder{1,6}{1,2} 0]; expt_vonleuder{1,6}{[2 3],2}.';
    expt_vonleuder{1,6}{[4 5],2}.'; expt_vonleuder{1,6}{[6 7],2}.';
    expt_vonleuder{1,6}{[8 9],2}.'; expt_vonleuder{1,6}{[8 9],2}.'];

sim_burke_outputs = sim_burke([2:3 5],outputs_idx_burke);
sim_burke_outputs = normalize(sim_burke_outputs);
sim_burke_pkg = sim_burke(2:5,30);
sim_burke_pkg = normalize(sim_burke_pkg);
sim_burke_rho = sim_burke(:,33);
sim_burke_rho = normalize(sim_burke_rho);

%% Visualize simulation vs. experimental data: Ramirez
expt_ramirez_plot = expt_ramirez_all{:,[2:4 6:7]};
expt_ramirez_plotsem = expt_ramirez_sem{:,[2:4 6:7]};
expt_ramirez_plotsem = expt_ramirez_plotsem ./ expt_ramirez_plot(1,:);
expt_ramirez_plot = normalize(expt_ramirez_plot);
xlabs_sim = ["Control" "LVC" "LVC+AT1Ri" "LVI" "LVI+AT1Ri"];
xlabs = ["Control" "LVC" "LVC+VAL" "LVI" "LVI+VAL"];
ymax = (max([sim_ramirez_all;expt_ramirez_plot])+max(expt_ramirez_plotsem))*1.1;

figure("Position",[50 25 425 750]);
for i = 1:length(sim_ramirez_all)
    subplot(length(sim_ramirez_all),2,(2*i)-1);
    bar(sim_ramirez_all(:,i)); hold on
    % ylim([0 max(sim_ramirez_all(:,i),[],'all')*1.1]);
    ylim([0 ymax(i)]);
    ylabel(speciesNames{outputs_idx_ramirez2(i)}); 
    set(gca,"YGrid",'on',"FontSize",10);
    if i == length(sim_ramirez_all)
        xticklabels(xlabs_sim); xtickangle(45);
    else
        xticklabels({});
    end
    if i == 1
        title("Simulation Data"); 
    end
    hold off
end

expt_ramirez_genes = ["Col1a1" "Col3a1" "Ctgf" "Postn" "Spp1"];
for i = 1:length(expt_ramirez_plot)
    subplot(length(expt_ramirez_plot),2,2*i);
    bar(expt_ramirez_plot(:,i),'FaceColor',[.3 .3 .3]); hold on
    errorbar(expt_ramirez_plot(:,i),expt_ramirez_plotsem(:,i),"LineStyle","none","Color","k");
    % ymax = (max(expt_ramirez_plot(:,i),[],'all')+max(expt_ramirez_plotsem(:,i),[],'all'))*1.05;
    ylim([0 ymax(i)]);
    ylabel(expt_ramirez_genes(i)); 
    set(gca,"YGrid",'on',"FontSize",10);
    if i == length(expt_ramirez_plot)
        xticklabels(xlabs); xtickangle(45); 
    else
        xticklabels({});
    end
    if i == 1
        title("Ramirez 2014"); 
    end
    hold off
end
saveas(gcf,strcat(plotdir,"5C_Ramirez_outputs.fig"));


%% Visualize experimental vs. simulation data: Von Leuder
%%%%%% Figure S5A: NEPi dose-responses comparison %%%%%%
expt_vonleuder_npi = expt_vonleuder{1}{:,2};
expt_vonleuder_npisem = expt_vonleuder{2}{:,2};
expt_vonleuder_npisem = expt_vonleuder_npisem ./ expt_vonleuder_npi(1,:);
expt_vonleuder_npi = normalize(expt_vonleuder_npi);
xlabs_sim = ["Control";"0 NEPi";strcat(string(perturblvls_vonleuder_npi(:,2))," NEPi")];
xlabs_expt = expt_vonleuder{1}{:,1};
xlabs_split = split(xlabs_expt(3:end),'+');
xlabs_expt(3:end) = xlabs_split(:,2); xlabs_expt(2) = {'0LBQ'};
ymax = max([sim_vonleuder_npi;expt_vonleuder_npi])*1.1;

figure("Position",[100 100 500 250]);
subplot(1,2,1);
bar(sim_vonleuder_npi);
ylabel('proCI Activation'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_sim); xtickangle(45); ylim([0 ymax]);
title("Model Prediction");
subplot(1,2,2);
bar(expt_vonleuder_npi,'FaceColor',[.3 .3 .3]); hold on
errorbar(expt_vonleuder_npi,expt_vonleuder_npisem,"LineStyle","none","Color","k");
ylabel('H^{3}-Proline Incorporation'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_expt); xtickangle(45); ylim([0 ymax]);
title("von Leuder 2015"); hold off
saveas(gcf,strcat(plotdir,"S5A_vonLeuder_NEPi.fig"));

%%%%%% Figure S5B: ARB dose-response comparison %%%%%%
expt_vonleuder_val = expt_vonleuder{3}{:,2};
expt_vonleuder_valsem = expt_vonleuder{4}{:,2};
expt_vonleuder_valsem = expt_vonleuder_valsem ./ expt_vonleuder_val(1,:);
expt_vonleuder_val = normalize(expt_vonleuder_val);
xlabs_sim = ["Control";"0 AT1Ri";strcat(string(-1*perturblvls_vonleuder_val(:,1))," AT1Ri")];
xlabs_expt = expt_vonleuder{3}{:,1};
xlabs_split = split(xlabs_expt(3:end),'+');
xlabs_expt(3:end) = xlabs_split(:,2); xlabs_expt(2) = {'0VAL'};
ymax = max([sim_vonleuder_val;expt_vonleuder_val])*1.1;

figure("Position",[100 100 500 250]);
subplot(1,2,1);
bar(sim_vonleuder_val);
ylabel('proCI Activity'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_sim); xtickangle(45); ylim([0 ymax]);
title("Model Prediction");
subplot(1,2,2);
bar(expt_vonleuder_val,'FaceColor',[.3 .3 .3]); hold on
errorbar(expt_vonleuder_val,expt_vonleuder_valsem,"LineStyle","none","Color","k");
ylabel('H^{3}-Proline Incorporation'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_expt); xtickangle(45); ylim([0 ymax]);
title("von Leuder 2015"); hold off
saveas(gcf,strcat(plotdir,"S5B_vonLeuder_ARB.fig"));

%%%%%% Figure 5A: ARB-NEPi combination comparison %%%%%%
expt_vonleuder_both_error = expt_vonleuder_both_error ./ expt_vonleuder_both_reordered(1);
expt_vonleuder_both_reordered = expt_vonleuder_both_reordered ./ expt_vonleuder_both_reordered(1);
xlabs_expt = ["Control","0 AT1Ri","0.01 AT1Ri","0.05 AT1Ri","0.15 AT1Ri","0.5 AT1Ri"];
legend_expt = ["-NEPi","+NEPi"];
ymax = max([sim_vonleuder_both_reordered;expt_vonleuder_both_reordered],[],'all')*1.1;
colors = flipud(brewermap(2,'Paired'));

figure("Position",[100 100 700 330]);
subplot(3,2,[1 3]);
bar(sim_vonleuder_both_reordered);
ylabel('proCI Activity'); 
set(gca,"YGrid",'on',"FontSize",10,"ColorOrder",colors);
xticklabels(xlabs_expt); xtickangle(45); ylim([0 ymax]);
title("Model Prediction");
subplot(3,2,5);
bar(sim_vonleuder_both_reordered,"Visible",false);
set(gca,"ColorOrder",colors,"FontSize",10);
axis off
legend(legend_expt,"Location","southoutside","NumColumns",2);

xlabs_expt = ["Control","0 VAL","0.03 VAL","0.1 VAL","0.3 VAL","1 VAL"];
legend_expt = ["-LBQ","+LBQ"];
colors = [.3 .3 .3; .6 .6 .6];
subplot(3,2,[2 4]);
hbar = bar(expt_vonleuder_both_reordered); hold on
for k1 = 1:size(expt_vonleuder_both_reordered,2)
    ctr(k1,:) = bsxfun(@plus, hbar(k1).XData, hbar(k1).XOffset');   % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = hbar(k1).YData;                                     % Individual Bar Heights
end
errorbar(ctr,ydt,expt_vonleuder_both_error.',"LineStyle","none","Color","k");
ylabel('H^{3}-Proline Incorporation'); 
set(gca,"YGrid",'on',"FontSize",10,"ColorOrder",colors);
xticklabels(xlabs_expt); xtickangle(45); ylim([0 ymax]);
title("von Leuder 2015"); hold off
subplot(3,2,6);
bar(expt_vonleuder_both_reordered,"Visible",false);
set(gca,"ColorOrder",colors,"FontSize",10);
legend(legend_expt,"Location","southoutside","NumColumns",2);
axis off
saveas(gcf,strcat(plotdir,"5A_vonLeuder_ARB-NEPi.fig"));


%% Visualize experimental vs. simulation data: Burke
%%%%%% Figure 5B: output comparison %%%%%%
expt_burke_plot = expt_burke{1}{:,2:end}.';
expt_burke_sem = expt_burke{2}{:,2:end}.';
legend_sim = ["T/A/A";strcat("T/A/A + ", ...
    string(abs(perturblvls_burke([1 3],1))),"AT1Ri + ", ...
    string(abs(perturblvls_burke([1 3],2))),"NEPi")];
xlabs_sim = expt_burke{1}.Properties.VariableNames(2:end);
legend_expt = expt_burke{1}{:,1};
ymax = (max([sim_burke_outputs;expt_burke_plot],[],'all')+max(expt_burke_sem,[],'all'))*1.1;
colors = redbluecmap(7);

figure("Position",[100 100 700 300]);
subplot(3,2,[1 3]);
bar(sim_burke_outputs.');
ylabel('Node Activity'); 
set(gca,"YGrid",'on',"FontSize",10,"ColorOrder",colors(1:3,:));
% legend(legend_sim,"Location","southoutside","NumColumns",1);
xticklabels(xlabs_sim); xtickangle(45); ylim([0 ymax]);
title("Model Prediction");

subplot(3,2,[2 4]);
bar(expt_burke_plot); hold on
% Set errorbar x positions
ngroups = size(expt_burke_plot, 1);
nbars = size(expt_burke_plot, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, expt_burke_plot(:,i), expt_burke_sem(:,i), 'k', 'LineStyle', 'none');
end
ylabel('mRNA Expression'); 
set(gca,"YGrid",'on',"FontSize",10,"ColorOrder",[.2 .2 .2; .5 .5 .5; .8 .8 .8]);
xticklabels(xlabs_sim); xtickangle(45); ylim([0 ymax]);
title("Burke 2019"); hold off

subplot(3,2,5);
bar(sim_burke_outputs.',"Visible",false);
set(gca,"ColorOrder",colors(1:3,:),"FontSize",10);
legend(legend_sim,"Location","southoutside","NumColumns",1);
axis off
subplot(3,2,6);
bar(sim_burke_outputs.',"Visible",false);
set(gca,"ColorOrder",[.2 .2 .2; .5 .5 .5; .8 .8 .8],"FontSize",10);
legend(legend_expt,"Location","southoutside","NumColumns",1);
axis off
saveas(gcf,strcat(plotdir,"5B_Burke_outputs.fig"));


%%%%%% Figure S5C: PKG activation comparison %%%%%%
expt_burke_plot = expt_burke{3}{:,2};
expt_burke_sem = expt_burke{4}{:,2};
xlabs_sim = ["T/A/A";strcat("T/A/A+", ...
    string(abs(perturblvls_burke(:,1))),"AT1Ri+", ...
    string(abs(perturblvls_burke(:,2))),"NEPi")];
xlabs_expt = expt_burke{3}{:,1};

figure("Position",[100 100 500 250]);
subplot(1,2,1);
bar(sim_burke_pkg);
ymax = max(sim_burke_pkg)*1.1;
ylabel('PKG Activity'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_sim); xtickangle(45); ylim([0 ymax]);
title("Model Prediction");

subplot(1,2,2);
bar(expt_burke_plot,'FaceColor',[.3 .3 .3]); hold on
errorbar(expt_burke_plot,expt_burke_sem,"LineStyle","none","Color","k");
ymax = (max(expt_burke_plot)+max(expt_burke_sem))*1.1;
ylabel('PKG Activity (FC)'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_expt); xtickangle(45); ylim([0 ymax]);
title("Burke 2019"); hold off
saveas(gcf,strcat(plotdir,"S5C_Burke_PKG.fig"));

%%%%%% Figure S5D: Rho activity comparison %%%%%%
expt_burke_plot = expt_burke{5}{:,2};
expt_burke_sem = expt_burke{6}{:,2};
xlabs_sim = ["T/A/A";strcat("T/A/A+", ...
    string(abs(perturblvls_burke(:,1))),"AT1Ri+", ...
    string(abs(perturblvls_burke(:,2))),"NEPi")];
xlabs_expt = expt_burke{5}{:,1};
ymax = max(sim_burke_rho)*1.1;
colors = get(gca,"ColorOrder");

figure("Position",[100 100 500 250]);
subplot(1,2,1);
bar(sim_burke_rho);
ymax = max(sim_burke_rho)*1.1;
ylabel('Rho Activity'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_sim); xtickangle(45); ylim([0 ymax]);
title("Model Prediction");

subplot(1,2,2);
bar(expt_burke_plot,'FaceColor',[.3 .3 .3]); hold on
errorbar(expt_burke_plot,expt_burke_sem,"LineStyle","none","Color","k");
ymax = (max(expt_burke_plot)+max(expt_burke_sem))*1.1;
ylabel('Rho (active)/Rho (total)'); set(gca,"YGrid",'on',"FontSize",10);
xticklabels(xlabs_expt); xtickangle(45); ylim([0 ymax]);
title("Burke 2019"); hold off
saveas(gcf,strcat(plotdir,"S5D_Burke_Rho.fig"));

