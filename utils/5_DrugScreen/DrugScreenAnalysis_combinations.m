% Drug Screen: results visualization
% 01.12.2019 JR
% 
% Purpose: combine low-/high-tension screens to identify perturbations that
% decrease matrix content in low tension and increase matrix content in
% high tension

%% Import + Process Data
homedir = "D:/Research/Aim2/ModelExpansion/1_1/rev4/";
sensdir = "DrugScreen/";
datadir = "drugscreen_data/2wks_final/";
cmapdir = "D:/Research/Aim2/BrewerMap-master/";
addpath(homedir);
addpath(strcat(homedir,sensdir));
addpath(strcat(homedir,sensdir,datadir));
addpath(cmapdir);

modelpath = 'snm_1_1_rev4.xlsx';
opts = detectImportOptions(modelpath);
snm = readtable(modelpath, opts);

% method 1: use all outputs, counting MMPs as anti-fibrotic and rest as
% pro-fibrotic
idx_ecm = find(strcmpi(snm.module,"ECM"));
idx_ecm_pro = idx_ecm(~contains(snm.ID(idx_ecm),"MMP"));
idx_ecm_anti = idx_ecm(contains(snm.ID(idx_ecm),"MMP"));
idx_combinations = nchoosek(1:length(snm.ID),2);
% method 2: limit outputs to those with direct remodeling capability
idx_ecm_pro2 = idx_ecm([3:5 10:13 19]);
idx_ecm_anti2 = idx_ecm([6:9 15:17]);
idx_ecm2 = [idx_ecm_pro2;idx_ecm_anti2];
% method 3: separate proteases/inhibs and add mutual inhibition
idx_ecm_pro3 = idx_ecm([10:13 19]);
idx_ecm_inhib = idx_ecm(3:5);

% remove combinations modulating outputs
remove_values = ismember(idx_combinations,idx_ecm);
remove_rows = (remove_values(:,1) == 1 | remove_values(:,2) == 1);
nodes_all = repmat(idx_combinations(~remove_rows,:),4,1);
len_c = length(idx_combinations(~remove_rows));

% import/filter data
makeFileName = @(x,y) strcat("combinations_tension",replace(string(x),".",""),"_",y,".mat");

tension = 0.1; 
choice = "KO";
act_lt_ko = load(makeFileName(tension,choice));
choice = "OE";
act_lt_oe = load(makeFileName(tension,choice));
choice = "KO-OE";
act_lt_kooe = load(makeFileName(tension,choice));
choice = "OE-KO";
act_lt_oeko = load(makeFileName(tension,choice));
act_lt_ecm_ko = act_lt_ko.act_delta(~remove_rows,idx_ecm2);
act_lt_ecm_oe = act_lt_oe.act_delta(~remove_rows,idx_ecm2);
act_lt_ecm_kooe = act_lt_kooe.act_delta(~remove_rows,idx_ecm2);
act_lt_ecm_oeko = act_lt_oeko.act_delta(~remove_rows,idx_ecm2);
act_lt = [act_lt_ko.act_delta(~remove_rows,:);act_lt_oe.act_delta(~remove_rows,:);
    act_lt_kooe.act_delta(~remove_rows,:);act_lt_oeko.act_delta(~remove_rows,:)];

% high-tension data
tension = 0.6;
choice = "KO";
act_ht_ko = load(makeFileName(tension,choice));
choice = "OE";
act_ht_oe = load(makeFileName(tension,choice));
choice = "KO-OE";
act_ht_kooe = load(makeFileName(tension,choice));
choice = "OE-KO";
act_ht_oeko = load(makeFileName(tension,choice));
act_ht_ecm_ko = act_ht_ko.act_delta(~remove_rows,idx_ecm2);
act_ht_ecm_oe = act_ht_oe.act_delta(~remove_rows,idx_ecm2);
act_ht_ecm_kooe = act_ht_kooe.act_delta(~remove_rows,idx_ecm2);
act_ht_ecm_oeko = act_ht_oeko.act_delta(~remove_rows,idx_ecm2);
act_ht = [act_ht_ko.act_delta(~remove_rows,:);act_ht_oe.act_delta(~remove_rows,:);
    act_ht_kooe.act_delta(~remove_rows,:);act_ht_oeko.act_delta(~remove_rows,:)];

%% Calculate overall matrix content change
% Method 1: calc. matrix content change metric from
% pro/anti/inhib outputs
lt_mcc = calculateMCC(act_lt,idx_ecm_pro2,idx_ecm_anti2,[]);
ht_mcc = calculateMCC(act_ht,idx_ecm_pro2,idx_ecm_anti2,[]);
% Method 2: rank matrix content changes in pro/anti/inhib categories
[lt_pro,lt_anti,lt_inhib] = calculatePAI(act_lt,idx_ecm_pro3,idx_ecm_anti2,idx_ecm_inhib);
[ht_pro,ht_anti,ht_inhib] = calculatePAI(act_ht,idx_ecm_pro3,idx_ecm_anti2,idx_ecm_inhib);
lt_scores = scorePAI(lt_pro,lt_anti,lt_inhib,"lt");
ht_scores = scorePAI(ht_pro,ht_anti,ht_inhib,"ht");

% Find mechano-adaptive perturbatons
% Method 1: ID perturbs from calculated MCC using threshold
threshold = 0;
idx_threshold = find(lt_mcc <= -threshold & ht_mcc >= threshold);

% Method 2: rank perturbs based on MCC, find top-ranked perturbs
[~,sort_lt_idx] = sort(lt_mcc(idx_threshold),"ascend");   % sort arrays
[~,sort_ht_idx] = sort(ht_mcc(idx_threshold),"descend");
sort_lt_idx(:,2) = 1:length(sort_lt_idx);       % add ranks
sort_ht_idx(:,2) = 1:length(sort_ht_idx);
rank_lt_idx = sortrows(sort_lt_idx,1);          % reorder to match lt/ht rows
rank_ht_idx = sortrows(sort_ht_idx,1);
score = @(x,y) length([x;y])-(x+y);
score_idx = [rank_lt_idx(:,1) score(rank_lt_idx(:,2),rank_ht_idx(:,2))];  % sum lt/ht ranks
sort_score_idx = sortrows(score_idx,2,"descend");         % sort by score
idx_adaptive_grp1 = idx_threshold(sort_score_idx(1:13,1)).';
idx_adaptive_grp2 = idx_threshold(sort_score_idx(14:156,1)).';

% Method 3: rank perturbs based on pro/anti/inhib scores independently,
% combine rankings for lt/ht and find top-ranked perturbs
score_both = [lt_scores(:,1) mean([lt_scores(:,2) ht_scores(:,2)],2)];
sort_score_both = sortrows(score_both,2,"descend");
idx_adaptive_3 = sort_score_both(1:50,1).';

% Method 3b: combine pro/anti/inhib rankings for lt/ht together, and find
% top-ranked perturbs
score_both2 = scoreAll(lt_pro,lt_anti,lt_inhib,ht_pro,ht_anti,ht_inhib);
sort_score_both2 = sortrows(score_both2,2,"descend");
idx_adaptive_3b = sort_score_both(1:50,1).';

% Find perturb targets/types for each hit
nodes_adaptive = nodes_all(idx_adaptive_grp1,:);     % filter node combinations
nodes_type = string(zeros(length(idx_adaptive_grp1),1));

for node = 1:length(idx_adaptive_grp1)
    if idx_adaptive_grp1(node)>len_c && idx_adaptive_grp1(node)<=2*len_c
        nodes_type(node) = "OE";
    elseif idx_adaptive_grp1(node)>2*len_c && idx_adaptive_grp1(node)<=3*len_c
        nodes_type(node) = "KO-OE";
    elseif idx_adaptive_grp1(node)>3*len_c && idx_adaptive_grp1(node)<=4*len_c
        nodes_type(node) = "OE-KO";
    else
        nodes_type(node) = "KO";
    end
end
ID_adaptive = snm.ID(nodes_adaptive(:,1));      % translate to node IDs
ID_adaptive(:,2) = snm.ID(nodes_adaptive(:,2));
changes_adaptive = lt_mcc(idx_adaptive_grp1);   % filter matrix changes
changes_adaptive(:,2) = ht_mcc(idx_adaptive_grp1);

% display data as a table
table_adaptive = cell2table(ID_adaptive);
table_adaptive(:,3:4) = array2table(changes_adaptive);
table_adaptive.Properties.VariableNames = ["Node1","Node2","LowTension","HighTension"];
table_adaptive.perturbation = nodes_type;

%% Plot Results
% plot overall matrix content change
% xdata = 1:length(all_lt_diff);
% [ydata_lt,idx_ydata] = sort(all_lt_diff,'ascend');
% ydata_ht = all_ht_diff(idx_ydata);
xdata = score_idx(:,2);
figure('position',[100 100 600 300]);
fig_lt = scatter(xdata,lt_mcc(idx_threshold),4); hold on
fig_ht = scatter(xdata,ht_mcc(idx_threshold),4);
set(gca,"FontSize",12,"XLim",[0 1000]);
xticks(0:100:1000);
legend('Low Tension','High Tension','Location','southwest');
xlabel('Mechano-Adaptive Score'); ylabel('Matrix Content Change');
annotation('rectangle',[0.76 0.35 0.07 0.5], ...
    'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--')
grid on

% plot individual output changes for mechano-adaptive candidates
map = brewermap(19,'*RdBu');
act_lt_adaptive_plot = [];
act_ht_adaptive_plot = [];
labels = [];
for i = 1:length(idx_adaptive_grp1)
    if idx_adaptive_grp1(i) > len_c && idx_adaptive_grp1(i) <= 2*len_c        % for OE idxs
        idx_adj = idx_adaptive_grp1(i) - len_c;
        act_lt_adaptive_plot = [act_lt_adaptive_plot;act_lt_ecm_oe(idx_adj,:)];
        act_ht_adaptive_plot = [act_ht_adaptive_plot;act_ht_ecm_oe(idx_adj,:)];
        label_adj = snm.ID(nodes_adaptive(i,1))+" OE + "+snm.ID(nodes_adaptive(i,2))+" OE";
        labels = [labels label_adj];
    elseif idx_adaptive_grp1(i) > 2*len_c && idx_adaptive_grp1(i) <= 3*len_c  % for KO-OE
        idx_adj = idx_adaptive_grp1(i) - 2*len_c;
        act_lt_adaptive_plot = [act_lt_adaptive_plot;act_lt_ecm_kooe(idx_adj,:)];
        act_ht_adaptive_plot = [act_ht_adaptive_plot;act_ht_ecm_kooe(idx_adj,:)];
        label_adj = snm.ID(nodes_adaptive(i,1))+" KD + "+snm.ID(nodes_adaptive(i,2))+" OE";
        labels = [labels label_adj];
    elseif idx_adaptive_grp1(i) > 3*len_c && idx_adaptive_grp1(i) <= 4*len_c  % for OE-KO
        idx_adj = idx_adaptive_grp1(i) - 3*len_c;
        act_lt_adaptive_plot = [act_lt_adaptive_plot;act_lt_ecm_oeko(idx_adj,:)];
        act_ht_adaptive_plot = [act_ht_adaptive_plot;act_ht_ecm_oeko(idx_adj,:)];
        label_adj = snm.ID(nodes_adaptive(i,1))+" OE + "+snm.ID(nodes_adaptive(i,2))+" KD";
        labels = [labels label_adj];
    else
        idx_adj = idx_adaptive_grp1(i);                                  % for KO idxs
        act_lt_adaptive_plot = [act_lt_adaptive_plot;act_lt_ecm_ko(idx_adj,:)];
        act_ht_adaptive_plot = [act_ht_adaptive_plot;act_ht_ecm_ko(idx_adj,:)];
        label_adj = snm.ID(nodes_adaptive(i,1))+" KD + "+snm.ID(nodes_adaptive(i,2))+" KD";
        labels = [labels label_adj];
    end
end

fig = figure("Position",[80 80 400 675]);
subplot("Position",[0.25 0.66 0.6 0.29]);
heatmap(act_lt_adaptive_plot.');
% set(gca,"ColorbarVisible","off");
set(gca,"XDisplayLabels",nan(length(labels),1));
% set(gca,"XDisplayLabels",nan(length(nodes_adaptive)));
set(gca,"YData",snm.ID(idx_ecm2));
set(gca,"ColorLimits",[-0.7 0.7]);
set(gca,"Colormap",map);
% xlabel("Drug X + Drug Y")
ylabel("Model Output");
annotation('textbox',[.825 .95 .18 .04],'String','\DeltaActivity','EdgeColor','none');
title("Low Tension");

% figure;
subplot("Position",[0.25 0.3 0.6 0.29]);
heatmap(act_ht_adaptive_plot.');
% set(gca,"ColorbarVisible","off");
set(gca,"XDisplayLabels",labels);
set(gca,"YData",snm.ID(idx_ecm2));
% set(gca,"YDisplayLabels",nan(length(idx_ecm),1));
set(gca,"ColorLimits",[-0.7 0.7]);
set(gca,"Colormap",map);
xlabel("Drug Combination");
ylabel("Model Output");
% annotation('textbox',[.85 .83 .3 .1],'String','\DeltaActivity','EdgeColor','none');
title("High Tension");
% sgtitle("Mechano-Adaptive Perturbations");



% changes in output activity for all perturbs
figure("Position",[250 460 560 330]);
subplot(6,1,1:5);
heatmap(transpose(act_lt(idx_threshold,idx_ecm2)));     % plot all outputs for perturbations
set(gca,"GridVisible",0);
set(gca,"XDisplayLabels",nan(size(idx_threshold,1),1));
set(gca,"YData",snm.ID(idx_ecm2));
set(gca,"ColorLimits",[-1 1]);
set(gca,"Colormap",map);
title("Low Tension");
subplot(6,1,6);
heatmap(transpose(lt_mcc(idx_threshold)));
set(gca,"GridVisible",0);
set(gca,"XDisplayLabels",nan(size(idx_threshold,1),1));
set(gca,"YData","\DeltaECM");
set(gca,"ColorLimits",[-0.7 0.7]);
set(gca,"Colormap",map);

figure("Position",[250 50 560 330]);
subplot(6,1,1:5);
heatmap(transpose(act_ht(idx_threshold,idx_ecm2)));     % plot all outputs for perturbations
set(gca,"GridVisible",0);
set(gca,"XDisplayLabels",nan(size(idx_threshold,1),1));
set(gca,"YData",snm.ID(idx_ecm2));
set(gca,"ColorLimits",[-0.7 0.7]);
set(gca,"Colormap",map);
title("High Tension");
subplot(6,1,6);
heatmap(transpose(ht_mcc(idx_threshold)));
set(gca,"GridVisible",0);
set(gca,"XDisplayLabels",nan(size(idx_threshold,1),1));
set(gca,"YData","\DeltaECM");
set(gca,"ColorLimits",[-1 1]);
set(gca,"Colormap",map);



% Calculate changes in matrix content (06.17.2020 JR)
% imports 'act' as mxn matrix (perturbs x nodes)
% uses 'pro' and 'anti' doubles to delineate pro-fibrotic + anti-fibrotic
function act_diff = calculateMCC(act,idx_pro,idx_anti,idx_inhib)
act_pro = act(:,idx_pro);
act_anti = act(:,idx_anti);
if isempty(idx_inhib) == 0
    act_inhib = act(:,idx_inhib);
    act_diff_anti = mean(act_anti,2) - mean(act_inhib,2);
else
    act_diff_anti = sum(act_anti,2);
end
act_diff = sum(act_pro,2) - act_diff_anti;
end

function [meas_pro,meas_anti,meas_inhib] = calculatePAI(act,idx_pro,idx_anti,idx_inhib)
act_pro = act(:,idx_pro);
act_anti = act(:,idx_anti);
act_inhib = act(:,idx_inhib);

meas_pro = mean(act_pro,2);
meas_anti = mean(act_anti,2);
meas_inhib = mean(act_inhib,2);
end

function score_idx = scorePAI(meas_pro,meas_anti,meas_inhib,type)
if type == "ht"
    dirs = ["descend";"ascend";"descend"];
elseif type == "lt"
    dirs = ["ascend";"descend";"ascend"];
else
    fprintf("Setting score type to 'lt'\n")
    dirs = ["ascend";"descend";"ascend"];
end
[~,sort_pro_idx] = sort(meas_pro,dirs(1));          % sort arrays
[~,sort_anti_idx] = sort(meas_anti,dirs(2));
[~,sort_inhib_idx] = sort(meas_inhib,dirs(3));
sort_pro_idx(:,2) = 1:length(sort_pro_idx);         % add ranks
sort_anti_idx(:,2) = 1:length(sort_anti_idx);
sort_inhib_idx(:,2) = 1:length(sort_inhib_idx);
rank_pro_idx = sortrows(sort_pro_idx,1);            % reorder to match lt/ht rows
rank_anti_idx = sortrows(sort_anti_idx,1);
rank_inhib_idx = sortrows(sort_inhib_idx,1);

scorefun = @(x,y,z) length(x)-mean([x y z],2);
scores = scorefun(rank_pro_idx(:,2),rank_anti_idx(:,2),rank_inhib_idx(:,2));
score_idx = [rank_pro_idx(:,1) scores];             % sum lt/ht ranks
end

function score_idx = scoreAll(lt_pro,lt_anti,lt_inhib,ht_pro,ht_anti,ht_inhib)
lt_dirs = ["ascend";"descend";"ascend"];
ht_dirs = ["descend";"ascend";"descend"];
[~,sort_lt_pro] = sort(lt_pro,lt_dirs(1));          % sort arrays
[~,sort_lt_anti] = sort(lt_anti,lt_dirs(2));
[~,sort_lt_inhib] = sort(lt_inhib,lt_dirs(3));
[~,sort_ht_pro] = sort(ht_pro,ht_dirs(1));          % sort arrays
[~,sort_ht_anti] = sort(ht_anti,ht_dirs(2));
[~,sort_ht_inhib] = sort(ht_inhib,ht_dirs(3));
sort_lt_pro(:,2) = 1:length(sort_lt_pro);         % add ranks
sort_lt_anti(:,2) = 1:length(sort_lt_anti);
sort_lt_inhib(:,2) = 1:length(sort_lt_inhib);
sort_ht_pro(:,2) = 1:length(sort_ht_pro);         % add ranks
sort_ht_anti(:,2) = 1:length(sort_ht_anti);
sort_ht_inhib(:,2) = 1:length(sort_ht_inhib);
rank_lt_pro = sortrows(sort_lt_pro,1);            % reorder to match lt/ht rows
rank_lt_anti = sortrows(sort_lt_anti,1);
rank_lt_inhib = sortrows(sort_lt_inhib,1);
rank_ht_pro = sortrows(sort_ht_pro,1);            % reorder to match lt/ht rows
rank_ht_anti = sortrows(sort_ht_anti,1);
rank_ht_inhib = sortrows(sort_ht_inhib,1);

scorefun = @(x,y,z,a,b,c) length(x)-mean([x y z a b c],2);
scores = scorefun(rank_lt_pro(:,2),rank_lt_anti(:,2),rank_lt_inhib(:,2),rank_ht_pro(:,2),rank_ht_anti(:,2),rank_ht_inhib(:,2));
score_idx = [rank_lt_pro(:,1) scores];             % sum lt/ht ranks

end