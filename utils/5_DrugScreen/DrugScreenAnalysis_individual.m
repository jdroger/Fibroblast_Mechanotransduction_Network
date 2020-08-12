% Drug Screen: results visualization
% 01.12.2019 JR
% 
% Purpose: combine low-/high-tension screens to identify perturbations that
% decrease matrix content in low tension and increase matrix content in
% high tension

%% Import + Process Data
homedir = "data\";
sensdir = "utils\5_DrugScreen";
datadir = "5_DrugScreen\";
cmapdir = "utils\BrewerMap-master\";
plotdir = "plots\";
addpath(homedir);
addpath(sensdir);
addpath(strcat(homedir,datadir));
addpath(cmapdir);

modelpath = 'snm_1_1_rev4.xlsx';
opts = detectImportOptions(modelpath);
snm = readtable(modelpath, opts);

% method 1: use all outputs, counting MMPs as anti-fibrotic and rest as
% pro-fibrotic
idx_ecm = find(strcmpi(snm.module,"ECM"));
idx_ecm_pro = idx_ecm(~contains(snm.ID(idx_ecm),"MMP"));
idx_ecm_anti = idx_ecm(contains(snm.ID(idx_ecm),"MMP"));
idx_individual = 1:109;
% method 2: limit outputs to those with direct remodeling capability
idx_ecm_pro2 = idx_ecm([3:5 10:13 19]);
idx_ecm_anti2 = idx_ecm([6:9 15:17]);
idx_ecm2 = [idx_ecm_pro2;idx_ecm_anti2];
% method 3: separate proteases/inhibs and add mutual inhibition
idx_ecm_pro3 = idx_ecm([10:13 19]);
idx_ecm_inhib = idx_ecm(3:5);


% remove perturbations modulating outputs
remove_values = ismember(transpose(1:length(snm.ID)),idx_ecm);
% remove_rows = (remove_values(:,1) == 1 | remove_values(:,2) == 1);
nodes_all = [idx_individual(~remove_values).';idx_individual(~remove_values).'];
len_c = length(idx_individual(~remove_values));

% import/filter data
makeFileName = @(x,y) strcat("individual_tension",replace(string(x),".",""),"_",y,".mat");

% low-tension data
tension = 0.1; 
choice = "KO";
act_lt_ko = load(makeFileName(tension,choice));
choice = "OE";
act_lt_oe = load(makeFileName(tension,choice));
act_lt_ecm_ko = act_lt_ko.act_delta(~remove_values,idx_ecm);
act_lt_ecm_oe = act_lt_oe.act_delta(~remove_values,idx_ecm);
act_lt = [act_lt_ko.act_delta(~remove_values,:);act_lt_oe.act_delta(~remove_values,:)];

% high-tension data
tension = 0.6;
choice = "KO";
act_ht_ko = load(makeFileName(tension,choice));
choice = "OE";
act_ht_oe = load(makeFileName(tension,choice));
act_ht_ecm_ko = act_ht_ko.act_delta(~remove_values,idx_ecm);
act_ht_ecm_oe = act_ht_oe.act_delta(~remove_values,idx_ecm);
act_ht = [act_ht_ko.act_delta(~remove_values,:);act_ht_oe.act_delta(~remove_values,:)];

%% Calculate overall matrix content change
% Method 1: calc. matrix content change metric from
% pro/anti/inhib outputs
lt_mcc = calculateMCC(act_lt,idx_ecm_pro2,idx_ecm_anti2,[]);
ht_mcc = calculateMCC(act_ht,idx_ecm_pro2,idx_ecm_anti2,[]);
% Method 2: rank matrix content changes in pro/anti/inhib categories
% [lt_pro,lt_anti,lt_inhib] = calculatePAI(act_lt,idx_ecm_pro3,idx_ecm_anti2,idx_ecm_inhib);
% [ht_pro,ht_anti,ht_inhib] = calculatePAI(act_ht,idx_ecm_pro3,idx_ecm_anti2,idx_ecm_inhib);
% lt_scores = scorePAI(lt_pro,lt_anti,lt_inhib,"lt");
% ht_scores = scorePAI(ht_pro,ht_anti,ht_inhib,"ht");


% Find mechano-adaptive perturbatons
% Method 1: ID perturbs from calculated MCC using threshold
threshold = 0;
idx_threshold = find(lt_mcc <= -threshold & ht_mcc >= threshold); 

% Method 2: rank perturbs based on MCC, find top-ranked perturbs
[~,sort_lt_idx] = sort(lt_mcc(idx_threshold),"ascend");     % sort arrays
[~,sort_ht_idx] = sort(ht_mcc(idx_threshold),"descend");
sort_lt_idx(:,2) = 1:length(sort_lt_idx);                   % add ranks
sort_ht_idx(:,2) = 1:length(sort_ht_idx);
rank_lt_idx = sortrows(sort_lt_idx,1);                      % reorder to match lt/ht rows
rank_ht_idx = sortrows(sort_ht_idx,1);
score = @(x,y) length([x;y])-(x+y);
score_idx = [rank_lt_idx(:,1) score(rank_lt_idx(:,2),rank_ht_idx(:,2))];  % sum lt/ht ranks
sort_score_idx = sortrows(score_idx,2,"descend");           % sort by score
idx_adaptive = idx_threshold(sort_score_idx(:,1)).';

% Method 3: rank perturbs based on pro/anti/inhib scores independently,
% combine rankings for lt/ht and find top-ranked perturbs
% score_both = [lt_scores(:,1) mean([lt_scores(:,2) ht_scores(:,2)],2)];
% sort_score_both = sortrows(score_both,2,"descend");
% idx_adaptive_3 = sort_score_both(1:50,1).';

% Method 3b: combine pro/anti/inhib rankings for lt/ht together, and find
% top-ranked perturbs
% score_both2 = scoreAll(lt_pro,lt_anti,lt_inhib,ht_pro,ht_anti,ht_inhib);
% sort_score_both2 = sortrows(score_both2,2,"descend");
% idx_adaptive_3b = sort_score_both(1:50,1).';

% Find perturb targets/types for each hit
nodes_adaptive = nodes_all(idx_adaptive,:);     % filter node combinations
nodes_type = string(zeros(length(idx_adaptive),1));

for node = 1:length(idx_adaptive)
    if idx_adaptive(node)>len_c
        nodes_type(node) = "OE";
    else
        nodes_type(node) = "KO";
    end
end
ID_adaptive = snm.ID(nodes_adaptive);      % translate to node IDs
changes_adaptive = lt_mcc(idx_adaptive);   % filter matrix changes
changes_adaptive(:,2) = ht_mcc(idx_adaptive);

% display data as a table
table_adaptive = cell2table(ID_adaptive);
table_adaptive(:,2:3) = array2table(changes_adaptive);
table_adaptive.Properties.VariableNames = ["Node","LowTension","HighTension"];
table_adaptive.perturbation = nodes_type;

%% Plot Results
%%%%%% Figure S6: Mechano-adaptive scoring of perturbaions %%%%%%
xdata = score_idx(:,2);
figure('position',[100 100 600 300]);
fig_lt = scatter(xdata,lt_mcc(idx_threshold),10); hold on
fig_ht = scatter(xdata,ht_mcc(idx_threshold),10);
set(gca,"FontSize",12,"YLim",[-1.5 1.5]);
% xticks(0:100:1000);
legend('Low Tension','High Tension','Location','southwest');
xlabel('Mechano-Adaptive Score'); ylabel('Matrix Content Change');
% annotation('rectangle',[0.76 0.35 0.07 0.5], ...
%     'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--')
grid on
saveas(gcf,strcat(plotdir,"S6A_DrugScreen_individual.fig"));


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
    act_diff_anti = mean(act_anti,2);
end
act_diff = mean(act_pro,2) - act_diff_anti;
end

% function [meas_pro,meas_anti,meas_inhib] = calculatePAI(act,idx_pro,idx_anti,idx_inhib)
% act_pro = act(:,idx_pro);
% act_anti = act(:,idx_anti);
% act_inhib = act(:,idx_inhib);
% 
% meas_pro = mean(act_pro,2);
% meas_anti = mean(act_anti,2);
% meas_inhib = mean(act_inhib,2);
% end
% 
% function score_idx = scorePAI(meas_pro,meas_anti,meas_inhib,type)
% if type == "ht"
%     dirs = ["descend";"ascend";"descend"];
% elseif type == "lt"
%     dirs = ["ascend";"descend";"ascend"];
% else
%     fprintf("Setting score type to 'lt'\n")
%     dirs = ["ascend";"descend";"ascend"];
% end
% [~,sort_pro_idx] = sort(meas_pro,dirs(1));          % sort arrays
% [~,sort_anti_idx] = sort(meas_anti,dirs(2));
% [~,sort_inhib_idx] = sort(meas_inhib,dirs(3));
% sort_pro_idx(:,2) = 1:length(sort_pro_idx);         % add ranks
% sort_anti_idx(:,2) = 1:length(sort_anti_idx);
% sort_inhib_idx(:,2) = 1:length(sort_inhib_idx);
% rank_pro_idx = sortrows(sort_pro_idx,1);            % reorder to match lt/ht rows
% rank_anti_idx = sortrows(sort_anti_idx,1);
% rank_inhib_idx = sortrows(sort_inhib_idx,1);
% 
% scorefun = @(x,y,z) length(x)-mean([x y z],2);
% scores = scorefun(rank_pro_idx(:,2),rank_anti_idx(:,2),rank_inhib_idx(:,2));
% score_idx = [rank_pro_idx(:,1) scores];             % sum lt/ht ranks
% end
% 
% function score_idx = scoreAll(lt_pro,lt_anti,lt_inhib,ht_pro,ht_anti,ht_inhib)
% lt_dirs = ["ascend";"descend";"ascend"];
% ht_dirs = ["descend";"ascend";"descend"];
% [~,sort_lt_pro] = sort(lt_pro,lt_dirs(1));          % sort arrays
% [~,sort_lt_anti] = sort(lt_anti,lt_dirs(2));
% [~,sort_lt_inhib] = sort(lt_inhib,lt_dirs(3));
% [~,sort_ht_pro] = sort(ht_pro,ht_dirs(1));          % sort arrays
% [~,sort_ht_anti] = sort(ht_anti,ht_dirs(2));
% [~,sort_ht_inhib] = sort(ht_inhib,ht_dirs(3));
% sort_lt_pro(:,2) = 1:length(sort_lt_pro);         % add ranks
% sort_lt_anti(:,2) = 1:length(sort_lt_anti);
% sort_lt_inhib(:,2) = 1:length(sort_lt_inhib);
% sort_ht_pro(:,2) = 1:length(sort_ht_pro);         % add ranks
% sort_ht_anti(:,2) = 1:length(sort_ht_anti);
% sort_ht_inhib(:,2) = 1:length(sort_ht_inhib);
% rank_lt_pro = sortrows(sort_lt_pro,1);            % reorder to match lt/ht rows
% rank_lt_anti = sortrows(sort_lt_anti,1);
% rank_lt_inhib = sortrows(sort_lt_inhib,1);
% rank_ht_pro = sortrows(sort_ht_pro,1);            % reorder to match lt/ht rows
% rank_ht_anti = sortrows(sort_ht_anti,1);
% rank_ht_inhib = sortrows(sort_ht_inhib,1);
% 
% scorefun = @(x,y,z,a,b,c) length(x)-mean([x y z a b c],2);
% scores = scorefun(rank_lt_pro(:,2),rank_lt_anti(:,2),rank_lt_inhib(:,2),rank_ht_pro(:,2),rank_ht_anti(:,2),rank_ht_inhib(:,2));
% score_idx = [rank_lt_pro(:,1) scores];             % sum lt/ht ranks
% 
% end