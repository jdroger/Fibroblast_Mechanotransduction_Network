function snm_1_1_rev4_sensitivitytoprankings
% 04-16-2020 JR
% Analysis/Visualization of sensitivity simulation results

% load/compile data
cmapdir = "D:/Research/Aim2/BrewerMap-master/";
filedir = "D:\Research\Aim2\ModelExpansion\1_1\rev4\SensitivityAnalysis_Zeigler\sensitivity_data\";
filename = 'snm_1_1_rev4_sensitivity_';
% load(strcat(filedir,filename,'speciesNames','_final.mat'),'speciesNames');
addpath(cmapdir);

tensioncases = [0.25 0.5 0.75];
sens_allcases = zeros(108,length(tensioncases));
sens_idx_allcases = zeros(10,length(tensioncases));
infl_allcases = zeros(108,length(tensioncases));
infl_idx_allcases = zeros(10,length(tensioncases));

for tension = 1:length(tensioncases)
    tensionname = replace(string(tensioncases(tension)),".","");
    load(strcat(filedir,filename,'toplists','_tension',tensionname,'_final.mat'), ...
        'sens_notens','sens_top_idx','infl_notens','infl_top_idx','speciesNames_notens');
    % compile separate vectors into matrices
    sens_allcases(:,tension) = sens_notens;
    sens_idx_allcases(:,tension) = sens_top_idx;
    infl_allcases(:,tension) = infl_notens;
    infl_idx_allcases(:,tension) = infl_top_idx;
end

% find unique nodes across all cases
sens_idx_unique = unique(sens_idx_allcases,'stable');
infl_idx_unique = unique(infl_idx_allcases,'stable');

% get values for unique nodes in each case
sens_allcases_unique = sens_allcases(sens_idx_unique,:);
infl_allcases_unique = infl_allcases(infl_idx_unique,:);

% set flags for top nodes in each stretch case
sens_flag = setFlags(sens_idx_unique,sens_idx_allcases);
infl_flag = setFlags(infl_idx_unique,infl_idx_allcases);

% plot values for each case
plotdir = "D:\Research\Aim2\ModelExpansion\1_1\rev4\SensitivityAnalysis_Zeigler\sensitivity_plots\";
map = brewermap(3,'Set2');

fig = figure('position',[100 100 650 450]);
for tension = 1:length(tensioncases)
    subplot(1,length(tensioncases),tension);
    ax = barh(flip(sens_allcases_unique(:,tension)));
    yticks(1:1:length(sens_idx_unique));
    xlim([0 max(sens_allcases_unique,[],'all')+0.5]);
    if tension == 1
        yticklabels(speciesNames_notens(flip(sens_idx_unique)));
    else
        yticklabels("");
    end
    % set colors based on flags
    % colors = get(gca,"ColorOrder");
    ax.FaceColor = 'flat';
    for c = 1:3
        for f = 1:length(tensioncases)
            ax.CData(flip(sens_flag)==f,c) = map(f,c); 
        end
    end
    xlabel("Knockdown Sensitivity");
    title(strcat("Tension = ",string(tensioncases(tension))));
    pos = get(gca,"Position");
    pos(2) = pos(2)+0.02;
    set(gca,"FontSize",10,"XGrid","on","Position",pos);
end
plotfilename = strcat(plotdir,'sens_allcases_final.fig');
saveas(fig,plotfilename)
plotfilename = strcat(plotdir,'sens_allcases_final.svg');
saveas(fig,plotfilename)

fig = figure('position',[100 100 650 350]);
for tension = 1:length(tensioncases)
    subplot(1,length(tensioncases),tension);
    ax = barh(flip(infl_allcases_unique(:,tension)));
    yticks(1:1:length(infl_idx_unique));
    xlim([0 max(infl_allcases_unique,[],'all')+0.5]);
    if tension == 1
        yticklabels(speciesNames_notens(flip(infl_idx_unique)));
    else
        yticklabels("");
    end
    % set colors based on flags
    % colors = get(gca,"ColorOrder");
    ax.FaceColor = 'flat';
    for c = 1:3
        for f = 1:length(tensioncases)
            ax.CData(flip(infl_flag)==f,c) = map(f,c); 
        end
    end
    xlabel("Knockdown Influence");
    title(strcat("Tension = ",string(tensioncases(tension))));
    pos = get(gca,"Position");
    pos(2) = pos(2)+0.02;
    set(gca,"FontSize",10,"XGrid","on","Position",pos);
end
plotfilename = strcat(plotdir,'infl_allcases_final.fig');
saveas(fig,plotfilename)
plotfilename = strcat(plotdir,'infl_allcases_final.svg');
saveas(fig,plotfilename)

fig = figure('position',[100 100 1000 300]);
ax = bar(sens_allcases_unique);
xticks(1:1:length(sens_idx_unique));
xticklabels(speciesNames_notens(sens_idx_unique));
xtickangle(45);
ylabel("Total Sensitivity");
legend("Tension = 0.25","Tension = 0.5","Tension = 0.75","Location","northeastoutside");
plotfilename = strcat(plotdir,'sens_allcases_horiz_final.fig');
saveas(fig,plotfilename)

fig = figure('position',[100 100 1000 300]);
ax = bar(infl_allcases_unique);
xticks(1:1:length(infl_idx_unique));
xticklabels(speciesNames_notens(infl_idx_unique));
xtickangle(45);
ylabel("Total Influence");
legend("Tension = 0.25","Tension = 0.5","Tension = 0.75","Location","northeastoutside");
plotfilename = strcat(plotdir,'infl_allcases_horiz_final.fig');
saveas(fig,plotfilename)


% Influencer-Sensitive Node Heatmaps
map2 = brewermap(19,'*RdBu');
idx_tension = 31;
% load act_delta data
act_allcases = zeros(10,10,length(tensioncases));
for tension = 1:length(tensioncases)
    tensionname = replace(string(tensioncases(tension)),".","");
    load(strcat(filedir,filename,'act_delta','_tension',tensionname,'_final.mat'), ...
        'act_delta');
    % remove tension node (07.29.2020 JR)
    act_delta(idx_tension,:) = [];
    act_delta(:,idx_tension) = [];
    % subset + compile matrices (for top 10 cases)
    act_allcases(:,:,tension) = act_delta(infl_idx_allcases(:,tension),sens_idx_allcases(:,tension));
end

% plot heatmaps
for tension = 1:length(tensioncases)
    pos = 400;
    if tension == length(tensioncases)
        pos = pos+32;
    end
    fig = figure('position',[100 100 pos 350]);
    fig_delta = heatmap(fig, ...
        speciesNames_notens(sens_idx_allcases(:,tension)), ...
        speciesNames_notens(infl_idx_allcases(:,tension)), ...
        act_allcases(:,:,tension));
    fig_delta.ColorLimits = [-1 1];
    fig_delta.Colormap = map2;
    if tension ~= length(tensioncases)
        fig_delta.ColorbarVisible = 0;
    else
        annotation('textbox',[.81 .9 .3 .1],'String','\DeltaActivity','EdgeColor','none');
    end
    xlabel("Measured Node");
    ylabel("Knocked-Down Node");
    title(strcat("Tension = ",string(tensioncases(tension))));
    set(gca,"FontSize",12);
    tensionname = replace(string(tensioncases(tension)),".","");
    plotfilename = strcat(plotdir,'act_delta_tension',tensionname,'_top_final.fig');
    saveas(fig_delta,plotfilename)
    plotfilename = strcat(plotdir,'act_delta_tension',tensionname,'_top_final.svg');
    saveas(fig_delta,plotfilename)
end


end

function flag = setFlags(unique,allcases)
flag = zeros(length(unique),1);
for node = 1:length(unique)
    inSet = (allcases == unique(node));
    if any(inSet(:,1) == 1)
        flag(node) = 1;
    elseif any(inSet(:,2) == 1)
        flag(node) = 2;
    elseif any(inSet(:,3) == 1)
        flag(node) = 3;
    end
end
end