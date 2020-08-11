function PerturbationAnalysis(tension)
% 04-16-2020 JR
% Analysis/Visualization of sensitivity simulation results

% load data
filedir = "data\3_PerturbationAnalysis\";
tensionname = replace(string(tension),".","");
load(strcat(filedir,'act_tension',tensionname,'.mat'),'act_delta');
load(strcat(filedir,'speciesNames.mat'),'speciesNames');

%%%%%% Figure S4: network-wide heatmap %%%%%%
plotdir = "plots\";

pos = 400;
if tension == 0.75
    pos = pos+32;
end
fig = figure('position',[100 100 pos 300]);
fig_delta = heatmap(fig,act_delta);
grid off
fig_delta.ColorLimits = [-1 1];
fig_delta.Colormap = redbluecmap();
fig_delta.XDisplayLabels = nan(1,length(fig_delta.XDisplayLabels));
fig_delta.YDisplayLabels = nan(length(fig_delta.YDisplayLabels),1);
if tension ~= 0.75
    fig_delta.ColorbarVisible = 0;
else
    annotation('textbox',[.81 .9 .3 .1],'String','\DeltaActivity','EdgeColor','none');
end
% fig_delta.FontSize = 7;
xlabel("Measured Node");
ylabel("Knocked-Down Node");
title(strcat("Tension = ",string(tension)));
plotfilename = strcat(plotdir,'S4_AllPerturbations_tension',tensionname,'.fig');
saveas(fig_delta,plotfilename)
% plotfilename = strcat(plotdir,'S4_AllPerturbations_tension',tensionname,'.svg');
% saveas(fig_delta,plotfilename)

% % output-only heatmap
% outputs_idx = [82,83,80,101,100,81,91,70,77,78,97,98,79,99,76,73,74,75,71,67,23];
% pos = 300;
% if tension == 0.75
%     pos = pos+24;
% end
% fig = figure('position',[100 100 pos 400]);
% fig_delta_out = heatmap(fig,speciesNames(outputs_idx),speciesNames,act_delta(:,outputs_idx));
% grid off
% fig_delta_out.ColorLimits = [-1 1];
% fig_delta_out.Colormap = redbluecmap();
% fig_delta_out.YDisplayLabels = nan(length(fig_delta_out.YDisplayLabels),1);
% if tension ~= 0.75
%     fig_delta_out.ColorbarVisible = 0;
% end
% % fig_delta_out.FontSize = 7;
% xlabel("Measured Node");
% ylabel("Knocked-Down Node");
% title(strcat("Tension = ",string(tension)));
% plotfilename = strcat(plotdir,'act_delta_outputs_tension_nolabs',tensionname,'_final.fig');
% saveas(fig_delta_out,plotfilename)
% plotfilename = strcat(plotdir,'act_delta_outputs_tension_nolabs',tensionname,'_final.svg');
% saveas(fig_delta_out,plotfilename)


% top 10 filtering: remove tension from results (07.29.2020 JR)
istension = strcmp(speciesNames,"tension");
speciesNames_notens = speciesNames(~istension);

% Calculate knockdown sensitivity + rank nodes
sens = sum(abs(act_delta),1);               % calculate sensitivity
sens_notens = sens(~istension);
[sens_sort,sens_sort_idx] = sort(sens_notens,'descend');
sens_top = sens_sort(1:10);
sens_top_idx = sens_sort_idx(1:10);

% % sensitivity barplot: top 10
% fig = figure('position',[100 100 400 300]);
% fig_sens = barh(flip(sens_top));
% fig_sens.FaceColor = '#4DBEEE';
% yticklabels(speciesNames_notens(flip(sens_top_idx)));
% xlabel("Total Sensitivity");
% title(strcat("Tension = ",string(tension)));
% plotfilename = strcat(plotdir,'sens_all_tension',tensionname,'_final.fig');
% saveas(fig_sens,plotfilename)
% plotfilename = strcat(plotdir,'sens_all_tension',tensionname,'_final.svg');
% saveas(fig_sens,plotfilename)

% Calculate knockdown influence + rank nodes
infl = sum(abs(act_delta),2);               % calculate influence
infl_notens = infl(~istension);
[infl_sort,infl_sort_idx] = sort(infl_notens,'descend');
infl_top = infl_sort(1:10);
infl_top_idx = infl_sort_idx(1:10);

% % influence barplot: top 10
% 
% fig = figure('position',[100 100 400 300]);
% fig_infl = barh(flip(infl_top));
% fig_infl.FaceColor = '#A2142F';
% yticklabels(speciesNames_notens(flip(infl_top_idx)));
% xlabel("Total Influence");
% title(strcat("Tension = ",string(tension)));
% plotfilename = strcat(plotdir,'infl_all_tension',tensionname,'_final.fig');
% saveas(fig_infl,plotfilename)
% plotfilename = strcat(plotdir,'infl_all_tension',tensionname,'_final.svg');
% saveas(fig_infl,plotfilename)


% export top sensitive/influential nodes as lists
save(strcat(filedir,'toplists_tension',tensionname,'.mat'), ...
    'sens_notens','sens_top_idx','infl_notens','infl_top_idx','speciesNames_notens');

end