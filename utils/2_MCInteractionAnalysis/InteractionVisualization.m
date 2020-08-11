%Compilaton/heatmap generation for simulation data
%Written by Jesse Rogers - 29 June 2017
%Edited by Jesse Rogers - 16 Oct 2019: tailored for SNM 1.1rev3

%% dAA matrix loading
%Matrix filename should be formatted with input name:
%'dSens_(input).mat'
%Line 7: Select number of input data sets to be used
imax = 9;
inputs = cell(1,imax);      %labels for x-axis (1 x 10 cell) (input nodes changed)
inputNames = {'AngII','TGFB','IL6','IL1','TNFa','NE','PDGF','ET1','NP'};
speciesNames = {'AngII','AT1R','AGT','ACE','NOX','ROS','ET1','ETAR','DAG','PKC','TRPC','NE','BAR','AC','cAMP','PKA','CREB','CBP','TGFB','TGFB1R','smad3','smad7','latentTGFB','BAMBI','PDGF','PDGFR','NP','NPRA','cGMP','PKG','tension','B1int','Rho','ROCK','Ca','calcineurin','NFAT','IL6','gp130','STAT','IL1','IL1RI','TNFa','TNFaR','NFKB','PI3K','Akt','p38','TRAF','ASK1','MKK3','PP1','JNK','abl','Rac1','MEKK1','MKK4','ERK','Ras','Raf','MEK1','FAK','epac','Factin','FA','cmyc','CTGF','proliferation','SRF','EDAFN','aSMA','AP1','TIMP1','TIMP2','PAI1','proMMP14','proMMP1','proMMP2','proMMP9','fibronectin','periostin','proCI','proCIII','B3int','Src','Grb2','p130Cas','YAP','MRTF','Gactin','TNC','mTORC1','mTORC2','p70S6K','EBP1','syndecan4','proMMP3','proMMP8','proMMP12','thrombospondin4','osteopontin','contractility','RhoGEF','RhoGDI','talin','vinculin','paxillin','MLC','AT2R',}; 
% dAAlo = zeros(length(speciesNames),imax);     %values for heatmap 1 
% dAAhi = zeros(length(speciesNames),imax);     %values for heatmap 2 
% rev = zeros(length(speciesNames),imax);       %reversal values for heatmap 2 

%opening datasets for each input
datadir = 'data\2_MCInteractionAnalysis\';
plotdir = 'plots\';
cmapdir = 'utils\BrewerMap-master\';
addpath(cmapdir);
map = brewermap(21,'PRGn');
map_rev = [1 0.47 0.24];

for k = 1:imax
    choice = inputNames{k};
    filename_dAUC = strcat(datadir,'dAUC_',choice,'.mat');
    filename_rev = strcat(datadir,'rev',choice,'.mat');
    load(filename_dAUC)
    load(filename_rev)
    % dAUClo(:,k) = dAUC(:,:,1).';
    dAUCplot(:,k,:) = permute(dAUC,[2 1 3]);
    for i = 1:size(dAUCplot,3)
        for j = 1:length(speciesNames)
            if rev(:,j,i) == -1 || rev(:,j,i) == 1
                dAUCplot(j,k,i) = NaN;                   % reassign reversed values
            end
        end
    end
    inputs(1,k) = cellstr(choice);
end

%% Heatmap Generation;

%setting colorbar limits (based on max/min values of each dataset)
% lolimit = max(dAAplot,[],'all');
limit = min(dAUCplot,[],'all');
reverselabel = "Reverse";
tensionlabel = 0.1:0.1:0.9;

%%%%%% Figure 3A-B: dAUC heatmaps at tension levels 0.5/0.9 %%%%%%
for i = [5 9]
% for i = 1:size(dAUCplot,2)
    pos = 250;
    if any(isnan(dAUCplot(:,:,i)),'all')
    % if any(isnan(dAUCplot(:,i,:)),'all')
        pos = pos+23.08;                   % adjust width for reverse-containing cases
    end
    fig = figure("position",[100+50*i 100 pos 400]);
    % subplot(1,2,1);
    fig1 = heatmap(fig,inputNames,speciesNames,dAUCplot(:,:,i));
    % fig1 = heatmap(fig,tensionlabel,speciesNames,permute(dAUCplot(:,i,:),[1 3 2]));
    grid off
    fig1.Title = strcat("Tension = ",string(0.1+0.1*(i-1)));
    % fig1.Title = inputs{i};
    fig1.YLabel = 'Measured Node';
    fig1.Colormap = map;
    fig1.ColorLimits = [limit -limit];
    fig1.YDisplayLabels = nan(length(fig1.YDisplayLabels),1);   % hide y labels
    % fig1.InnerPosition = [0.1 0.1 0.56 0.8];                  % for subplot
    fig1.MissingDataLabel = reverselabel;
    fig1.MissingDataColor = map_rev;
    if i == size(dAUCplot,3)
    % if any(isnan(dAUCplot(:,i,:)),'all')
        annotation('textbox',[.71 .88 .3 .1],'String','\DeltaAUC','EdgeColor','none','FontSize',8);
    end
    % fig1.FontSize = 8;
    % fig1.ColorbarVisible = 'off';
    filename_plot = strcat('heatmap_dAUC_tension',replace(string(0.1+0.1*(i-1)),".",""));
    % filename_plot = strcat('heatmap_fback075_input',inputs{i});
    saveas(fig1,strcat(plotdir,filename_plot,'.fig'))
    % saveas(fig1,strcat(plotdir,filename_plot,'.svg'))
end

% Distribution plots: between inputs
% figure("Position",[80 80 350 700]);
% map2 = brewermap(9,'Spectral');
% for i = [2 5 9]
%     % subplot(5,2,i);
%     figure("Position",[100+20*i 300 400 300]);
%     inputlvls = 1:9;
%     for j = 1:length(inputlvls)
%         [kde,xi] = ksdensity(dAUCplot(:,inputlvls(j),i),"Bandwidth",0.01);
%         % [counts,bins] = histcounts(dAUCplot(:,j,i), ...
%         %     "Normalization","pdf", ...
%         %     "BinWidth",0.05, ...
%         %     "BinLimits",[-0.6 0.6]);
%         % plot(bins(2:end),counts);
%         plot(xi,kde,"Color",map2(j,:),"LineWidth",1.5);
%         % ylim([0 1]);
%         if j==1
%             hold on
%         end
%     end
%     xline(-0.05,'--'); xline(0.05,'--');
%     annotation('arrow',[0.4 0.2],[0.8 0.8],"Color",map(5,:), ...
%         "HeadLength",5,"HeadWidth",5,"HeadStyle","plain","LineWidth",1);
%     annotation('arrow',[0.63 0.83],[0.8 0.8],"Color",map(16,:), ...
%         "HeadLength",5,"HeadWidth",5,"HeadStyle","plain","LineWidth",1);
%     annotation('textbox',[0.2 0.76 0.15 0.13],'String','Dampening', ...
%         'EdgeColor','none','Color',map(5,:));
%     annotation('textbox',[0.62 0.76 0.15 0.13],'String','Amplification', ...
%         'EdgeColor','none','Color',map(16,:));
%     xlim([-0.6 0.6]); xticks(-0.6:0.2:0.6);
%     xlabel("\DeltaAUC"); ylabel("Probability Density");
%     title(strcat("Tension = ",string(0.1+0.1*(i-1))));
%     legend(inputNames,"Location","southoutside","NumColumns",3);
%     hold off
% end

% Distribution plots: between tension
map2 = brewermap(size(dAUCplot,3),'Blues');
% map2 = [0.07 0.62 1;1 0.41 0.16];
tensionlvls = [2:9];
% figure("Position",[160 80 350 700]);
for i = 1:size(dAUCplot,2)
    figure("Position",[100+20*i 300 400 200]);
    for j = 1:length(tensionlvls)
        % [counts,bins] = histcounts(dAAplot(:,i,j), ...
        %     "Normalization","pdf", ...
        %     "BinWidth",0.05, ...
        %     "BinLimits",[-0.6 0.6]);
        % plot(bins(2:end),counts,"Color",map2(j,:));
        % ylim([0 1]);
        [kde,xi] = ksdensity(dAUCplot(:,i,tensionlvls(j)),"Bandwidth",0.01);
        plot(xi,kde,"LineWidth",1.5,"Color",map2(j+1,:));
        if j==1
            hold on
        end
    end
    xline(-0.05,'--'); xline(0.05,'--');
    annotation('arrow',[0.35 0.15],[0.75 0.75], ...
        "HeadLength",5,"HeadWidth",5,"HeadStyle","plain","LineWidth",1);
    annotation('arrow',[0.5 0.7],[0.75 0.75], ...
        "HeadLength",5,"HeadWidth",5,"HeadStyle","plain","LineWidth",1);
    annotation('textbox',[0.15 0.76 0.15 0.13],'String','Dampening', ...
        'EdgeColor','none');
    annotation('textbox',[0.49 0.76 0.15 0.13],'String','Amplification', ...
        'EdgeColor','none');
    xlim([-0.6 0.6]); xticks(-0.6:0.2:0.6);
    xlabel("\DeltaAUC"); ylabel("Probability Density");
    set(gca,"Colormap",map2);
    title(inputNames{i});
    c = colorbar;
    c.Label.String = "Tension";
    %legend(["Tension = 0.3";"Tension = 0.5"], ...
    %    "Location","southoutside","NumColumns",2);
    hold off
    filename_plot = strcat('dist_dAUC_',inputNames{i});
    saveas(gca,strcat(plotdir,filename_plot,'.fig'))
end

% K-S tests + B-H correction: between tension
p_inputs_25 = [];
p_inputs_29 = [];
p_inputs_59 = [];
for row = 1:size(dAAplot,2)
    [~,p25] = kstest2(dAAplot(:,row,2),dAAplot(:,row,5));   % compare tension 0.2/0.5
    [~,p29] = kstest2(dAAplot(:,row,2),dAAplot(:,row,9));   % compare tension 0.2/0.9
    [~,p59] = kstest2(dAAplot(:,row,5),dAAplot(:,row,9));   % compare tension 0.2/0.9
    p_inputs_25 = [p_inputs_25;p25];
    p_inputs_29 = [p_inputs_29;p29];
    p_inputs_59 = [p_inputs_59;p59];
end
p_all = [p_inputs_25;p_inputs_29;p_inputs_59];
fdr_all = mafdr(p_all,"BHFDR","true");
fwer_all = p_all*length(p_all);
len_p = length(p_inputs_25);
fdr_inputs_25 = fdr_all(1:len_p);
fdr_inputs_29 = fdr_all(len_p+1:2*len_p);
fdr_inputs_59 = fdr_all(2*len_p+1:end);
