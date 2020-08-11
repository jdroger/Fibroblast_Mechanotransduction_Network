%Analysis of Netflux simulation data
%Written by Jesse Rogers - Edited 29 June 2017
%% Activation matrix loading
%Matrix filename should be formatted with input name:
%'activation_(input).mat'
homedir = 'D:\Research\Aim2\ModelExpansion\1_1\rev4\SensitivityAnalysis_PrelimData\';
datadir = 'simulation_data\';
plotdir = 'simulation_plots\';
imax = 9;
inputNames = {'AngII','TGFB','IL6','IL1','TNFa','NE','PDGF','ET1','NP'};

for k = 1:imax
    choice = inputNames{k};
    filename = strcat(homedir,datadir,choice,'_final.mat');
    load(filename)
    speciesNames = {'AngII','AT1R','AGT','ACE','NOX','ROS','ET1','ETAR','DAG','PKC','TRPC','NE','BAR','AC','cAMP','PKA','CREB','CBP','TGFB','TGFB1R','smad3','smad7','latentTGFB','BAMBI','PDGF','PDGFR','NP','NPRA','cGMP','PKG','tension','B1int','Rho','ROCK','Ca','calcineurin','NFAT','IL6','gp130','STAT','IL1','IL1RI','TNFa','TNFaR','NFKB','PI3K','Akt','p38','TRAF','ASK1','MKK3','PP1','JNK','abl','Rac1','MEKK1','MKK4','ERK','Ras','Raf','MEK1','FAK','epac','Factin','FA','cmyc','CTGF','proliferation','SRF','EDAFN','aSMA','AP1','TIMP1','TIMP2','PAI1','proMMP14','proMMP1','proMMP2','proMMP9','fibronectin','periostin','proCI','proCIII','B3int','Src','Grb2','p130Cas','YAP','MRTF','Gactin','TNC','mTORC1','mTORC2','p70S6K','EBP1','syndecan4','proMMP3','proMMP8','proMMP12','thrombospondin4','osteopontin','contractility','RhoGEF','RhoGDI','talin','vinculin','paxillin','MLC','AT2R',}; 

    %% Data Normalization
    %Create new matrix by subtracting all response values from responses at
    %input of 0
    activation = real(activation);
    normalized = zeros(size(activation));
    doses = (0:0.01:1).';
    base = activation(1,:,:);
    for i = 1:size(activation,1)
        normalized(i,:,:) = activation(i,:,:) - base;
    end
    %% Normalized Plots
    % color = [''];
    fig = figure('Position',[100 100 500 300]);
    for j = 1:length(speciesNames)              %for proCI only: i=82
        legendlabels = cell(1,size(activation,3));
        for i = 1:size(activation,3)            %loop through tension vals
            %line 28: lines only, line 29: shaded lines
            ax = plot(doses(:),normalized(:,j,i),'LineWidth',2.0);
            % fig = area(doses(:),normalized(:,j,i),'LineWidth',1.0,'FaceColor',color(i,:),'EdgeColor',color(i,:),'FaceAlpha',0.5);
            hold on
            legendlabels{i} = strcat("Tension = ",string(0.1*i));
        end
        legend(legendlabels,'Location','northeastoutside');
        xlabel(strcat(choice,' Dose'));
        ylabel(strcat(speciesNames{j},' Normalized Activation'));
        plotfile = strcat(choice,'_',speciesNames{j},'_normalized');
        % title(speciesNames(j));
        hold off
        plotfilename = strcat(homedir,plotdir,choice,'\',plotfile,'_final.fig');
        saveas(ax,plotfilename)
    end
    %% Data Integration
    %Integrate for input range of 0 to 1 for each curve
    AA = zeros(1,size(activation,2),size(activation,3));
    for j = 1:size(activation,3)
        for i = 1:size(activation,2)
            %totalresponse(1,i,j) = trapz(normalized(:,i,j))*0.01;
            %unit spacing = 0.01 (so that maximum value = 1)
            %Activity Area(AA) = sum of the differences between Ai for i
            %concentration and A0 at each point in dose-response curve.
            %Summed value is then divided by 100 for a maximum value of 1.
            %AA=1 indicates maximum response, AA=0 indicates no response. 
            %A0=0 since curves are already normalized to 0.
            AA(1,i,j) = sum(normalized(:,i,j))*0.01;
        end
    end
    fig2 = figure('Position',[100 100 400 300]);
    ax = scatter(AA(:,:,1),AA(:,:,4),"filled");
    grid on
    ax.MarkerFaceAlpha = 0.5;
    xline(0); yline(0);
    title(strcat(choice,": Tension = 0.4"));
    plotfile = strcat(choice,'_final_scatter');
    plotfilename = strcat(homedir,plotdir,choice,'\',plotfile,'.fig');
    saveas(fig2,plotfilename);
    
    ax = scatter(AA(:,:,1),AA(:,:,end),"filled");
    grid on
    ax.MarkerFaceAlpha = 0.5;
    xline(0); yline(0);
    title(strcat(choice,": Tension = 0.9"));
    plotfile = strcat(choice,'_final_scatter');
    plotfilename = strcat(homedir,plotdir,choice,'\',plotfile,'.fig');
    saveas(fig2,plotfilename);
    
    %% Sensitivity Comparisons
    dAA = zeros(size(AA));
    revAA = zeros(size(AA));
    % revAA = zeros(1,length(speciesNames));
    threshold = 0.05;
    for j = 1:size(AA,3)
        dAA(:,:,j) = abs(AA(:,:,j)) - abs(AA(:,:,1));    % calc change in AA
        % dAAhi = abs(AA(:,:,end)) - abs(AA(:,:,1));  % change in AA: high stretch
        % dAA = cat(3,dAAlo,dAAhi);
        
        % flag for reversals - flag=0 unless:
        % if AA(tension)<-0.05 && AA(basal)>0.05, flag = -1 
        % if AA(tension)>0.05 && AA(basal)<-0.05, flag = 1
        for i = 1:size(AA,2)
            if AA(1,i,1)>threshold && AA(1,i,j)<(-1*threshold)
                revAA(1,i,j) = -1;
            elseif AA(1,i,1)<(-1*threshold) && AA(1,i,j)>threshold
                revAA(1,i,j) = 1;
            end
        end
    end
    
    %% Data Export
    filename_AA = strcat(homedir,datadir,'dAA_',choice,'_final.mat');
    save(filename_AA,'dAA')
    filename_rev = strcat(homedir,datadir,'rev_AA_',choice,'_final.mat');
    save(filename_rev,'revAA')
end