% SNM Drug Screen: Combinatory Perturbations
% 04.21.2020 JR
% Script wraps drugscreen() fcn in order to run using HPC. 
% 
% Edited 05.07.2020 JR: Added input levels using idealized curves (Zeigler
% MB 2020) + Added multiple screens

%% define drugscreen arguments:
datadir = "data\5_DrugScreen\";
screentype = "combinations";    % screen type: combinations (string)
tension_all = [0.6];            % tension level (doubles)
choice_all = ["KO"];            % knock-out or over-expression (strings)
peak = 0.6;                     % peak input level for model (for Input_12_19)
[InputCurves,~,inputNode,~,~] = InputCurve_12_19(peak,peak);
t = 2*7*24;                     % time in h
t0_analysis = 168;              % 168 h added as baseline
t_analysis = t+t0_analysis;     % matches InputCurve script
inputs_baseline = InputCurves(:,t0_analysis);
inputs_infarct = InputCurves(:,t_analysis);

inputs_fraction = 1;
inputs_remote = (inputs_fraction*(inputs_infarct-inputs_baseline))+inputs_baseline;
inputs_all = [inputs_remote inputs_infarct];

%% perform drug screens
fprintf("\n<<<<<< Starting Screen: t=%dwks | LVC=%.2f >>>>>>\n",t/7/24,inputs_fraction)
parpool()

for choice = choice_all
    for sim = 1:length(tension_all)
        inputs = inputs_all(:,sim);
        tension = tension_all(sim);
        act_delta = DrugScreenSimulations(screentype,inputs,inputNode,tension,choice);
        % save results
        tensionname = replace(string(tension),".","");
        filename = strcat(datadir,screentype,"_tension",tensionname,"_",choice,".mat");
        save(filename,'act_delta');
    end
end

delete(gcp)
