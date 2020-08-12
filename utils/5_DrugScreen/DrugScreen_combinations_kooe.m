% SNM Drug Screen: Combinatory Perturbations
% 04.21.2020 JR
% Script wraps drugscreen() fcn in order to run using HPC. 
% 
% Edited 05.07.2020 JR: Added input levels using idealized curves (Zeigler
% MB 2020) + Added multiple screens

%% define drugscreen arguments:
datadir = "data\5_DrugScreen\";
screentype = "combinations";        % screen type: combinations (string)
tension_all = [0.1 0.6];                % tension level (doubles)
choice_all = ["KO-OE"];             % knock-out or over-expression (strings)
combinations = nchoosek(1:109,2);   % to break down into smaller simulations
c_len = round(length(combinations)/4);
c_1 = combinations(1:c_len,:);
c_2 = combinations(c_len+1:2*c_len,:);
c_3 = combinations(2*c_len+1:3*c_len,:);
c_4 = combinations(3*c_len+1:3.5*c_len,:);
c_5 = combinations(3.5*c_len+1:end,:);
peak = 0.6;                         % peak input level for model (for Input_12_19)
[InputCurves,~,inputNode,~,~] = InputCurve_12_19(peak,peak);
t = 2*7*24;                         % time in h
t0_analysis = 168;                  % 168 h added as baseline
t_analysis = t+t0_analysis;         % matches InputCurve script
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
        act_delta_1 = DrugScreenSimulation(screentype,inputs,inputNode,tension,choice,c_1);
        act_delta_2 = DrugScreenSimulation(screentype,inputs,inputNode,tension,choice,c_2);
        act_delta_3 = DrugScreenSimulation(screentype,inputs,inputNode,tension,choice,c_3);
        act_delta_4 = DrugScreenSimulation(screentype,inputs,inputNode,tension,choice,c_4);
        act_delta_5 = DrugScreenSimulation(screentype,inputs,inputNode,tension,choice,c_5);
        act_delta = [act_delta_1;act_delta_2;act_delta_3;act_delta_4;act_delta_5];
        % save results
        tensionname = replace(string(tension),".","");
        filename = strcat(datadir,screentype,"_tension",tensionname,"_",choice,".mat");
        save(filename, 'act_delta');
    end
end

delete(gcp)
