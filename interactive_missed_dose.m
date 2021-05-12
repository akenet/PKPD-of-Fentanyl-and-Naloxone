% Systems Pharmacology and Personalized Medicine
% Adam Kenet, Shiker Nair, Lydia Fozo, Amy van Ee
% May 2021
% Final Project

%%%%%%%%%%%%%%   THREE COMPARTMENT MODEL   %%%%%%%%%%%%%%%
%     Compartment 1 -- Blood                             % 
%     Compartment 2 -- Body (scarcely perfused, deep)    % 
%     Compartment 3 -- Brain (highly perfused, shallow)  % 
%     Compartment 4 -- Clearance (imaginary)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Needs files:
%        fentanyl_naloxone_missed.m



function [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = interactive_missed_dose(inpt_missed_dose, inpt_late)
% INTERACTIVE_MISSED_DOSE [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = interactive_missed_dose(inpt_missed_dose, inpt_late)
%   Function for interative missed dose scenario.
%   Called by R code with user inputs as parameters.

%   INPUT:
%       inpt_missed_dose  -- [unitless] -- (user input) which fentanyl dose was missed
%       inpt_late         -- [min]      -- (user input) how many minutes late was the fentanyl dose taken

%   OUTPUT:
%       T        -- [min]  -- times corresponding to the concentrations in Y
%       Y(1:9)   -- [mg/L] -- values output from the ODEs in concentrations
%       Y(10:11) -- [mg]   -- values output from the ODEs in amounts

%       Balance_f -- [unitless] -- mass balance of fentanyl each row is a time point
%       Balance_n -- [unitless] -- mass balance of naloxone each row is a time point

%       AUC_f -- [unitless] -- area under the time versus fentanyl concentration curve in the blood
%       AUC_n -- [unitless] -- area under the time versus naloxone concentration curve in the blood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% References for parameter values
% Fentanyl Drug Bank            https://go.drugbank.com/drugs/DB00813
% Ryan et al. 2018              https://www.futuremedicine.com/doi/full/10.2217/pmt-2017-0060
% Cascone et al. 2013           https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3829787/pdf/tm7_p18.pdf
% Scheibe et al. 1984           https://www.jbc.org/article/S0021-9258(18)90693-9/pdf
% Papathanasiou et al. 2019     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6676012/
% Yassen et al. 2007            https://link.springer.com/content/pdf/10.2165/00003088-200746110-00004.pdf
% Dowling et al. 2008           https://pubmed.ncbi.nlm.nih.gov/18641540/
% Moss et al. 2020              https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7297366/pdf/pone.0234683.pdf
% Livingston et al. 2018        https://elifesciences.org/articles/32499
% Pedersen et al. 2019          https://www.sciencedirect.com/science/article/pii/S0028390819301108
% Encinas et al. 2013           https://link.springer.com/content/pdf/10.1007/s40272-013-0029-1.pdf

%% PARAMETERS
p = struct();

p.weight = 70; % average weight of a male [kg]

p.plasma_f = 0.2; % fraction of fentanyl not bound to plasma protein (80% bound) [unitless]

original_abs = 0.5;
p.abs_n = original_abs; % bioavalability of naloxone (percent of naloxone dose absorbed) [unitless] % from Ryan et al. 2018

p.Vd_brain = 1.260 * 0.01; % subvolume of the brain that both fentanyl and naloxone interact with [L] % estimated value

p.Vd_f1 = 13;         % volume of distribution of fentanyl in blood [L] % from Encinas et al. 2013
p.Vd_f2 = 295;        % volume of distribution of fentanyl in body  [L] % from Encinas et al. 2013
p.Vd_f3 = p.Vd_brain; % volume of distribution of fentanyl in brain [L] % estimated value

p.Vd_n1 = 0.408;      % volume of distribution of naloxone in blood [L] % from Papathanasiou et al. 2019
p.Vd_n2 = 1.637;      % volume of distribution of naloxone in body  [L] % from Papathanasiou et al. 2019
p.Vd_n3 = p.Vd_brain; % volume of distribution of naloxone in brain [L] % estimated value

p.k12_f = 0.373;                               % transport rate of fentanyl from blood to body [1/min] % from Cascone et al. 2013
p.k12_n = 29.8/(60*p.Vd_n1)*(p.Vd_n2/p.Vd_n1); % transport rate of naloxone from blood to body [1/min] % from Rowlings et al. 2008 % converted to [1/min] from [L/min] and included volume correction 

p.k21_f = 0.103;             % transport rate of fentanyl from body to blood [1/min] % from Cascone et al. 2013
p.k21_n = 29.8/(60*p.Vd_n2); % transport rate of naloxone from body to blood [1/min] % from Rowlings et al. 2008 % converted to [1/min] from [L/min]

p.k13_f = 0.0367;      % transport rate of fentanyl from blood to brain [1/min] % from Cascone et al. 2013
p.k13_n = 0.8*p.k13_f; % transport rate of naloxone from blood to brain [L/min] % estimated value

p.k31_f = 0.0124;      % transport rate of fentanyl from brain to blood [1/min] % from Cascone et al. 2013
p.k31_n = 1.2*p.k31_f; % transport rate of naloxone from brain to blood [1/min] % estimated value

p.kcl_f = 42/(60*p.Vd_f1);  % clearance rate of fentanyl from blood [L/min] % from Fentanyl Drug Bank   % changed to 1/min
p.kcl_n = 91/(60*p.Vd_n1);  % clearance rate of naloxone from blood [L/min] % from Dowlings et al. 2008 % changed to 1/min

p.kon_f = 0.026; % binding rate of fentanyl to mu opioid receptor [1/(min*nM)] % from Moss et al. 2020 
p.kon_n = 1.073; % binding rate of naloxone to mu opioid receptor [1/(min*nM)] % from Moss et al. 2020

p.koff_f = 0.035;  % dissociation rate of fentanyl from mu opioid receptor[1/min] % from Livingston et al. 2018
p.koff_n = 1.196 ; % dissociation rate of naloxone from mu opioid receptor[1/min] % from Pedersen et al. 2019

original_mOR = 0.1;
p.amt_mOR = original_mOR; % amount of mOR receptors in brain [nmol] % estimated value


%% MISSED DOSE
dose_f = 4;  % [mg] -- dose of fentanyl
num_f  = 6; % [unitless] -- number of doses of fentanyl given
missed_dose = inpt_missed_dose; % (user input) which dose was missed

time_start = 0;
time_end = 180;

interactive_bool = true; % for interactive plot
min_late = inpt_late; % (user input) how many minutes late was the drug taken

% simulation
[T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_missed(p, dose_f, num_f, missed_dose, time_start, time_end, interactive_bool, min_late);

% save data
save_data('missed_dose_run', T, Y, Balance_f, Balance_n, AUC_f, AUC_n); 


%%% base case
min_late = 0; % taken on time

% simulation
[T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_missed(p, dose_f, num_f, missed_dose, time_start, time_end, interactive_bool, min_late);

% save data
save_data('missed_dose_base', T, Y, Balance_f, Balance_n, AUC_f, AUC_n); 
disp('done')
end


%% Function to save data
function save_data(scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n)

outputdata_T = strcat('int_T_(',scenario,').mat');
outputdata_Y = strcat('int_Y_(',scenario,').mat');
outputdata_Balance_f = strcat('int_Balance_f_(',scenario,').mat');
outputdata_Balance_n = strcat('int_Balance_n_(',scenario,').mat');
outputdata_AUC_f = strcat('int_AUC_f_(',scenario,').mat');
outputdata_AUC_n = strcat('int_AUC_n_(',scenario,').mat');

save(outputdata_T, 'T');
save(outputdata_Y, 'Y');
save(outputdata_Balance_f, 'Balance_f');
save(outputdata_Balance_n, 'Balance_n');
save(outputdata_AUC_f, 'AUC_f');
save(outputdata_AUC_n, 'AUC_n');

end





