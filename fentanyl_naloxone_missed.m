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
%        fentanyl_naloxone_main.m


function [T_f, Y_f, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_missed(p, dose_f, num_f, missed_dose, time_start, time_end, interactive_bool, min_late, q)
% FENTANYL_NALOXONE_MISSED [T_f, Y_f, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_missed(p, dose_f, num_f, missed_dose, time_start, time_end, interactive_bool, min_late, q)
%   Function for adminstration of fentanyl and then naloxone (for survival scenarios).
%   Simulates one fentanyl administration followed by repeated naloxone administration.
%   NOTE: modification of fentanyl_naloxone_fxnx

%   INPUT:
%       p -- stucture with parameter values

%       dose_f -- [mg]       -- amount of one dose of fentanyl
%       time_f -- [min]      -- time between repeated fentanyl administrations

%       missed_dose  -- [unitless] -- which fentanyl dose was missed

%       time_start -- [min]      -- starting time of simulation (default=0)
%       time_end   -- [min]      -- ending time of simulation 

%       interactive_bool -- [boolean] -- true when function is called by interactive visualization code
%       min_late         -- [min]     -- how many minutes late was the fentanyl dose take (only for interactive visualization)

%       q -- [unitless] -- takes values of 1, 2, 3, 4 for missed dose scenario in driver (acts as 'min_late' for non-interactive code)


%   OUTPUT:
%       T_f        -- [min]  -- times corresponding to the concentrations in Y
%       Y_f(1:9)   -- [mg/L] -- values output from the ODEs in concentrations
%       Y_f(10:11) -- [mg]   -- values output from the ODEs in amounts

%       Balance_f -- [unitless] -- mass balance of fentanyl each row is a time point
%       Balance_n -- [unitless] -- mass balance of naloxone each row is a time point

%       AUC_f -- [unitless] -- area under the time versus fentanyl concentration curve in the blood
%       AUC_n -- [unitless] -- area under the time versus naloxone concentration curve in the blood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


time_array = linspace(time_start, time_end, num_f+1);


if (nargin == 7) || (~interactive_bool)   % default case w/o interactive bool
    time_array(missed_dose) = time_array(missed_dose) + (time_array(2)-time_array(1))*(q/5);
end
if (nargin == 8) && (interactive_bool) % for interactive plot
    time_array(missed_dose) = time_array(missed_dose) + min_late;
end
disp(time_array)

T_f=[];
Y_f=[];

p.dose_f = dose_f; 
y0 = zeros(11,1);
y0(1) = (p.dose_f/p.Vd_f1);
y0(7) = p.amt_mOR/p.Vd_brain;
p.dose_n = 0; 
p.i_n = 0;

for i = 1:(length(time_array)-1) 
        
        % set dose number
        p.i_f = i; % keep track of fentanyl administrations
        
        % run simulation
        [Ti,Yi,Balance_f,Balance_n,AUC_f,AUC_n] = fentanyl_naloxone_main([time_array(i), time_array(i+1)], y0, p);
        
        % reset initial conditions
        y0 = Yi(end,:);
        y0(1) = y0(1) + (p.dose_f/p.Vd_f1);   % add the next fentanyl dose as a concentration

        
        % combine values for fentanyl administration
        Y_f = vertcat(Y_f(1:end-1,:),Yi);
        T_f = vertcat(T_f(1:end-1,:),Ti);
        
end
end