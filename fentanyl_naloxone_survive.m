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


function [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end, y0)
% FENTANYL_NALOXONE_SURVIVE [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end, y0)
%   Function for adminstration of fentanyl and then naloxone (for survival scenarios).
%   Simulates one fentanyl administration followed by repeated naloxone administration.
%   NOTE: modification of fentanyl_naloxone_fxnx

%   INPUT:
%       p -- stucture with parameter values

%       dose_f -- [mg]       -- amount of one dose of fentanyl
%       time_f -- [min]      -- time between repeated fentanyl administrations
%       num_f  -- [unitless] -- number of doses of fentanyl given

%       delay  -- [min]      -- time between last fentanyl administration and first naloxone administration

%       dose_n -- [mg]       -- amount of one dose of naloxone
%       time_n -- [min]      -- time between repeated naloxone administrations
%       num_n  -- [unitless] -- number of doses of naloxone given

%       time_start -- [min]      -- starting time of simulation (default=0)
%       time_end   -- [min]      -- ending time of simulation 
%       y0         -- [mg or mM] -- initial conditions when starting simulation (default is zeros except for inital fentanyl and free mOR receptors)


%   OUTPUT:
%       T        -- [min]  -- times corresponding to the concentrations in Y
%       Y(1:9)   -- [mg/L] -- values output from the ODEs in concentrations
%       Y(10:11) -- [mg]   -- values output from the ODEs in amounts

%       Balance_f -- [unitless] -- mass balance of fentanyl each row is a time point
%       Balance_n -- [unitless] -- mass balance of naloxone each row is a time point

%       AUC_f -- [unitless] -- area under the time versus fentanyl concentration curve in the blood
%       AUC_n -- [unitless] -- area under the time versus naloxone concentration curve in the blood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 10   % default case w/o y0
    time_start = 0;
    
    % initial conditions
    p.dose_f = dose_f;  % dose of fentanyl [mg]
    p.dose_n = 0; % no naloxone when administering fentanyl
    p.i_n = 0; % no naloxone when administering fentanyl
    
    y0 = zeros(11,1);
    y0(1) = p.dose_f/p.Vd_f1;      % add the first fentanyl dose as a concentration
    y0(7) = p.amt_mOR/p.Vd_brain;  % concentration of free mu opiod receptor in brain [nM]
    
end

%% Fentanyl

% initial conditions
p.dose_f = dose_f;  %dose of fentanyl [mg]
p.dose_n = 0; % no naloxone when administering fentanyl
p.i_n = 0; % no naloxone when administering fentanyl


y0 = zeros(11,1);
y0(1) = p.dose_f/p.Vd_f1;      % add the first fentanyl dose as a concentration
y0(7) = p.amt_mOR/p.Vd_brain;  % concentration of free mu opiod receptor in brain [nM]


% initialize arrays
T_f=[];
Y_f=[];


if num_f == 1
    % set dose number
    p.i_f = 1; % only one fentanyl administrations
    tspan_f = [time_start, time_start+delay];
    % run simulation
    [Ti,Yi,Balance_f,Balance_n,AUC_f,AUC_n] = fentanyl_naloxone_main(tspan_f,y0,p);
    
    % combine values for fentanyl administration
    Y_f = Yi;
    T_f = Ti;
    
elseif num_f > 1
    % repeated fentanyl doses
    for i = 1:num_f-1 % administer all doses except for last dose
        
        % set dose number
        p.i_f = i; % keep track of fentanyl administrations
        
        % run simulation
        [Ti,Yi,Balance_f,Balance_n,AUC_f,AUC_n] = fentanyl_naloxone_main(tspan_f,y0,p);
        
        % reset initial conditions
        y0 = Yi(end,:);
        y0(1) = y0(1) + (p.dose_f/p.Vd_f1);    % add the next fentanyl dose as a concentration
        tspan_f = [Ti(end), Ti(end) + time_f]; % shift the time range for next administration
        
        % combine values for fentanyl administration
        Y_f = vertcat(Y_f(1:end-1,:),Yi);
        T_f = vertcat(T_f(1:end-1,:),Ti);
        
    end

else
    disp("Error: Invalid number of fentanyl doses")
end


%% Naloxone

% initial conditions
p.dose_n = dose_n;                 % dose of naloxone NOT adjusted by bioavalability [mg]
p.dose_n_adj = p.dose_n * p.abs_n; % dose of naloxone     adjusted by bioavalability [mg]

y0 = Y_f(end,:); % set initial conditions at time of first naloxone administration to values "delay" minutes after last fentanyl administration
y0(4) = y0(4) + (p.dose_n_adj/p.Vd_n1); % add the first naloxone dose as a concentration

% t_array = linspace(T_f(end), time_end, num_n+1);
t_array = [linspace(T_f(end), T_f(end)+(num_n*time_n), num_n+1)];
t_array(num_n+1) = time_end; 



% initialize arrays
T_n=[];
Y_n=[];

% repeated naloxone doses
if num_n > 0
    for i = 1:length(t_array)-1 % number of doses of naloxone to administer
        
        % set dose number
        p.i_f = num_f; % fentanyl has been administered this many times already
        p.i_n = i;     % keep track of naloxone administrations
        
        % run simulation
        [Ti,Yi,Balance_f,Balance_n,AUC_f,AUC_n] = fentanyl_naloxone_main([t_array(i), t_array(i+1)],y0,p);
        
        % reset initial conditions
        y0 = Yi(end,:);
        y0(4) = y0(4) + (p.dose_n_adj/p.Vd_n1); % add the next naloxone dose as a concentration
        
        % combine values for naloxone administration
        Y_n = vertcat(Y_n(1:end-1,:),Yi);
        T_n = vertcat(T_n(1:end-1,:),Ti);
        
    end
end

% combine values for entire simulation
Y = vertcat(Y_f(1:end-1,:),Y_n);
T = vertcat(T_f(1:end-1,:),T_n);
end