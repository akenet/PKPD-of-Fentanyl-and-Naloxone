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
%        fentanyl_naloxone_eqns.m


function [T,Y,Balance_f,Balance_n,AUC_f,AUC_n] = fentanyl_naloxone_main(tspan,y0,p)
% FENTANYL_NALOXONE_MAIN [T,Y,Balance,AUC] = fentanyl_naloxone_main(tspan,y0,p)
%   Solves the system of ODEs given in the fentanyl_naloxone_eqns file based
%   on the parameters given. Inputs are the length of time, initial
%   conditions, and parameters; outputs are the time, concentrations, mass
%   balance, and area under the curve from the ODE solver.
%
%   INPUT:
%       tspan -- [min] vector of times to use for the beginning and end of the system
%       y0    -- [mg or mM] -- initial conditions when starting simulation (default is zeros except for inital fentanyl and free mOR receptors)
%       p     -- stucture with parameter values

%   OUTPUT:
%       T        -- [min]  -- times corresponding to the concentrations in Y
%       Y(1:9)   -- [mg/L] -- values output from the ODEs in concentrations
%       Y(10:11) -- [mg]   -- values output from the ODEs in amounts

%       Balance_f -- [unitless] -- mass balance of fentanyl each row is a time point
%       Balance_n -- [unitless] -- mass balance of naloxone each row is a time point

%       AUC_f -- [unitless] -- area under the time versus fentanyl concentration curve in the blood
%       AUC_n -- [unitless] -- area under the time versus naloxone concentration curve in the blood

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Solver

% Options for ODE Solver
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,...
    'InitialStep', 1e-2);

% Run ODE Solver
[T,Y] = ode15s(@fentanyl_naloxone_eqns, tspan, y0, options, p);


%% Mass Balance of Fentanyl
TotalD_f(:,1) = Y(:,1)*p.Vd_f1; % amount of free  fentanyl in blood [mg]
TotalD_f(:,2) = Y(:,2)*p.Vd_f2; % amount of free  fentanyl in body  [mg]
TotalD_f(:,3) = Y(:,3)*p.Vd_f3; % amount of free  fentanyl in brain [mg]
TotalD_f(:,4) = Y(:,8)*p.Vd_f3; % amount of bound fentanyl in brain [mg]

DrugIn_f = p.dose_f * p.i_f; % total amount of fentanyl going into system [mg]
DrugOut_f = Y(:,10);         % total amount of fentanyl cleared from the system [mg]

Balance_f = DrugIn_f - DrugOut_f - TotalD_f(:,1) - TotalD_f(:,2) - TotalD_f(:,3) - TotalD_f(:,4);

if any(abs(Balance_f) > 1e-6)
    fprintf('Mass of fentanyl is not balanced: %g\n', Balance_f)
end

%% Mass Balance of Naloxone
TotalD_n(:,1) = Y(:,4)*p.Vd_n1; % amount of free  naloxone in blood [mg]
TotalD_n(:,2) = Y(:,5)*p.Vd_n2; % amount of free  naloxone in body  [mg]
TotalD_n(:,3) = Y(:,6)*p.Vd_n3; % amount of free  naloxone in brain [mg]
TotalD_n(:,4) = Y(:,9)*p.Vd_n3; % amount of bound naloxone in brain [mg]

DrugIn_n = p.dose_n * p.i_n * p.abs_n; % total amount of naloxone going into system [mg]

DrugOut_n = Y(:,11);                   % total amount of naloxone cleared from the system [mg]

Balance_n = DrugIn_n - DrugOut_n - TotalD_n(:,1) - TotalD_n(:,2) - TotalD_n(:,3) - TotalD_n(:,4);

if any(abs(Balance_n) > 1e-6)
    fprintf('Mass of naloxone is not balanced: %g\n', Balance_n)
end

%% AUC
AUC_f = trapz(T,Y(:,3)); % AUC of free fentanyl in brain
AUC_n = trapz(T,Y(:,6)); % AUC of free naloxone in brain
end