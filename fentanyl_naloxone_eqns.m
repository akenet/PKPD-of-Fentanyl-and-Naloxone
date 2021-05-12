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

% FENTANYL_NALOXONE_EQNS
%   Differential equations for three compartment model of fentanyl and naloxone


function dydt = fentanyl_naloxone_eqns(~,y,p)
%% Fentanyl
% y(1 )-- fent_bl -- concentration of free fentanyl in blood [nM]
% y(2) -- fent_bo -- concentration of free fentanyl in body  [nM]
% y(3) -- fent_br -- concentration of free fentanyl in brain [nM]

%% Naloxone
% y(4) -- nalox_bl -- concentration of free naloxone in blood [nM]
% y(5) -- nalox_bo -- concentration of free naloxone in body  [nM]
% y(6) -- nalox_br -- concentration of free naloxone in brain [nM]

%% mu Opiod Receptor
% y(7) -- muR_free -- concentration of free mu opiod receptor in brain [nM]
% y(8) -- muR_f -- concentration of mu opiod receptor bound to fentanyl in brain [nM]
% y(9) -- muR_n -- concentration of mu opiod receptor bound to naloxone in brain [nM]

%% Clearance
% y(10) -- cl_f -- amount of cleared fentanyl [mg]
% y(11) -- cl_n -- amount of cleared naloxone [mg]


%% Parameters (see driver for units and values)
% p.plasma_f -- percent of fentanyl not bound to plasma protein
% p.dose_f   -- dose of fentanyl
% p.dose_n   -- dose of naloxone
% p.abs_n    -- percent of naloxone dose absorbed

% p.i_f -- number of doses of fentanyl
% p.i_n -- number of doses of naloxone

% p.k12_f -- transport rate of fentanyl from blood to body 
% p.k12_n -- transport rate of naloxone from blood to body

% p.k21_f -- transport rate of fentanyl from body to blood 
% p.k21_n -- transport rate of naloxone from body to blood 

% p.k13_f -- transport rate of fentanyl from blood to brain 
% p.k13_n -- transport rate of naloxone from blood to brain 

% p.k31_f -- transport rate of fentanyl from brain to blood 
% p.k31_n -- transport rate of naloxone from brain to blood 

% p.kcl_f -- clearance rate of fentanyl from blood 
% p.kcl_n -- clearance rate of naloxone from blood 

% p.kon_f -- binding rate of fentanyl to mu opioid receptor 
% p.kon_n -- binding rate of naloxone to mu opioid receptor 

% p.koff_f -- dissociation rate of fentanyl from mu opioid receptor 
% p.koff_n -- dissociation rate of naloxone from mu opioid receptor 

% p.Vd_f1 -- volume of distribution of fentanyl of blood
% p.Vd_f2 -- volume of distribution of fentanyl of body
% p.Vd_f3 -- volume of distribution of fentanyl of brain

% p.Vd_n1 -- volume of distribution of naloxone of blood
% p.Vd_n2 -- volume of distribution of naloxone of body
% p.Vd_n3 -- volume of distribution of naloxone of brain


%% Equations
dydt = zeros(11,1);

% concentration of free fentanyl in blood [nM] % only account for dose when it is continuous
dydt(1) = - p.k12_f*y(1)*p.plasma_f                   - p.k13_f*y(1)*p.plasma_f ...
          + p.k21_f*y(2)*(p.Vd_f2/p.Vd_f1) + p.k31_f*y(3)*(p.Vd_f3/p.Vd_f1) ...
          - p.kcl_f*y(1)*p.plasma_f;
      
% concentration of free fentanyl in body [nM]
dydt(2) =  p.k12_f*y(1)*p.plasma_f*(p.Vd_f1/p.Vd_f2) - p.k21_f*y(2);


% concentration of free fentanyl in brain [nM]
dydt(3) =   p.k13_f*y(1)*p.plasma_f*(p.Vd_f1/p.Vd_f3) - p.k31_f*y(3) ...
          + p.koff_f*y(8) - p.kon_f*y(3)*y(7);
            
% concentration of free naloxone in blood [nM]
dydt(4) = -p.k12_n*y(4)                   - p.k13_n*y(4)  ...
          + p.k21_n*y(5)*(p.Vd_n2/p.Vd_n1) + p.k31_n*y(6)*(p.Vd_n3/p.Vd_n1) ...
          - p.kcl_n*y(4);
      
% concentration of free naloxone in body [nM]
dydt(5) =  p.k12_n*y(4)*(p.Vd_n1/p.Vd_n2) - p.k21_n*y(5);


% concentration of free naloxone in brain [nM]
dydt(6) =   p.k13_n*y(4)*(p.Vd_n1/p.Vd_n3) - p.k31_n*y(6) ...
          + p.koff_n*y(9) - p.kon_n*y(6)*y(7);

% concentration of free mu opiod receptor in brain [nM]
dydt(7) =   p.koff_f*y(8) - p.kon_f*y(3)*y(7) ...
          + p.koff_n*y(9) - p.kon_n*y(6)*y(7);

% concentration of mu opiod receptor bound to fentanyl in brain [nM]
dydt(8) = p.kon_f*y(3)*y(7) - p.koff_f*y(8);

% concentration of mu opiod receptor bound to naloxone in brain [nM]
dydt(9) = p.kon_n*y(6)*y(7) - p.koff_n*y(9);

% amount of cleared fentanyl [mg]
dydt(10) = p.kcl_f*y(1)*p.plasma_f*p.Vd_f1; 

% amount of cleared naloxone [mg]
dydt(11) = p.kcl_n*y(4)*p.Vd_n1;
end
