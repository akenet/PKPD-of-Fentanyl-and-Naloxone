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
%        fentanyl_naloxone_fxnx.m
%        fentanyl_naloxone_missed.m
%        fentanyl_naloxone_survive.m


% FENTANYL_NALOXONE_DRIVER
%  Different scenarios of the model: repeated dosing, missed dosing,
%  survival, population pharmacokinetics, sensitivity analysis


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

NumberOfSubjects = 100; % distribution parameter

%% SCENARIOS

%                         f25n0 -- 25 repeated doses of fentanyl; no naloxone
%                          f5n3 --  5 repeated doses of fentanyl; 3 repeated doses of naloxone
%                   missed_dose -- 4th dose (of 6) is missed and taken at a later time
%               s2_add_naloxone -- survival simulation 1: add naloxone
%    s3_inc_delay_half_naloxone -- survival simulation 2: increase delay, half naloxone
%            s4_double_naloxone -- survival simulation 3: double naloxone
%                     pop_smoke -- population pk: smoking
%                      pop_nose -- population pk: naloxone absorption
%                     pop_smose -- (need to run pop_smoke and pop_nose first) population pk: smoking and naloxone absorption
%                   sensitivity -- sensitivity analysis


root = pwd; % root directory

saveFlag = true;  % Do you want to save the data?


scenarios = ["f25n0","f5n3","missed_dose","s2_add_naloxone","s3_inc_delay_half_naloxone","s4_double_naloxone","pop_smoke","pop_nose","pop_smose","sensitivity"]; % which simulation(s) to run

for run = 1:length(scenarios)
    scenario = scenarios(run)
    
    switch scenario

        case "f25n0"
            %%%% 1 mg of fentanyl administered every 45 minutes for a total of 25 doses
            %%%% 0 minutes between last fentanyl administration and first naloxone administration
            %%%% 0 mg of naloxone administered every 0 minutes for a total of 0 doses
            
            dose_f = 1;  % [mg]       -- amount of one dose of fentanyl
            time_f = 45; % [min]      -- time between repeated fentanyl administrations
            num_f  = 25; % [unitless] -- number of doses of fentanyl given
            
            delay  = 0;  % [min]      -- time between last fentanyl administration and first naloxone administration
            
            dose_n = 0;  % [mg]       -- amount of one dose of naloxone
            time_n = 0;  % [min]      -- time between repeated naloxone administrations
            num_n  = 0;  % [unitless] -- number of doses of naloxone given
            
            [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_fxnx(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n);
            
            if(saveFlag); save_data(root, scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end
            
        case "f5n3"
            %%%% 1 mg of fentanyl administered every 45 minutes for a total of 5 doses
            %%%% 7 minutes between last fentanyl administration and first naloxone administration
            %%%% 4 mg of naloxone administered every 3 minutes for a total of 3 doses
            
            dose_f = 1;  % [mg]       -- amount of one dose of fentanyl
            time_f = 45; % [min]      -- time between repeated fentanyl administrations
            num_f  = 5;  % [unitless] -- number of doses of fentanyl given
            
            delay  = 7;  % [min]      -- time between last fentanyl administration and first naloxone administration
            
            dose_n = 4;  % [mg]       -- amount of one dose of naloxone
            time_n = 3;  % [min]      -- time between repeated naloxone administrations
            num_n  = 3;  % [unitless] -- number of doses of naloxone given
            
            [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_fxnx(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n);
            
            if(saveFlag); save_data(root, scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end
            
        case "missed_dose"
            % 4th dose (of 6) is missed and taken at a later time
            
            %%%% i=0 (base) -- 4 mg of fentanyl administered every 30 minutes for a total of 6 doses
            %%%% i=1        -- 4th dose taken 6  minutes late
            %%%% i=2        -- 4th dose taken 12 minutes late
            %%%% i=3        -- 4th dose taken 18 minutes late
            %%%% i=4        -- 4th dose taken 24 minutes late

            dose_f = 4;       % [mg]       -- dose of fentanyl
            num_f  = 6;       % [unitless] -- number of doses of fentanyl given
            missed_dose = 4;  % which dose was missed
            
            % length of simulation
            time_start = 0;
            time_end = 180;
            
            for i=0:4 % loop through scenarios
                
                [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_missed(p, dose_f, num_f, missed_dose, time_start, time_end, false, 0, i);
                
                % change file name
                if(i==0); scenario = "missed_dose_base"; end
                if(i==1); scenario = "missed_dose_1"; end
                if(i==2); scenario = "missed_dose_2"; end
                if(i==3); scenario = "missed_dose_3"; end
                if(i==4); scenario = "missed_dose_4"; end
                
                if(saveFlag); save_data(root, scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end

            end
            
        case "s2_add_naloxone"
            % Survival Simulation 1
            
            %%%% 1 dose of fentanyl (2 mg) given at t=0
            %%%% 1 minute later, 3 doses of naloxone (4 mg) given 3 min apart
          
            
            % initial concentration of free receptor
            y0_7 = p.amt_mOR/p.Vd_brain; 
            
            % length of simulation
            time_start = 0;
            time_end = 30;
            
            
            % BASE CASE -- NO NALOXONE
            dose_f = 2; % [mg]       -- amount of one dose of fentanyl
            time_f = 0; % [min]      -- time between repeated fentanyl administrations
            num_f  = 1; % [unitless] -- number of doses of fentanyl given
            
            delay  = time_end;  % [min]  -- no naloxone given so run simulation until time_end
            
            % no naloxone given
            dose_n = 0;  % [mg]       -- amount of one dose of naloxone
            time_n = 0;  % [min]      -- time between repeated naloxone administrations
            num_n  = 0;  % [unitless] -- number of doses of naloxone given
            
            [base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
            
            % save base case data
            if(saveFlag); save_data(root, strcat(scenario,'_base'), base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end
            

            % SIMULATION -- ADD NALOXONE
            delay  = 1;  % [min]      -- time between last fentanyl administration and first naloxone administration
            dose_n = 4;  % [mg]       -- amount of one dose of naloxone
            num_n  = 3;  % [unitless] -- number of doses of naloxone given
            time_n = 3;  % [min]      -- time between repeated naloxone administrations
            
            [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
            
            % save survive data
            if(saveFlag); save_data(root, scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end

        case "s3_inc_delay_half_naloxone"
            % Survival Simulation 2
            
            %%%% 1 dose of fentanyl (2 mg) given at t=0
            %%%% 10 minutes later, 3 doses of naloxone (2 mg) given 3 min apart
          
            
            % initial concentration of free receptor
            y0_7 = p.amt_mOR/p.Vd_brain; 
            
            % length of simulation
            time_start = 0;
            time_end = 30;
            
            
            % BASE CASE -- NO NALOXONE
            dose_f = 2; % [mg]       -- amount of one dose of fentanyl
            time_f = 0; % [min]      -- time between repeated fentanyl administrations
            num_f  = 1; % [unitless] -- number of doses of fentanyl given
            
            delay  = time_end;  % [min]  -- no naloxone given so run simulation until time_end
            
            % no naloxone given
            dose_n = 0;  % [mg]       -- amount of one dose of naloxone
            time_n = 0;  % [min]      -- time between repeated naloxone administrations
            num_n  = 0;  % [unitless] -- number of doses of naloxone given
            
            [base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
            
            % save base case data
            if(saveFlag); save_data(root, strcat(scenario,'_base'), base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end
            

            % SIMULATION -- ADD NALOXONE
            delay  = 10; % [min]      -- time between last fentanyl administration and first naloxone administration % changed from sim 1
            dose_n = 2;  % [mg]       -- amount of one dose of naloxone  % changed from sim 1
            num_n  = 3;  % [unitless] -- number of doses of naloxone given
            time_n = 3;  % [min]      -- time between repeated naloxone administrations
            
            [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
            
            % save survive data
            if(saveFlag); save_data(root, scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end

        case "s4_double_naloxone"
            % Survival Simulation 3
            
            %%%% 1 dose of fentanyl (2 mg) given at t=0
            %%%% 10 minutes later, 3 doses of naloxone (4 mg) given 3 min apart
          
            
            % initial concentration of free receptor
            y0_7 = p.amt_mOR/p.Vd_brain; 
            
            % length of simulation
            time_start = 0;
            time_end = 30;
            
            
            % BASE CASE -- NO NALOXONE
            dose_f = 2; % [mg]       -- amount of one dose of fentanyl
            time_f = 0; % [min]      -- time between repeated fentanyl administrations
            num_f  = 1; % [unitless] -- number of doses of fentanyl given
            
            delay  = time_end;  % [min]  -- no naloxone given so run simulation until time_end
            
            % no naloxone given
            dose_n = 0;  % [mg]       -- amount of one dose of naloxone
            time_n = 0;  % [min]      -- time between repeated naloxone administrations
            num_n  = 0;  % [unitless] -- number of doses of naloxone given
            
            [base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
            
            % save base case data
            if(saveFlag); save_data(root, strcat(scenario,'_base'), base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end
            

            % SIMULATION -- ADD NALOXONE
            delay  = 10; % [min]      -- time between last fentanyl administration and first naloxone administration
            dose_n = 4;  % [mg]       -- amount of one dose of naloxone  % changed from sim 2
            num_n  = 3;  % [unitless] -- number of doses of naloxone given
            time_n = 3;  % [min]      -- time between repeated naloxone administrations
            
            [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
            
            % save survive data
            if(saveFlag); save_data(root, scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n); end

        case "pop_smoke"
            % Population PK -- Smoking
           
            
            % Change concentration of receptors distribution
            % heavier smoker -- range from 8 - 14% reduction in mOR availability
            
            %%%%%%%%% create smoking distributions %%%%%%%%%
            
            % randomly generate percentage of mOR availability reduced by between 0 and 14%
            % 14 of every 100 American adults are smokers
            
            rng('default') % reproducibility
            
            % max and min percentages of total mOR availability
            min = 0.86;
            max = 1;
            
            % get amount of mOR from a uniform distribution
            smokers_mOR = ((max - min).*rand(NumberOfSubjects, 1) + min) .* original_mOR;
            
            
            
            %%%%%%%%% run simulation (same as s3_inc_delay_half_naloxone) %%%%%%%%%
            
            %%%% 1 dose of fentanyl (2 mg) given at t=0
            %%%% 10 minutes later, 3 doses of naloxone (2 mg) given 3 min apart
            
            
            % initialize array to store AUCs
            popsmoke_AUC = [];
            popsmoke_surv = [];
            
            % for loop to get survival AUC for each patient
            for i = 1:NumberOfSubjects
                
                % set new initial mOR amount for input
                p.amt_mOR = smokers_mOR(i);
                
                % initial concentration of free receptor
                y0_7 = p.amt_mOR/p.Vd_brain;
                
                % length of simulation
                time_start = 0;
                time_end = 30;
                
                
                % BASE CASE -- NO NALOXONE
                dose_f = 2; % [mg]       -- amount of one dose of fentanyl
                time_f = 0; % [min]      -- time between repeated fentanyl administrations
                num_f  = 1; % [unitless] -- number of doses of fentanyl given
                
                delay  = time_end;  % [min]  -- no naloxone given so run simulation until time_end
                
                % no naloxone given
                dose_n = 0;  % [mg]       -- amount of one dose of naloxone
                time_n = 0;  % [min]      -- time between repeated naloxone administrations
                num_n  = 0;  % [unitless] -- number of doses of naloxone given
                
                [base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
                
                % receptor occupancy
                base_mu_occ_f = Y(:,8) ./ y0_7;
                
                dead_base_auc = trapz(base_T,base_mu_occ_f);
                
                
                % SIMULATION -- ADD NALOXONE
                delay  = 10; % [min]      -- time between last fentanyl administration and first naloxone administration
                dose_n = 2;  % [mg]       -- amount of one dose of naloxone
                num_n  = 3;  % [unitless] -- number of doses of naloxone given
                time_n = 3;  % [min]      -- time between repeated naloxone administrations
                
                [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
                
                % receptor occupancy
                mu_occ_f = Y(:,8) ./ y0_7;
                dead_test_auc = trapz(T,mu_occ_f);
                
                % store data
                popsmoke_AUC(i) = dead_test_auc; % AUC
                popsmoke_surv(i) = 1 - (dead_test_auc/dead_base_auc); % survival
                
            end
            
            popsmoke_AUC = transpose(popsmoke_AUC);
            popsmoke_surv = transpose(popsmoke_surv);
            
            %%% save data for plotting in R
            if(saveFlag)
                new_dir = strcat(root,'/Data/',scenario);
                mkdir(new_dir)
                cd(new_dir)
                
                save("smokers_mOR.mat", 'smokers_mOR');
                save("popsmoke_AUC.mat", 'popsmoke_AUC');
                save("popsmoke_surv.mat", 'popsmoke_surv');
                
                cd(root)
            end

        case "pop_nose"
            % Population PK -- Naloxone Absorption
            

            %%%%%%%%% create nose absorption distributions %%%%%%%%%
            
            rng('default') % reproducibility
            
            % percentages
            mean = 1;      
            SD = 0.1;

            % randomly generate weights from Normal distribution
            noses = normrnd(mean, SD, [1, NumberOfSubjects]);


            %%%%%%%%% run simulation (same as s3_inc_delay_half_naloxone) %%%%%%%%%
            
            %%%% 1 dose of fentanyl (2 mg) given at t=0
            %%%% 10 minutes later, 3 doses of naloxone (2 mg) given 3 min apart
            
            
            % initialize array to store AUCs
            popnose_AUC = [];
            popnose_surv = [];

            % for loop to get survival AUC for each patient
            for i = 1:NumberOfSubjects

                % set new nose abs for input
                p.abs_n = noses(i) * original_abs;
                % reset just in case
                p.amt_mOR = original_mOR;
                
                % initial concentration of free receptor
                y0_7 = p.amt_mOR/p.Vd_brain;
                
                % length of simulation
                time_start = 0;
                time_end = 30;

                
                % BASE CASE -- NO NALOXONE
                dose_f = 2; % [mg]       -- amount of one dose of fentanyl
                time_f = 0; % [min]      -- time between repeated fentanyl administrations
                num_f  = 1; % [unitless] -- number of doses of fentanyl given
                
                delay  = time_end;  % [min]  -- no naloxone given so run simulation until time_end
                
                % no naloxone given
                dose_n = 0;  % [mg]       -- amount of one dose of naloxone
                time_n = 0;  % [min]      -- time between repeated naloxone administrations
                num_n  = 0;  % [unitless] -- number of doses of naloxone given

                [base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);

                % receptor occupancy
                base_mu_occ_f = Y(:,8) ./ y0_7;

                dead_base_auc = trapz(base_T,base_mu_occ_f);

                % SIMULATION -- ADD NALOXONE
                delay  = 10; % [min]      -- time between last fentanyl administration and first naloxone administration
                dose_n = 2;  % [mg]       -- amount of one dose of naloxone
                num_n  = 3;  % [unitless] -- number of doses of naloxone given
                time_n = 3;  % [min]      -- time between repeated naloxone administrations
                
                [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);

                % receptor occupancy
                mu_occ_f = Y(:,8) ./ y0_7;
                dead_test_auc = trapz(T, mu_occ_f);
                
                % store data
                popnose_AUC(i) = dead_test_auc; % AUC
                popnose_surv(i) = 1 - (dead_test_auc/dead_base_auc); % survival
                
            end
            
            noses = transpose(noses);
            popnose_AUC = transpose(popnose_AUC);
            popnose_surv = transpose(popnose_surv);
            
            %%% save data for plotting in R
            if(saveFlag)
                new_dir = strcat(root,'/Data/',scenario);
                mkdir(new_dir)
                cd(new_dir)
                
                save("nose_abs.mat", 'noses');
                save("popnose_AUC.mat", 'popnose_AUC');
                save("popnose_surv.mat", 'popnose_surv');
                
                cd(root)
            end

        case "pop_smose"
            %%%%%%%%% NEED TO RUN POP_SMOKE AND POP_NOSE FIRST %%%%%%%%%
            
            % Population PK -- Smoking and Naloxone Absorption

            
            
            %%%%%%%%% use distributions from prev cases %%%%%%%%%
            
           
            %%%%%%%%% run simulation (same as s3_inc_delay_half_naloxone) %%%%%%%%%
            
            %%%% 1 dose of fentanyl (2 mg) given at t=0
            %%%% 10 minutes later, 3 doses of naloxone (2 mg) given 3 min apart
            
            
            % initialize array to store AUCs
            popsmose_AUC = [];
            popsmose_surv = [];

            % for loop to get survival AUC for each patient
            for i = 1:NumberOfSubjects

                % set conditions for each person
                p.amt_mOR = smokers_mOR(i);
                p.abs_n = noses(i);
                
                % initial concentration of free receptor
                y0_7 = p.amt_mOR/p.Vd_brain;
                
                % length of simulation
                time_start = 0;
                time_end = 30;

                
                % BASE CASE -- NO NALOXONE
                dose_f = 2; % [mg]       -- amount of one dose of fentanyl
                time_f = 0; % [min]      -- time between repeated fentanyl administrations
                num_f  = 1; % [unitless] -- number of doses of fentanyl given
                
                delay  = time_end;  % [min]  -- no naloxone given so run simulation until time_end
                
                % no naloxone given
                dose_n = 0;  % [mg]       -- amount of one dose of naloxone
                time_n = 0;  % [min]      -- time between repeated naloxone administrations
                num_n  = 0;  % [unitless] -- number of doses of naloxone given

                [base_T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);

                % receptor occupancy
                base_mu_occ_f = Y(:,8) ./ y0_7;

                dead_base_auc = trapz(base_T,base_mu_occ_f);

                % SIMULATION -- ADD NALOXONE
                delay  = 10; % [min]      -- time between last fentanyl administration and first naloxone administration
                dose_n = 2;  % [mg]       -- amount of one dose of naloxone
                num_n  = 3;  % [unitless] -- number of doses of naloxone given
                time_n = 3;  % [min]      -- time between repeated naloxone administrations
                
                [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);

                % receptor occupancy
                mu_occ_f = Y(:,8) ./ y0_7;
                dead_test_auc = trapz(T,mu_occ_f);
                
                % store data
                popsmose_AUC(i) = dead_test_auc; % AUC
                popsmose_surv(i) = 1 - (dead_test_auc/dead_base_auc); % survival

            end
            
            popsmose_AUC = transpose(popsmose_AUC);
            popsmose_surv = transpose(popsmose_surv);
            
            %%% save data for plotting in R
            if(saveFlag)
                new_dir = strcat(root,'/Data/',scenario);
                mkdir(new_dir)
                cd(new_dir)
                
                save("smose_values.mat", 'smokers_mOR', 'noses');
                save("smose_AUC.mat", 'popnose_AUC');
                save("smose_surv.mat", 'popnose_surv');
                
                cd(root)
            end

        case "sensitivity"
            % Sensitivity Analysis
            
            % ----- Baseline -----
            
            %%%%%%%%% run simulation (same as s3_inc_delay_half_naloxone) %%%%%%%%%
            
            %%%% 1 dose of fentanyl (2 mg) given at t=0
            %%%% 10 minutes later, 3 doses of naloxone (2 mg) given 3 min apart
            
            % initial concentration of free receptor
            y0_7 = p.amt_mOR/p.Vd_brain;
            
            % length of simulation
            time_start = 0;
            time_end = 30;
            
            dose_f = 2; % [mg]       -- amount of one dose of fentanyl
            time_f = 0; % [min]      -- time between repeated fentanyl administrations
            num_f  = 1; % [unitless] -- number of doses of fentanyl given
            
            delay  = 10; % [min]      -- time between last fentanyl administration and first naloxone administration
            dose_n = 2;  % [mg]       -- amount of one dose of naloxone
            num_n  = 3;  % [unitless] -- number of doses of naloxone given
            time_n = 3;  % [min]      -- time between repeated naloxone administrations
            
            % note: b for "baseline"
            [bT, bY, bBalance_f, bBalance_n, bAUC_f, bAUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
            
            baseAUC = [bAUC_f, bAUC_n];
            
            
            % ----- Sensitivity -----
            
            % percentage change parameters by
            delta = 0.1;
            
            % list parameters we want to test
            params = {'p.Vd_brain', 'p.amt_mOR', 'p.plasma_f', ...
                'p.dose_f', 'p.dose_n', 'p.abs_n', ...
                'p.k12_f', 'p.k12_n', ...
                'p.k21_f', 'p.k21_n', ...
                'p.k13_f', 'p.k13_n', ...
                'p.k31_f', 'p.k31_n', ...
                'p.kcl_f', 'p.kcl_n', ...
                'p.kon_f', 'p.kon_n', ...
                'p.koff_f', 'p.koff_n' ...
                'p.Vd_f1', 'p.Vd_f2', ...
                'p.Vd_n1', 'p.Vd_n2'};
            
            % initialize output matrix
            % rows = parameters, columns = fentanyl (1) and naloxone (2)
            AUC = zeros(length(params), 2);
            
            delay  = 10; % [min]      -- time between last fentanyl administration and first naloxone administration
            dose_n = 2;  % [mg]       -- amount of one dose of naloxone
            num_n  = 3;  % [unitless] -- number of doses of naloxone given
            time_n = 3;  % [min]      -- time between repeated naloxone administrations
            
            p.dose_f = dose_f;
            p.dose_n = dose_n;
            
            % For each parameter to be varied
            for z = 1:length(params)
                
                % Change varied parameter by delta
                eval(sprintf('%s = %s * (delta+1);', params{z}, params{z}))
                %eval(sprintf('disp(%s)', params{z}))
                
                
                [T, Y, Balance_f, Balance_n, AUC_f, AUC_n] = fentanyl_naloxone_survive(p, dose_f, time_f, num_f, delay, dose_n, time_n, num_n, time_start, time_end);
                
                % Select the output from the central compartment
                AUC(z,1) = AUC_f;
                AUC(z,2) = AUC_n;
                
                % Change parameter back to original value
                eval(sprintf('%s = %s / (delta+1);', params{z}, params{z}))
                
            end
            
            % Calculate change in AUC
            deltaAUC_f = AUC(:,1) - baseAUC(1);
            deltaAUC_n = AUC(:,2) - baseAUC(2);
            
            % Convert to percentage
            percentAUC_f = deltaAUC_f./baseAUC(1);
            percentAUC_n = deltaAUC_n./baseAUC(2);
            
            % Percentage change in each parameter
            percentP = delta;
            
            % Percent change in AUC divided by percent change in parameter
            AUCchange_f = percentAUC_f/percentP;
            AUCchange_n = percentAUC_n/percentP;
            
            % save data
            if(saveFlag)
                new_dir = strcat(root,'/Data/',scenario);
                mkdir(new_dir);
                cd(new_dir)
                
                save("sensitivity.mat", 'AUCchange_f', 'AUCchange_n');
                
                cd(root)
            end
            
    end
end


%% Function to save data
function save_data(root, scenario, T, Y, Balance_f, Balance_n, AUC_f, AUC_n)

new_dir = strcat(root,'/Data/',scenario);
mkdir(new_dir)
cd(new_dir)

outputdata_T = strcat('T_(',scenario,').mat');
outputdata_Y = strcat('Y_(',scenario,').mat');
outputdata_Balance_f = strcat('Balance_f_(',scenario,').mat');
outputdata_Balance_n = strcat('Balance_n_(',scenario,').mat');
outputdata_AUC_f = strcat('AUC_f_(',scenario,').mat');
outputdata_AUC_n = strcat('AUC_n_(',scenario,').mat');

save(outputdata_T, 'T');
save(outputdata_Y, 'Y');
save(outputdata_Balance_f, 'Balance_f');
save(outputdata_Balance_n, 'Balance_n');
save(outputdata_AUC_f, 'AUC_f');
save(outputdata_AUC_n, 'AUC_n');

cd(root)

end
