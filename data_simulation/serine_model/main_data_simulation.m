clc; clear;
addpath('..\general_functions');
%% User-defined inputs
% Model file for serine
isotopmer_model = 'model_serine.xlsx';

% specify version to save the simulated data
version = 'serine_model';
 
% Model file for parameter bounds
param_bound_file = strcat('parameter_bounds_', version, '.xlsx');

% Scaling factor for exchange flux based on Wiechert's bidirectional modeling
F = 200; 

% Number of simulations. Always a factor of 100.
nsim = 50000;

%% Generate equations for serine reactions
% Read the model information
[rxn_serine, input_metabs, metab_char, metab_sym, metab_size, ...
    fwd_rxn_idx, bkd_rxn_idx, stoi_full, AM_full, slope_metab, ...
    eq_unbalance, AMeq, Seq, metab_eq, metab_to_remove, pool_metabs, ...
    unlabeled_metabs, balance_metabs, n_flux_param, n_c_param, ...
    Aeq, M0, varFlag] = import_stoich_AM(isotopmer_model);
disp('Stoichiometry and atom-mapping of serine model imported from files');

% Generate ODEs for serine model
[remove_idx,nterms] = ODE_fn_gen(stoi_full, metab_char, AM_full, ...
    metab_size, unlabeled_metabs, input_metabs, balance_metabs, ...
    n_flux_param, n_c_param, fwd_rxn_idx, bkd_rxn_idx, F, ...
    metab_eq, AMeq, Seq, eq_unbalance);

rxn_full = rxn_serine;
metabs = balance_metabs;

%% Read the reaction bounds 
[flux_bounds, rxn_id] = xlsread(param_bound_file, 'flux', '', 'basic');
[pool_bounds, pool_id] = xlsread(param_bound_file, 'pool', '', 'basic');

% Re-arrange the bounds
[~,index] = ismember(rxn_full, rxn_id);
for i = 1:numel(index) - 1
    if(index(i + 1) == index(i))
        index(i + 1) = index(i + 1) + 1;
    end
end
flux_bounds = flux_bounds(index, :);

writecell(rxn_id(index), strcat('rxn_ids_', version, '.txt'), 'Delimiter', 'tab');
writematrix(flux_bounds, strcat('flux_bounds_', version, '.txt'), 'Delimiter', 'tab');

clear index
% order of metabs
[~,index] = ismember(metabs, pool_id);
pool_bounds = pool_bounds(index, :);
pool_id(index)
writecell(pool_id(index), strcat('pool_ids_', version, '.txt'), 'Delimiter', 'tab');
writematrix(pool_bounds, strcat('pool_bounds_', version, '.txt'), 'Delimiter', 'tab');
%% Sample from random flux distributions between the bounds
flux0 = repmat(flux_bounds(:, 1), 1, nsim) + ...
    rand(size(flux_bounds, 1), nsim) .* repmat(flux_bounds(:, 2) - ...
    flux_bounds(:, 1), 1, nsim);
pool0 = repmat(pool_bounds(:, 1), 1, nsim) + ...
    rand(size(pool_bounds, 1), nsim) .* repmat(pool_bounds(:, 2) - ...
    pool_bounds(:, 1), 1, nsim);

% define bounds for input metabolite MIDs based on patient data
serp10 = 0.2 + (0.7-0.2) * rand(1, nsim);
serp2 = 0.02 * rand(1, nsim);
serp3 = 0.01 * rand(1, nsim);

pgc3 = 0.02 + (0.25 - 0.02) * rand(1, nsim); 
pgc1 = 0.03 * rand(1, nsim);
pgc2 = 0.015 * rand(1, nsim);

pg3 = 0.02 + (0.17 - 0.02) * rand(1, nsim); 
pg1 = 0.06 * rand(1, nsim);
pg2 = 0.035 * rand(1, nsim);

%% Ensure that the fluxes are stoichiometrically balanced
% serine fluxes
plasma_uptake = find(ismember(rxn_serine, 'SERp == SERg'));
ser_synthesis = find(ismember(rxn_serine, 'PGg == SERg'));
tme_uptake = find(ismember(rxn_serine, 'SERc == SERg'));

plasma_uptake_c = find(ismember(rxn_serine, 'SERp == SERc'));
ser_synthesis_c = find(ismember(rxn_serine, 'PGc == SERc'));

% PGg to total
ratio_values = linspace(0.01, 0.99, 100);
ratio_values = repmat(ratio_values, 1, nsim / 100);
ratio_values = ratio_values(randperm(length(ratio_values))); % shuffle

fluxnew = zeros(numel(rxn_serine), nsim);
options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10000, ...
    'ConstraintTolerance', 1e-4);
exitFlagSer = zeros(nsim, 1);
fvalSer = zeros(nsim, 1);
outputSer = [];

% Add constraints for uniform distribution of glucose-derived serine synthesis, 
% TME-derived serine uptake, and plasma serine uptake
% This file ran three times each time with one of the following
% constraints. 
Aend = zeros(nsim,size(Aeq, 2));

% simulate data for relative glucose-derived serine synthesis flux in
% glioma
Aend(:, ser_synthesis) = 1;
Aend(:, plasma_uptake) = ratio_values ./ (ratio_values - 1);
Aend(:, tme_uptake) = ratio_values ./ (ratio_values - 1);

% simulate data for relative plasma serine uptake flux in glioma
% Aend(:, plasma_uptake) = 1;
% Aend(:, ser_synthesis) = ratio_values./(ratio_values-1);
% Aend(:, tme_uptake) = ratio_values./(ratio_values-1);

% simulate data for relative glucose-derived serine synthesis flux in
% cortex
% Aend(:, ser_synthesis_c) = 1;
% Aend(:, plasma_uptake_c) = ratio_values./(ratio_values-1);

parfor (i = 1:nsim, 10)
    [fluxnew(:,i), fvalSer(i), exitFlagSer(i), output] = fmincon(@(x) ...
        sumsqr(1 - x ./ flux0(1:numel(rxn_serine), i)), ...
        flux0(1:numel(rxn_serine), i), [], [], [Aeq; Aend(i, :)], ... 
        zeros(size(Aeq, 1) + 1, 1), flux_bounds(1:numel(rxn_serine), 1), ...
        flux_bounds(1:numel(rxn_serine), 2), [], options) ;
    outputSer = [outputSer, output];
end

outputSer = struct2table(outputSer);
successful_sims_ser = outputSer.constrviolation < 1e-4 & exitFlagSer >= 0;
sum(~successful_sims_ser)
flux0(1:numel(rxn_serine), :) = fluxnew;

%% check if any flux estimation does not satisfy a constraint
notSuccessful_ser = find(~successful_sims_ser);

ratio_values(notSuccessful_ser)

%% rerun for not successful instances
i = 1;
while (i <= length(notSuccessful_ser))
    j = notSuccessful_ser(i);
    flux0_2 = flux_bounds(1:numel(rxn_serine),1) + rand(length(rxn_serine), 1) .* ...
        (flux_bounds(1:numel(rxn_serine),2)-flux_bounds(1:numel(rxn_serine),1));

    [fluxnew(:,j), fvalSer(j), exitFlagSer(j), output] = fmincon(@(x)sumsqr(1 - x ./ flux0_2), ...
        flux0_2, [], [], [Aeq; Aend(i, :)], zeros(size(Aeq, 1) + 1, 1), ...
        flux_bounds(1:numel(rxn_serine), 1), ...
        flux_bounds(1:numel(rxn_serine), 2), [], options) ;
    if (output.constrviolation < 1e-4 && exitFlagSer(j) >= 0)
        i = i + 1;
        continue;
    end
end
flux0(1:numel(rxn_serine), :) = fluxnew;

%% save fluxes to later plot the distributions
writematrix(flux0, strcat('flux_', version, '.txt'), 'Delimiter', 'tab');

%% Simulate

% Create the matrix to insert circulating metabolite IDV
J = zeros(nterms + numel(remove_idx), nterms);
j = 1;
for i = 1:size(J, 1)
    if(~ismember(remove_idx, i))
        J(i, j) = 1;
        j = j + 1;
    end
end

%% Create a table to save simulation parameters
simulation_parameters = array2table(zeros(0, size(flux0, 1) + size(pool0, 1) + 9));
names = [rxn_full', metabs, "SERp1", "SERp2", "SERp3", "PGc1", "PGc2", "PGc3", "PGg1", "PGg2", "PGg3"];
names(bkd_rxn_idx) = strcat(names(bkd_rxn_idx), '_reverse');
simulation_parameters.Properties.VariableNames = names;

% Create table to save the metabolite MIDs
k = 3; % starts from the third column, the first two columns are index and time
pg_index = [k:k+3];
k = k + 3 + 1;
serine_index = [k:k+3];
k = k + 3 + 1;
glycine_index = [k:k+2];
k = k + 2 + 1;
mthf_index = [k:k+1];
k = k + 1 + 1;

pg_c_index = [k:k+3];
k = k + 3 + 1;
serine_c_index = [k:k+3];
k = k + 3 + 1;
glycine_c_index = [k:k+2];
k = k + 2 + 1;
mthf_c_index = [k:k+1];
k = k + 1;

simulated_data = array2table(zeros(0, k));
clear names;
names = ["index", "time"];
names = [names, strcat(repmat("PG", 1, 4), ["0", "1", "2", "3"])];
names = [names, strcat(repmat("SER", 1, 4), ["0", "1", "2", "3"])];
names = [names, strcat(repmat("GLY", 1, 3), ["0", "1", "2"])];
names = [names, strcat(repmat("MTHF", 1, 2), ["0", "1"])];
names = [names, strcat(repmat("PGc", 1, 4), ["0", "1", "2", "3"])];
names = [names, strcat(repmat("SERc", 1, 4), ["0", "1", "2", "3"])];
names = [names, strcat(repmat("GLYc", 1, 3), ["0", "1", "2"])];
names = [names, strcat(repmat("MTHFc", 1, 2), ["0", "1"])];
simulated_data.Properties.VariableNames = names;

varFlag2 = varFlag;
[~, ip_idx] = ismember(input_metabs, metab_char); 
varFlag2(ip_idx) = 1;
clear m
m = find(ismember(metab_char, 'PGg'));
M_pg = uIDVindex(m, metab_size, varFlag2, 0);
map_pg = idv_to_mid(3);
clear m
m = find(ismember(metab_char, 'SERg'));
M_serine = uIDVindex(m, metab_size, varFlag2, 0);
map_serine = idv_to_mid(3);
clear m
m = find(ismember(metab_char, 'GLYg'));
M_glycine = uIDVindex(m, metab_size, varFlag2, 0);
map_glycine = idv_to_mid(2);
clear m
m = find(ismember(metab_char, 'MTHFg'));
M_mthf = uIDVindex(m, metab_size, varFlag2, 0);

clear m
m = find(ismember(metab_char, 'PGc'));
M_pg_c = uIDVindex(m, metab_size, varFlag2, 0);
map_pg_c = idv_to_mid(3);
clear m
m = find(ismember(metab_char, 'SERc'));
M_serine_c = uIDVindex(m, metab_size, varFlag2, 0);
map_serine_c = idv_to_mid(3);
clear m
m = find(ismember(metab_char, 'GLYc'));
M_glycine_c = uIDVindex(m, metab_size, varFlag2, 0);
map_glycine_c = idv_to_mid(2);
clear m
m = find(ismember(metab_char, 'MTHFc'));
M_mthf_c = uIDVindex(m, metab_size, varFlag2, 0);
clear m

tspan = [0:0.1:4];
num_sim = sum(tspan >= 0);
for j = 1:nsim/10
    clear M v c M_save save_matrix

    parfor (i = (10 * (j - 1) + 1):(10 * j), 10)
        v = flux0(:, i);
        c = pool0(:, i);
        M{i} = simulate(v, c, n_flux_param, n_c_param,...
            serp10(i), serp2(i), serp3(i), pgc1(i), pgc2(i), pgc3(i), ...
            pg1(i), pg2(i), pg3(i), M0, J, ...
            tspan, metab_char, metab_size, varFlag);
    end
    
    for i = (10 * (j - 1) + 1):(10 * j)
        clear Mi

        % Save the simulation parameters
        v = flux0(:, i);
        c = pool0(:, i);
        simulation_parameters(i, :) = num2cell([v', c', serp10(i), serp2(i), ...
            serp3(i), pgc1(i), pgc2(i), pgc3(i), pg1(i), pg2(i), pg3(i)]);
    
        % Save the simulation data
        save_matrix = zeros(num_sim, k);
        Mi = M{i};
        M_save = Mi(tspan >= 0, :);
        save_matrix(:, 1) = i;
        save_matrix(:, 2) = tspan(tspan >= 0)';
        save_matrix(:, pg_index) = M_save(:, M_pg(1):M_pg(2)) * map_pg';
        save_matrix(:, serine_index) = M_save(:, M_serine(1):M_serine(2)) * map_serine';
        save_matrix(:, glycine_index) = M_save(:, M_glycine(1):M_glycine(2)) * map_glycine';
        save_matrix(:, mthf_index) = M_save(:, M_mthf(1):M_mthf(2));

        save_matrix(:, pg_c_index) = M_save(:, M_pg_c(1):M_pg_c(2)) * map_pg_c';
        save_matrix(:, serine_c_index) = M_save(:, M_serine_c(1):M_serine_c(2)) * map_serine_c';
        save_matrix(:, glycine_c_index) = M_save(:, M_glycine_c(1):M_glycine_c(2)) * map_glycine_c';
        save_matrix(:, mthf_c_index) = M_save(:, M_mthf_c(1):M_mthf_c(2));
                
        simulated_data(size(simulated_data, 1) + 1:size(simulated_data, 1) + num_sim, :) = num2cell(save_matrix);
    end

end

% Save the files
writetable(simulation_parameters, strcat('simulation_parameters_', version, '.csv'));
writetable(simulated_data, strcat('simulated_data_', version, '.csv'));
save(strcat('data_', version, '.mat'));

