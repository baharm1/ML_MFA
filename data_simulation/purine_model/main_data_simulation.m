clc; clear;
addpath('..\general_functions');
%% User-defined inputs
% Model file for serine
isotopmer_model = 'model_serine.xlsx';

% Model file for purines
isotopologue_model = 'model_purine.xlsx';

% specify version to save the simulated data
version = 'gmp_denovo_glioma';
 
% Model file for parameter bounds
param_bound_file = strcat('parameter_bounds_', version, '.xlsx');

% Scaling factor for exchange flux based on Wiechert's bidirectional modeling
F = 200; 

% Number of simulations. Always a factor of 100.
nsim = 50000;

%% Generate equations for serine reactions
% Read the model information
[rxn_serine, input_metabs, metab_char, metab_sym, metab_size,...
    fwd_rxn_idx, bkd_rxn_idx, stoi_full, AM_full, slope_metab, ...
    eq_unbalance, AMeq, Seq, metab_eq, metab_to_remove, pool_metabs, ...
    unlabeled_metabs, balance_metabs, n_flux_param, n_c_param, ...
    Aeq, M0, varFlag] = import_stoich_AM(isotopmer_model);
disp('Stoichiometry and atom-mapping of serine model imported from files');

% Generate ODEs for serine model
[remove_idx, nterms] = ODE_fn_gen(stoi_full, metab_char, AM_full, ...
    metab_size, unlabeled_metabs, input_metabs, balance_metabs, ...
    n_flux_param, n_c_param, fwd_rxn_idx, bkd_rxn_idx, F, ...
    metab_eq, AMeq, Seq, eq_unbalance);

%% Read model information for purine reactions
[Sp, purine_list, purine_rxns] = read_purine_model(isotopologue_model);

% combine the reaction list
rxn_full = [rxn_serine; purine_rxns];
metabs = [balance_metabs, purine_list];

%% Read the reaction bounds 
[flux_bounds, rxn_id] = xlsread(param_bound_file, 'flux', '', 'basic');
[pool_bounds, pool_id] = xlsread(param_bound_file, 'pool', '', 'basic');

% Re-arrange the bounds
[~, index] = ismember(rxn_full, rxn_id);
for i = 1:numel(index) - 1
    if(index(i + 1) == index(i))
        index(i + 1) = index(i + 1) + 1;
    end
end
flux_bounds = flux_bounds(index,:);

% order of reactions, pay attention to purine reactions. used in simulate function
rxn_id(index) % equals to rxn_full
writecell(rxn_id(index), strcat('rxn_ids_', version, '.txt'), 'Delimiter', 'tab');
writematrix(flux_bounds, strcat('flux_bounds_', version, '.txt'), 'Delimiter', 'tab');

clear index
% order of metabs, pay attention to purines. used in simulate function
[~, index] = ismember(metabs, pool_id);
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
ser1 = 0.04 + (0.2 - 0.04) * rand(1, nsim);
ser2 = 0.08 * rand(1, nsim);
ser3 = 0.12 * rand(1, nsim);

r5p1 = 0.06 * rand(1, nsim);
r5p2 = 0.08 * rand(1, nsim);
r5p3 = 0.12 * rand(1, nsim);
r5p4 = 0.03 * rand(1, nsim);
r5p5 = 0.08 * rand(1, nsim);

%% Ensure that the fluxes are stoichiometrically balanced
% serine fluxes
fluxnew = zeros(numel(rxn_serine), nsim);

options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10000, ...
    'ConstraintTolerance', 1e-4);
exitFlagSer = zeros(nsim, 1);
fvalSer = zeros(nsim, 1);
outputSer = [];

parfor (i = 1:nsim, 10)
    [fluxnew(:,i), fvalSer(i), exitFlagSer(i), output] = fmincon(@(x) ...
        sumsqr(1 - x ./ flux0(1:numel(rxn_serine), i)), ...
        flux0(1:numel(rxn_serine), i), [], [], Aeq, ...
        zeros(size(Aeq, 1), 1), flux_bounds(1:numel(rxn_serine), 1), ...
        flux_bounds(1:numel(rxn_serine), 2), [], options) ;
    outputSer = [outputSer, output];
end

outputSer = struct2table(outputSer);
successful_sims_ser = outputSer.constrviolation < 1e-4 & exitFlagSer >= 0;
sum(~successful_sims_ser)
flux0(1:numel(rxn_serine),:) = fluxnew;

%% check if any flux estimation does not satisfy a constraint
notSuccessful_ser = find(~successful_sims_ser);

%% rerun for not successful instances
i = 1;
while (i <= length(notSuccessful_ser))
    j = notSuccessful_ser(i);
    flux0_2 = flux_bounds(1:numel(rxn_serine), 1) + rand(length(rxn_serine), 1) .* ...
        (flux_bounds(1:numel(rxn_serine), 2) - flux_bounds(1:numel(rxn_serine), 1));

    [fluxnew(:, j), fvalSer(j), exitFlagSer(j), output] = fmincon(@(x)sumsqr(1 - x ./ flux0_2), ...
        flux0_2, [], [], Aeq, zeros(size(Aeq, 1), 1), ...
        flux_bounds(1:numel(rxn_serine), 1), ...
        flux_bounds(1:numel(rxn_serine), 2), [], options) ;
    if (output.constrviolation < 1e-4 && exitFlagSer(j) >= 0)
        i = i + 1;
        continue;
    end
end
flux0(1:numel(rxn_serine), :) = fluxnew;

% purine fluxes
% Also ensure that the ratio of de-novo GMP synthesis is uniformly distributed
% add constraints for uniform distribution of GMP synthesis
ratio_values = linspace(0.01, 0.99, 100); % cannot set to 1
ratio_values = repmat(ratio_values, 1, nsim / 100);
ratio_values = ratio_values(randperm(length(ratio_values))); % shuffle

fluxnew = zeros(size(flux0, 1) - numel(rxn_serine), nsim);
exitFlagPur = zeros(nsim, 1);
fvalPur = zeros(nsim, 1);
outputPur = [];

Aend = zeros(nsim, size(Sp, 2));

GMP_salvage = find(ismember(purine_rxns, 'GUA + R5P == GMP'));
GMP_denovo = find(ismember(purine_rxns, 'IMP == GMP'));

Aend(:, GMP_denovo) = 1;
Aend(:, GMP_salvage) = (ratio_values ./ (ratio_values - 1))';

parfor (i = 1:nsim, 10)
    [fluxnew(:, i), fvalPur(i), exitFlagPur(i), output] = fmincon(@(x) ...
        sumsqr(1 - x ./ flux0(numel(rxn_serine) + 1:end, i)), ...
        flux0(numel(rxn_serine) + 1:end, i), [], [], [Sp; Aend(i, :)], ...
        zeros(size(Sp, 1) + 1, 1), flux_bounds(numel(rxn_serine) + 1:end, 1), ...
        flux_bounds(numel(rxn_serine) + 1:end, 2), [], options);
    outputPur = [outputPur, output];    
end
outputPur = struct2table(outputPur);
successful_sims_pur = outputPur.constrviolation < 1e-4 & exitFlagPur >= 0;
sum(~successful_sims_pur)

flux0(numel(rxn_serine) + 1:end, :) = fluxnew;

%% check if any flux estimation does not satisfy a constraint
notSuccessful = find(~successful_sims_pur);

ratio_values(notSuccessful)

%% rerun for not successful instances
i = 1;
while (i <= length(notSuccessful))
    j = notSuccessful(i);
    flux0_2 = flux_bounds(numel(rxn_serine) + 1:end, 1) + rand(length(purine_rxns), 1) .* ...
        (flux_bounds(numel(rxn_serine) + 1:end, 2) - flux_bounds(numel(rxn_serine) + 1:end, 1));

    [fluxnew(:, j), fvalPur(j), exitFlagPur(j), output] = fmincon(@(x)sumsqr(1 - x ./ flux0_2), ...
        flux0_2, [], [], [Sp; Aend(j, :)], zeros(size(Sp, 1) + 1, 1), ...
        flux_bounds(numel(rxn_serine) + 1:end, 1), ...
        flux_bounds(numel(rxn_serine) + 1:end, 2), [], options);

    if (output.constrviolation < 1e-6 && exitFlagPur(j) >= 0)
        i = i + 1;
        continue;
    end
end
flux0(numel(rxn_serine) + 1:end, :) = fluxnew;

%% save fluxes to later plot the distributions
writematrix(flux0, strcat('flux_', version, '.txt'), 'Delimiter', 'tab');

%% Simulate
% Create M0
% add isotopologues of purines to isotopomers of serine
M0p = zeros(6, 1); % purine model considers 5 carbons for each metabolite M+0 to M+5
M0p(1) = 1; % set M+0 to 1
% adding isotopologues of IMP, GMP, GDP, GUANOSINE, INOSINE, AMP (balanced metabolites)
M0 = [M0; repmat(M0p, 6, 1)];

% Create the matrix to insert input metabolite IDV
J = zeros(nterms + numel(remove_idx), nterms);
j = 1;
for i = 1:size(J, 1)
    if(~ismember(remove_idx, i))
        J(i, j) = 1;
        j = j + 1;
    end
end

%% Create a table to save simulation parameters
simulation_parameters = array2table(zeros(0, size(flux0, 1) + size(pool0, 1) + 8));
names = [rxn_full', metabs, "SER1", "SER2", "SER3", "R5P1", "R5P2", "R5P3", "R5P4", "R5P5"];
names(bkd_rxn_idx) = strcat(names(bkd_rxn_idx), '_reverse');
simulation_parameters.Properties.VariableNames = names;

% Create table to save the metabolite MIDs
k = 3; % starts from the third column, the first two columns are index and time
glycine_index = [k:k+2];
k = k + 2 + 1;
cthf_index = [k:k+1];
k = k + 1 + 1;
amp_index = [k:k+5];
k = k + 5 + 1;
gdp_index = [k:k+5];
k = k + 5 + 1;
gmp_index = [k:k+5];
k = k + 5 + 1;
guanosine_index = [k:k+5];
k = k + 5 + 1;
imp_index = [k:k+5];
k = k + 5 + 1;
inosine_index = [k:k+5];
k = k + 5;

simulated_data = array2table(zeros(0, k));
clear names;
names = ["index", "time", strcat(repmat("GLY", 1, 3), ["0", "1", "2"])];
names = [names, strcat(repmat("MTHF", 1, 2), ["0", "1"])];
names = [names, strcat(repmat("AMP", 1, 6), ["0", "1", "2", "3", "4", "5"])];
names = [names, strcat(repmat("GDP", 1, 6), ["0", "1", "2", "3", "4", "5"])];
names = [names, strcat(repmat("GMP", 1, 6), ["0", "1", "2", "3", "4", "5"])];
names = [names, strcat(repmat("GUO", 1, 6), ["0", "1", "2", "3", "4", "5"])];
names = [names, strcat(repmat("IMP", 1, 6), ["0", "1", "2", "3", "4", "5"])];
names = [names, strcat(repmat("INO", 1, 6), ["0", "1", "2", "3", "4", "5"])];
simulated_data.Properties.VariableNames = names;

varFlag2 = varFlag;
[~, ip_idx] = ismember(input_metabs, metab_char); 
varFlag2(ip_idx) = 1;
clear m
m = find(ismember(metab_char, 'GLYg'));
M_glycine = uIDVindex(m, metab_size, varFlag2, 0);
map_glycine = idv_to_mid(2);
clear m
m = find(ismember(metab_char, 'MTHFg'));
M_cthf = uIDVindex(m, metab_size, varFlag2, 0);
clear m

tspan = [0:0.1:4];
num_sim = sum(tspan >= 0);
for j = 1:nsim/10
    clear M v c M_save save_matrix

    parfor (i = (10 * (j - 1) + 1):(10 * j), 10)
        v = flux0(:, i);
        c = pool0(:, i);
        M{i} = simulate(v, c, n_flux_param, n_c_param, ...
            ser1(i), ser2(i), ser3(i), r5p1(i), r5p2(i), r5p3(i), r5p4(i), r5p5(i), ...
            M0, J, tspan, metab_char, metab_size, varFlag);
    end
    
    for i = (10 * (j - 1) + 1):(10 * j)
        clear Mi

        % Save the simulation parameters
        v = flux0(:, i);
        c = pool0(:, i);
        simulation_parameters(i, :) = num2cell([v', c', ...
            ser1(i), ser2(i), ser3(i), r5p1(i), r5p2(i), r5p3(i), r5p4(i), r5p5(i)]);
    
        % Save the simulation data
        save_matrix = zeros(num_sim, k);
        Mi = M{i};
        M_save = Mi(tspan >= 0, :);
        save_matrix(:, 1) = i;
        save_matrix(:, 2) = tspan(tspan >= 0)';
		% save purine MIDs, 6 balanced metabolites and 6 MIDs --> 36 total
        save_matrix(:, end - 36 + 1:end) = M_save(:, end - 36 + 1:end);
        save_matrix(:, glycine_index) = M_save(:, M_glycine(1):M_glycine(2)) * map_glycine';
        save_matrix(:, cthf_index) = M_save(:, M_cthf(1):M_cthf(2));
      
        simulated_data(size(simulated_data, 1) + 1:size(simulated_data, 1) + num_sim, :) = num2cell(save_matrix);
    end

end

% Save the files
writetable(simulation_parameters, strcat('simulation_parameters_', version, '.csv'));
writetable(simulated_data, strcat('simulated_data_', version, '.csv'));
save(strcat('data_', version, '.mat'));
