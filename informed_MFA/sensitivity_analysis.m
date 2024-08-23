clc; clear;
addpath('.\general_functions');
ML_files = dir("output_files");
ML_files = struct2table(ML_files);
ML_files = ML_files(3:end, :);
ML_patient_sites = ML_files.name;
%%
for ps_ML = 1:2 %1:length(ML_patient_sites)

ps_dir = "output_files\" + ML_patient_sites(ps_ML);
load(ps_dir + "\data.mat");

%% Perform parameter sensitivity analysis by Monte Carlo Sampling

% Generate sampled data
mid_mc = mcb_sampler(nmc, MID_metab, MID, SD);

% Optimize for sampled data
for i = 1:nmc
    clear mc_objective
    mc_objective = @(x)calc_obj(x, MID_metab, mid_mc(:, i), SD, ...
        metab_char, metab_size, nflux, varFlag);
    [xmc(:, i), fmc(i), exitmc(i)] = knitro_nlp(mc_objective, xopt_min, ...
        A, b, Aeq, beq, lb, ub, get_NL_cons, [], opt_options, []);
    save(ps_dir + '/data_mc.mat');
end

% Determine and save the flux confidence bounds
alpha = 100 * (1 - conf) / 2;
ubflux = zeros(nflux, 1);
lbflux = ubflux;
for i = 1:nflux
    flux_values = xmc(i, :);
    ubflux(i) = prctile(flux_values, 100 - alpha);
    lbflux(i) = prctile(flux_values, alpha);
end
bounds_table = table(rxn_full, xopt_min(1:nflux), lbflux, ubflux);
bounds_table.Properties.VariableNames = ["Reaction", "Flux", "Lower_CI", "Upper_CI"];
writetable(bounds_table, ps_dir + '/flux_results_ML.xlsx', "Sheet", "flux_bounds");

xmc_table = table(rxn_full, xmc(1:nflux, :));
xmc_table.Properties.VariableNames = ["Reaction", "Flux"];
writetable(xmc_table, ps_dir + '/flux_results_ML.xlsx', "Sheet", "flux");

save(ps_dir +'/data_mc.mat');

%%
xmc_idvs = xmc(nflux+1:end, :);
metab_list = unique(MID_metab, 'stable');

MID_names = [];
xmc_mids = [];
for i = 1:numel(metab_list)
    
    clear m index sim_idvs map sim_mids values error
    % get index of the metabolite in the list
    m = find(ismember(metab_char, metab_list(i)));
    % Get the simulated IDV values
    index = uIDVindex(m,metab_size, varFlag, 0);
    sim_idvs = xmc_idvs(index(1):index(2), :);
    % Convert IDV to MID
    map = create_mapping_matrix(metab_size(m));
    sim_mids = map * sim_idvs;
    MID_names = [MID_names; cellstr(strcat(char(metab_list(i)), ...
        num2str((0:metab_size(m))')))];
    xmc_mids = [xmc_mids; sim_mids];

end  

mid_table = table(MID_names, xmc_mids);
mid_table.Properties.VariableNames = ["MID", "Value"];
writetable(mid_table, ps_dir + '/flux_results_ML.xlsx', "Sheet", "MID");

save(ps_dir + '/data_mc.mat');

end

function map = create_mapping_matrix(s)

map = zeros(s + 1, 2 ^ s);
convert = sum(dec2bin(0:(2 ^ s - 1)) == '1', 2);
for k = 1:s + 1
    map(k, convert == k - 1) = 1;
end
end