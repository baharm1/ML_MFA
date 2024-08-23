function [obj] = calc_obj(x, MID_metab, MID, SD, metab_char, metab_size, ...
    nflux, varFlag)
    
obj = 0;
x = x(nflux + 1:end);

% Estimate the objective for the metabolites
% with experimentally measured MIDs
metab_list = unique(MID_metab, 'stable');
sim_mid_data = zeros(size(MID));
idx1 = 1;
for i = 1:numel(metab_list)
    
    clear m index sim_idvs map sim_mids
    % get index of the metabolite in the list
    m = find(ismember(metab_char, metab_list(i)));
    % Get the simulated IDV values
    index = uIDVindex(m, metab_size, varFlag, 0);
    sim_idvs = x(index(1):index(2));
    % Convert IDV to MID
    map = create_mapping_matrix(metab_size(m));
    sim_mids = map * sim_idvs;
    % Save the MID values
    idx2 = idx1 + numel(sim_mids) - 1;
    sim_mid_data(idx1:idx2) = sim_mids;
    idx1 = idx2 + 1;
end
if(idx2 ~= numel(MID))
    error('ERROR: The size of the simulated MID vector is different from the size of the experimental MID vector.');
end

obj = obj + sumsqr((MID - sim_mid_data) ./ SD);
end

function map = create_mapping_matrix(s)

map = zeros(s + 1, 2 ^ s);
convert = sum(dec2bin(0:(2 ^ s - 1)) == '1', 2);
for k = 1:s + 1
    map(k, convert == k - 1) = 1;
end
end