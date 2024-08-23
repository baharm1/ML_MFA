function mid_mc = mcb_sampler(nmc, MID_metab, MID, SD)

% Create random variables sampled from a truncated normal distribution
% This is the inverse transform method
% Sources: http://web.michaelchughes.com/research/sampling-from-truncated-normal + wikipedia

% Truncated normal distribution between 0 and 1
alpha = -MID ./ SD;
beta = (1 - MID) ./ SD;

a = normcdf(alpha);
b = normcdf(beta);
u = a + (b - a) .* unifrnd(0, 1, numel(MID), nmc);
x = norminv(u);

% Sample from the distribution
mid_mc = repmat(MID, 1, nmc) + repmat(SD, 1, nmc) .* x;

% MIDs should sum up to 1
% Normalize the data 
metabs = unique(MID_metab, 'stable');
for i = 1:numel(metabs)
   clear metab_idx 
   metab_idx = find(ismember(MID_metab, metabs(i)));
   metab_data = mid_mc(metab_idx, :);
   mid_mc(metab_idx, :) = metab_data ./ (repmat(sum(metab_data), ...
       size(metab_data, 1), 1));

end
