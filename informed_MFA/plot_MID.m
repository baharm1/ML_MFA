function out = plot_MID(x, metab_size, nflux, varFlag, metab_char, ...
    unlabeled_metabs, MID_metab, MID, SD)

out = [];

x = x(nflux + 1:end);

global filepath_stored_gl
mkdir([filepath_stored_gl, '/plots'])

% Plot MIDs of metabolites with measured experimental values
metab_list = unique(MID_metab, 'stable');
for i = 1:numel(metab_list)
    
    clear m index sim_idvs map sim_mids values error
    % get index of the metabolite in the list
    m = find(ismember(metab_char, metab_list(i)));
    % Get the simulated IDV values
    index = uIDVindex(m, metab_size, varFlag, 0);
    sim_idvs = x(index(1):index(2));
    % Convert IDV to MID
    map = create_mapping_matrix(metab_size(m));
    sim_mids = map * sim_idvs;
    writematrix(sim_mids, [filepath_stored_gl, '/flux_results.xlsx'], ...
        "Sheet", char(metab_list(i)));

    % Create plots
    values = [MID(ismember(MID_metab, metab_list(i)))'; sim_mids'];
    error = [SD(ismember(MID_metab, metab_list(i)))'; zeros(1, metab_size(m) + 1)];
    figure;
    b = bar(0:metab_size(m), values, 'grouped');
    hold on
    pos = zeros(size(values));
    for k = 1:size(pos, 1)
        pos(k, :) = b(k).XEndPoints;
    end
    errorbar(pos, values, error, 'Color', 'black', 'LineStyle', 'none'); 
    title(metab_list(i))
    ylabel('MID')
    xlabel('M+i')
    legend('Experimental', 'Simulated')
    hold off
    filepath = [filepath_stored_gl, '/plots/MID_', char(metab_char(m)), '.png'];
    saveas(gcf, filepath);
    close
   
end


% Plot MIDs of labeled metabolites without experimental data
clear metab_list
metab_list = setdiff(metab_char, [unique(MID_metab)', unlabeled_metabs]);
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
    writematrix(sim_mids, [filepath_stored_gl, '/flux_results.xlsx'], ...
        "Sheet", char(metab_list(i)));
    % Create plots
    figure;
    bar(0:metab_size(m), sim_mids);
    hold on
    title(metab_list(i))
    ylabel('MID')
    xlabel('M+i')
    legend('Simulated')
    hold off
    filepath = [filepath_stored_gl, '/plots/MID_', char(metab_char(m)), '.png'];
    saveas(gcf, filepath);
    close
   
end

end

function map = create_mapping_matrix(s)

map = zeros(s + 1, 2 ^ s);
convert = sum(dec2bin(0:(2 ^ s - 1)) == '1', 2);
for k = 1:s+1
    map(k, convert == k - 1) = 1;
end
end
