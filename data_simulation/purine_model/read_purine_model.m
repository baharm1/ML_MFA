function [Sp,purine_list,purine_rxns] = read_purine_model(file2)

[~,rxns] = xlsread(file2, 'rxns', '', 'basic');
[~,input_metabs] = xlsread(file2, 'input_metabs', '', 'basic');

% Create stoichiometric matrix and metabolite list

% Convert the reactions to symbolic form
rxns = rxns(:, 2);
for i = 1:size(rxns, 1)
    text2(i) = str2sym(rxns(i, 1));
end

% Extract metabolites
metab_sym = symvar(text2);

% Create stoichiometric matrix
S = equationsToMatrix(text2);
S = double(-S');

% only keep the metabs to be balanced
index = ismember(string(metab_sym), string(input_metabs));
S = S(~index, :);
Sp = S;

% Save model properties
purine_list = string(metab_sym(~index));
purine_rxns = rxns;

end