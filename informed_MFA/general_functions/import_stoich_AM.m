function [Sfull, rxn_full, irrev_rxn, exc_rxn, rev_rxn, metab_char, ...
    AM_full, metab_size, input_metabs, metabs_remove] = import_stoich_AM(file1)

%% Importing stoichiometry
[~, text] = xlsread(file1, 'rxns', '', 'basic');
% Remove empty columns
remove_col = sum(~cellfun(@isempty, text));
text = text(:, remove_col > 0);
rxn_flag = cat(1, text{:, 1});
irrev_rxn = find(rxn_flag == 'I');
exc_rxn = find(rxn_flag == 'E');
rev_rxn = find(rxn_flag == 'R');
rxn = text(:, 2);

%% Convert equations from text file into MATLAB symbolic form
n_rxn = size(rxn, 1);
for i = 1:n_rxn
    text2(i) = str2sym(rxn(i));
end

% Extract unique symbolic variables
metab_sym = symvar(text2);
n_met = length(metab_sym);

% Convert symbolic array of variables to strings
for i = 1:n_met
    metab_char{i} = char(metab_sym(i));
end

% Convert symbolic equations to stoichiometric matrix
nRev2 = length(rev_rxn) * 2;
fwd_rxn_idx = [1:2:nRev2-1];
bkd_rxn_idx = [2:2:nRev2];
flag_stoi_matrix = equationsToMatrix(text2);
flag_stoi_matrix = double(-flag_stoi_matrix');
stoi_rev = flag_stoi_matrix(:, rev_rxn);
rxn_rev = rxn(rev_rxn);
stoi_irrev = flag_stoi_matrix(:, irrev_rxn);
rxn_irrev = rxn(irrev_rxn);
stoi_cell(:, fwd_rxn_idx) = stoi_rev;
rxn_cell(fwd_rxn_idx, 1) = rxn_rev;
stoi_cell(:, bkd_rxn_idx) = -stoi_rev;
rxn_cell(bkd_rxn_idx, 1) = rxn_rev;
stoi_full = [stoi_cell, stoi_irrev, flag_stoi_matrix(:, exc_rxn)];
rxn_full = [rxn_cell; rxn_irrev; rxn(exc_rxn)];

Sfull = stoi_full;

%% Remove input metabolites and unbalanced metabolites from the 
% stoichiometric matrix

[~, input_metabs] = xlsread(file1, 'input_metab', '', 'basic');
[~, metabs_remove] = xlsread(file1, 'metab_to_remove', '', 'basic');

%% Importing and converting atom mapping 

AM_info = text(:, 3);
length_AM = length(AM_info);
length_rxn = length(rxn);

if (length_AM ~= length_rxn)
    disp(['Error: Array lengths do not match, [AM,RXN] = ',...
        num2str([length_AM, length_rxn])]);
end

%% Sort stoichiometric information to consolidate with atom-transition info
clear metab_idx
for jj = 1:length(rxn)
    metab_idx{jj, 1} = [];
end
for ii = length(metab_char):-1:1
flag_metab = metab_char{ii};
    for jj = 1:length(rxn)
        flag_rxn = rxn{jj};
        flag_idx = strfind(flag_rxn, flag_metab);
        if (~isempty(flag_idx))
            clear flag_ii
            flag_ii = repmat(ii, size(flag_idx));
            flag_pair = [flag_idx; flag_ii];
            metab_idx{jj} = [metab_idx{jj}, flag_pair];
        end
    end
end

rxn_lengths = cellfun(@length, metab_idx);
max_metab = max(rxn_lengths) + 1;

%% Creating base character array for storing atom-map information in the
% format suitable for AMMgen2 and IMMgen3
for kk = 1:max_metab
        for jj = 1:length(rxn)
            AM(kk, :, jj) = '000000000000';
        end
end

for jj = 1:length(rxn)
%     rxn{jj}
    clear sorted_flag flag_AM flag_split flag_AM_cat
    sorted_flag = (sortrows(metab_idx{jj}'));
    str_idx = num2str(sorted_flag(:, 2),'%02d');
    flag_AM = AM_info{jj};
    flag_split = strsplit(flag_AM, ',')';
    flag_AM_cat = strcat(str_idx, flag_split);
    
    for kk = 1:size(flag_AM_cat, 1)
        
        flag_strlen = length(flag_AM_cat{kk});
        AM(kk, 1:flag_strlen, jj) = char(flag_AM_cat{kk});
    end
end

AM_rev = AM(:, :, rev_rxn);
AM_irrev = AM(:, :, irrev_rxn);
AM_cell(:, :, fwd_rxn_idx) = AM_rev;
AM_cell(:, :, bkd_rxn_idx) = AM_rev;
AM_full = cat(3, AM_cell, AM_irrev);
%% Extract metabolite lengths 

clear text
[text_metab_size, text_metab_list] = xlsread(file1, 'metab_size', '', 'basic');
[~, flag_size_idx] = ismember(metab_char, text_metab_list);
metab_size = (text_metab_size(flag_size_idx'))';

