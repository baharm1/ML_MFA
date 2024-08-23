function [rxn_full, input_metabs, metab_char, metab_sym, metab_size,...
    fwd_rxn_idx,bkd_rxn_idx, stoi_full,AM_full,slope_metab,...
    eq_unbalance,AMeq,Seq,metab_eq,metab_to_remove,pool_metabs,...
    unlabeled_metabs,balance_metabs, n_flux_param, n_c_param, Aeq,M0,varFlag] = import_stoich_AM(file1)

%% Import the data
[~,rxns] = xlsread(file1,'rxns','','basic');
[metab_size,text_metab_list] = xlsread(file1,'metab_size','','basic');
[~,input_metabs] = xlsread(file1,'input_metab','','basic');
[~,unlabeled_metabs] = xlsread(file1,'unlabeled','','basic');
[~,eq_unbalance] = xlsread(file1,'EQ','','basic');
[~,metab_to_remove] = xlsread(file1,'remove','','basic');
try
    [~,pool_metabs] = xlsread(file1,'pool','','basic');
catch
    pool_metabs = [];
end

%% Extract the stoitiometric matrix and reaction information

% Remove equilibrium reactions from the reaction list
eq_rxns = rxns(cat(1,rxns{:,1}) == 'Q',:);
rxns =  rxns(cat(1,rxns{:,1}) ~= 'Q',:);


% Extract reaction type information
rxn_flag = cat(1,rxns{:,1});
irrev_rxn = find(rxn_flag=='I');
exc_rxn = find(rxn_flag=='E');
rev_rxn = find(rxn_flag=='R');
rxn = rxns(:,2);

% Convert reaction equations to symbolic form
n_rxn = size(rxn,1);
for i = 1:n_rxn
    text2(i) = str2sym(rxn(i));
end

% Extract metabolite names from equations
metab_sym = symvar(text2);
metab_char = cellstr(string(metab_sym)); % Convert symbolic array of variables to strings

% Convert symbolic equations to stoichiometric matrix
nRev2 = length(rev_rxn)*2;
fwd_rxn_idx = [1:2:nRev2-1];
bkd_rxn_idx = [2:2:nRev2];
flag_stoi_matrix = equationsToMatrix(text2);
flag_stoi_matrix = double(-flag_stoi_matrix');
stoi_rev = flag_stoi_matrix(:,rev_rxn);
stoi_irrev = flag_stoi_matrix(:,irrev_rxn);
stoi_exc = flag_stoi_matrix(:,exc_rxn);
stoi_cell(:,fwd_rxn_idx) = stoi_rev;
stoi_cell(:,bkd_rxn_idx) = -stoi_rev;
stoi_full = [stoi_cell,stoi_irrev,stoi_exc];


% Save the reactions in order of the stoichiometric matrix
rxn_rev = rxn(rev_rxn);
rxn_irrev = rxn(irrev_rxn);
rxn_exc = rxn(exc_rxn);
rxn_cell(fwd_rxn_idx,1) = rxn_rev;
rxn_cell(bkd_rxn_idx,1) = rxn_rev;
rxn_full = [rxn_cell;rxn_irrev;rxn_exc];


%% Extract stoichiometric and reaction information for equilibrium reactions

if(~isempty(eq_rxns))

    % Convert reaction equations to symbolic form
    n_rxn = size(eq_rxns,1);
    clear text2
    for i = 1:n_rxn
        text2(i) = str2sym(eq_rxns(i,2));
    end
    
    % Extract metabolite names from equations
    metab_eq = cellstr(string(symvar(text2))); 
    
    % Convert symbolic equations to stoichiometric matrix
    Seq = double(-equationsToMatrix(text2)');
else
    Seq = [];
    metab_eq = [];
end

%% Extract metabolite sizes

text_metab_list = text_metab_list';
[~,flag_size_idx] = ismember(metab_char,text_metab_list);
metab_size = (metab_size(flag_size_idx'))';

%% Rearrange input metabolite names and create slope_metab
[~,input_idx] = ismember(input_metabs,metab_char);
input_idx = sort(input_idx);
input_metabs = metab_char(input_idx);
slope_metab = [];
for k = input_idx'
    if(~ismember(metab_char(k),unlabeled_metabs))
        slope_metab = [slope_metab; repmat(string(metab_char(k)),2^metab_size(k),1)];
    end
end
input_metabs = input_metabs';

%% Import atom mapping information

AM_info = rxns(:,3);

% Sort stoichiometric information to consolidate with atom-transition info
metab_idx = cell(length(rxn),1);
for ii = length(metab_char):-1:1
flag_metab = metab_char{ii};
    for jj = 1:length(rxn)
        flag_rxn = rxn{jj};
        flag_idx = strfind(flag_rxn,flag_metab);
        if (~isempty(flag_idx))
            clear flag_ii
            flag_ii = repmat(ii,size(flag_idx));
            flag_pair = [flag_idx;flag_ii];
            metab_idx{jj} = [metab_idx{jj},flag_pair];
        end
    end
end

rxn_lengths = cellfun(@length,metab_idx);
max_metab = max(rxn_lengths)+1;

% Creating base character array for storing atom-map information in the
% format suitable for AMMgen2 and IMMgen3
for kk = 1:max_metab
        for jj = 1:length(rxn)
            AM(kk,:,jj) = '000000000000';
        end
end

for jj = 1:length(rxn)
    
    clear sorted_flag flag_AM flag_split flag_AM_cat
    sorted_flag = (sortrows(metab_idx{jj}'));
    str_idx = num2str(sorted_flag(:,2),'%02d');
    flag_AM = AM_info{jj};
    flag_split = strsplit(flag_AM,',')';
    flag_AM_cat = strcat(str_idx,flag_split);
    
    for kk = 1:size(flag_AM_cat,1)
        
        flag_strlen = length(flag_AM_cat{kk});
        AM(kk,1:flag_strlen,jj) = char(flag_AM_cat{kk});
    end
end

AM_rev = AM(:,:,rev_rxn);
AM_irrev = AM(:,:,irrev_rxn);
AM_cell(:,:,fwd_rxn_idx) = AM_rev;
AM_cell(:,:,bkd_rxn_idx) = AM_rev;
AM_full = cat(3,AM_cell,AM_irrev);

%% Atom mapping information for EQ reactions

if(isempty(metab_eq))
    AMeq = [];

else

clear AM_info flag_rxn flag_idx rxns rxn_lengths max_metab AM
AM_info = eq_rxns(:,3);
rxns = eq_rxns(:,2);

% Sort stoichiometric information to consolidate with atom-transition info
clear metab_idx
metab_idx = cell(length(rxns),1);
for ii = length(metab_eq):-1:1
flag_metab = metab_eq{ii};
    for jj = 1:length(rxns)
        flag_rxn = rxns{jj};
        flag_idx = strfind(flag_rxn,flag_metab);
        if (~isempty(flag_idx))
            clear flag_ii
            flag_ii = repmat(ii,size(flag_idx));
            flag_pair = [flag_idx;flag_ii];
            metab_idx{jj} = [metab_idx{jj},flag_pair];
        end
    end
end

rxn_lengths = cellfun(@length,metab_idx);
max_metab = max(rxn_lengths)+1;

% Creating base character array for storing atom-map information in the
% format suitable for AMMgen2 and IMMgen3
for kk = 1:max_metab
        for jj = 1:length(rxns)
            AM(kk,:,jj) = '000000000000';
        end
end

for jj = 1:length(rxns)
    
    clear sorted_flag flag_AM flag_split flag_AM_cat
    sorted_flag = (sortrows(metab_idx{jj}'));
    str_idx = num2str(sorted_flag(:,2),'%02d');
    flag_AM = AM_info{jj};
    flag_split = strsplit(flag_AM,',')';
    flag_AM_cat = strcat(str_idx,flag_split);
    
    for kk = 1:size(flag_AM_cat,1)
        
        flag_strlen = length(flag_AM_cat{kk});
        AM(kk,1:flag_strlen,jj) = char(flag_AM_cat{kk});
    end
end
AMeq = AM;
end
%% Create stoichiometric matrix
Sfull = stoi_full;
Stemp = Sfull;
Stemp(:,bkd_rxn_idx) = zeros(size(Sfull,1),numel(bkd_rxn_idx)); % Set the coeff for backward exchange reactions to zero

% Remove metabolites that are not mass balanced in the model
metab_to_remove = [input_metabs;metab_to_remove;eq_unbalance];
[~,metab_remove_idx] = ismember(metab_to_remove,metab_char);
metab_remove_idx = metab_remove_idx(metab_remove_idx>0);
Stemp(metab_remove_idx,:)=[];
balance_metabs = metab_char(~ismember(metab_char,metab_to_remove));

% Save the number of parameters
Aeq = Stemp;
n_flux_param = size(Stemp,2);
n_c_param = size(Aeq,1);

%% Create M0

varFlag = zeros(1,length(metab_size));
[~,unlabeled_idx] = ismember(unlabeled_metabs, metab_char); 
varFlag(unlabeled_idx) = 1;
varFlag(ismember(metab_char,eq_unbalance)) = 1;
M0 = zeros(sum((1-varFlag).*2.^metab_size),1);
for m = 1:numel(metab_size)
    index = uIDVindex(m,metab_size,varFlag,0);
    if(index(1) > 0)
        if(ismember(metab_char(m),input_metabs))
         M0(index(1):index(2)) = -1; 
        elseif(index(1) > 0)
         M0(index(1)) = 1;
        end
    end
end
M0(M0 == -1) = [];

