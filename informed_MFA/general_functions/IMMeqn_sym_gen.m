function [num_nlcons] = IMMeqn_sym_gen(Sfull, metab_char, AM_full, ...
    metab_size, nflux, varFlag, mx, input_metabs, metabs_remove)
global timestamp_gl filepath_stored_gl
% Total number of metabolites (tnM) and reactions (tnR) from stoichiometric table S
tnM = size(Sfull, 1);
balance_metabs = setdiff(metab_char, [input_metabs; metabs_remove]');

%% Generating major matrix with all IMMs of the model

[IMMfull, IMMindex] = modelIMM_gen(metab_size, varFlag, Sfull, AM_full);

%% Symbolic variable generation of flux and IDV vectors
flux = sym('v', [nflux 1]); % Flux vector is of the form {vFwd,vBkd}

x = [];
for i = 1:tnM
    if (varFlag(i) == 0)
        tempVar = sym(metab_char{i}, [2 ^ metab_size(i) 1]);
        x = [x; tempVar];
        clear tempVar
    end
end

%% IMM equation generation
disp('Creating f_NL_const.m ...')
uIDVsize = sum((1 - varFlag) .* (2 .^ metab_size));

% Function vector for isotopomer balance equations
A = zeros(uIDVsize, 1);
A = sym(A);

for k = 1:tnM

if (varFlag(k) == 0)
    indexP = uIDVindex(k, metab_size, varFlag, 0);
    
    stoimet = Sfull(k, :);
    rxnind1 = find(stoimet > 0);
    rxnInd = 0;            
    
    % Loop for reactions that produce metabolite 'k'
    for i = rxnind1   

        stoirxn = Sfull(:, i)';
        atom = AM_full(:, :, i);
        rxnind2 = find(stoirxn < 0);
        flag = 0;            
        [atomrow, ~]=size(atom);

        % Finding the index of the product metabolite in the 'AM_full'
        for l = 1:atomrow
            if (k == str2double(atom(l, 1:2)))
                met2 = l;
                break;
            end
        end
        
        % Loop to generate IMM for reactant 'j'->'k'
        while (flag == 0)
        
        % Initializing the rxnIDV vector
        rxnIDV = ones(2 ^ metab_size(k), 1);    
        for j = rxnind2

            stoi_j = abs(stoirxn(j));
            j_counter = 1; % Required for metabolites producing symmetric molecules
            while(stoi_j - j_counter + 1 > 0)
                rxnInd = rxnInd + 1;

                IMM = referIMM(rxnInd, k, j, IMMfull, IMMindex, metab_size, ...
                    varFlag);

                % Product of flux j->k and IDVj
                if varFlag(j) == 0
                    indexR = uIDVindex(j, metab_size, varFlag, 0);
                    vIDV = x(indexR(1):indexR(2));
                else
                    indexR = kIDVindex(j, metab_size, varFlag);
                    vIDV = mx(indexR(1):indexR(2));
                end

                if sum(sum(IMM)) ~= 0
                    rxnIDV = rxnIDV .* (IMM * vIDV);
                end
                j_counter = j_counter + 1;
            end
        end
        t1 = atom(met2, 3:2 + metab_size(k));
        t2 = flip(atom(met2 + 1,3:2 + metab_size(k)));
        if (str2double(atom(met2 + 1, 1:2)) == k && (strcmp(t1, t2) == 0))
            met2 = met2 + 1;
        else
            flag = 1;
        end
        A(indexP(1):indexP(2)) = A(indexP(1):indexP(2)) + flux(i) * rxnIDV;
        end
    end
    
    % Loop for reactions that consume metabolite 'k'
    stoimet2 = Sfull(k,:);
    rxnind3 = find(stoimet2 < 0);
    indexPc = uIDVindex(k, metab_size, varFlag, 0);
    
    for i = rxnind3
        % fprintf([num2str(i),': ',dispRxn(i,Sfull,metNAMe),'\n'])
        rxnIDV = stoimet2(i) * flux(i) * x(indexPc(1):indexPc(2));
        A(indexP(1):indexP(2)) = A(indexP(1):indexP(2)) + rxnIDV;
    end
    
end
end

%% Only keep the metabolites that are to be balanced

for i = 1:numel(metabs_remove)
    m = find(ismember(metab_char, metabs_remove(i)));
    index = uIDVindex(m, metab_size, varFlag, 0);
    if(index(1) > 0)
        A(index(1):index(2)) = 0;
    end
end

for i = 1:numel(input_metabs)
    m = find(ismember(metab_char, input_metabs(i)));
    index = uIDVindex(m, metab_size, varFlag, 0);
    if(index(1) > 0)
        A(index(1):index(2)) = 0;
    end
end
A_old = A;
A(A == 0) = [];
num_nlcons = numel(A);

%% Generate M-file for IMM equations
filepath = 'f_NL_const.m';

matlabFunction(A, 'file', filepath, 'vars', {flux, x});
copyfile(filepath, filepath_stored_gl)
disp([filepath, ' successfully created and copied to ', filepath_stored_gl])

% Store IMM equations in text file
filepath_text = 'f_NL_const.txt';
fid = fopen(filepath_text, 'w');
for k = 1:tnM
    if(varFlag(k) == 0 && ismember(metab_char(k), balance_metabs))
        fprintf(fid, '%s\n', metab_char{k});
        A_idx = uIDVindex(k, metab_size, varFlag, 0);
        for fw = A_idx(1):A_idx(2)
            fw_flag = char(A_old(fw));
            fprintf(fid, '%s\n', fw_flag);
        end
    end
end
fclose(fid);
copyfile(filepath_text, filepath_stored_gl)

end

%% Subroutines specific to parent function
function [IMMfull, IMMindex] = modelIMM_gen(metSize, varFlag, stoimatrix, ...
    atommatrix)
% Function to store all combinations of product-reactant IMMs
% Each metabolite's IMM for all the reactions it is involved in 
% Reactions will be listed column-wise
% Pair-IMM for each reaction will be stored column-wise

[nM,~] = size(stoimatrix);
% nM is the number of mass-balanced metabolites

B = zeros(nM, 1);
for k = 1:nM
    stoimet = stoimatrix(k, :);
    rxnind1 = find(stoimet > 0);
    B(k, 1) = length(rxnind1);
end

maxRxn = max(B);
A = zeros(nM, maxRxn);
IMMindex = zeros(nM, maxRxn);
indexA = [0 0];
rowInd = 0;
% Apply isotope balance for metabolite k
for k = 1:nM

    if (varFlag(k) == 0)
    
        indexA(1) = 1 + rowInd;
        indexA(2) = indexA(2) + 2 ^ metSize(k);
        rowInd = rowInd + 2 ^ metSize(k);
        stoimet = stoimatrix(k, :);
        rxnind1 = find(stoimet > 0);
        colIndA = 0;
        colIndC = 0;
        colIndIMM = 0;
        indexB = [0 0];
        % Iterate over the reactions that produce metabolite k
        for i = rxnind1

            colIndA = colIndA + 1;
            stoirxn = stoimatrix(:, i)';
            atom = atommatrix(:, :, i);
            rxnind2 = find(stoirxn < 0);
            flag = 0;            
            [atomrow, ~]=size(atom);
            tempIMMsize = 0;
            for l = 1:atomrow
                if (k == str2double(atom(l, 1:2)))
                    met2 = l;
                    break;
                end
            end
            while (flag == 0)
                % Iterate over the meatobiltes j that produce metabolite k
                % through reaction i
                for j = rxnind2

                    stoi_j = abs(stoirxn(j));
                    j_counter = 1;
                    while(stoi_j - j_counter + 1 > 0)
                        if (colIndIMM == 0)
                            indexB(1) = indexB(1) + 1;
                        else
                            indexB(1) = 1 + colIndIMM;
                        end

                        indexB(2) = indexB(1) + 2 ^ metSize(j) - 1; 
                        tempIMMsize = tempIMMsize + 2 ^ metSize(j);
                        met1 = stoi_j - j_counter;
                        
                        IMM = IMMgen(metSize, i, stoimatrix, atommatrix, ...
                            j, k, met1, met2);
                        IMMfull(indexA(1):indexA(2), indexB(1):indexB(2)) = IMM;
                        colIndC = colIndC + 1;
                        colIndIMM = colIndIMM + 2 ^ metSize(j);
                        IMMindex(k, colIndC) = 2 ^ metSize(j);
                        j_counter = j_counter + 1;
                    end
                end

                A(k, colIndA) = length(rxnind2);
                t1 = atom(met2, 3:2 + metSize(k));
                t2 = flip(atom(met2 + 1,3:2 + metSize(k)));
                if (str2double(atom(met2 + 1, 1:2)) == k && (strcmp(t1, t2) == 0))
                    colIndA = colIndA + 1;
                    met2 = met2 + 1;
                else
                    flag = 1;
                end
               
            end
        end
    end 
end

fprintf('\nIMMs generated\n')
end

function IMM = IMMgen(met, rxn, stoimatrix, atommatrix, metR, metP, met1, met2)
% This function generates the Isotopomer Matrix Map for 
% given reaction-product pair of metabolites met1 and met2
% Inputs: 
% 1. Matrix 'met' contains the # of labeled atoms corresponding to the
% metabolite index
% 2. Reaction label index
% 3. Stoichiometric matrix for reaction
% 4. atom transition information matrix
% 5,6. Labels of the reactant and product corresponding to their index
% 7. Flag variable to account for multiple instances of same product in the
% reaction
% This function, written by Abhinav Achreja is a subfunction for the IMM
% method based on Schmidt et al (1997)
% Version 1.2, modified on 03-31-2013
% Update: Fixed issues with IMM of condensation reactions with reactants
% contributing less number of atoms than present in the product
% Update: 10/30/2013: Added proper generation of IMMs for symmetric product
% metabolites

stoi = stoimatrix(:, rxn)';
atom = atommatrix(:, :, rxn);

length1 = met(metR);
length2 = met(metP);

% Vectors containing all possible labeling patterns of reactant(IDV1)
% and product(IDV2)
IDV1 = zeros(2 ^ length1, 1);

for i = 0:(2 ^ length1 - 1)
    bin = dec2bin(i, length1);
    for j = 1:length1
        IDV1(i + 1, j) = str2double(bin(j));
    end
end

IMM = zeros(2 ^ length2, 2 ^ length1);
AMM = AMMgen(met, stoi, atom, metR, metP, met1, met2);
m = ones(1, length2);

% Vector identifying the atoms in the product molecule which do not come
% from chosen reactant molecule
for i = 1:length2
    if (sum(AMM(i, :) ~= 0))
        m(i) = 0;
    end
end

% Generating indices of all possible product IDVs which can come from the
% labeled/unlabeled reactant

extraP = sum(m);
k = zeros(extraP + 1,1);
k(1) = 1;

if (extraP > 0)
nExtraP = 2 ^ extraP - 1;
extraIDV = zeros(nExtraP, extraP);

for i = 1:nExtraP
    bin = dec2bin(i, extraP);
    for j = 1:extraP
        extraIDV(i, j) = str2double(bin(j));
    end
end

extraPIDV = zeros(nExtraP, length2);

for i = 1:nExtraP
    u = 1;
    for j = 1:length2
        if (m(j) == 1)
            extraPIDV(i, j) = extraIDV(i, u);
            u = u + 1;
        end
    end
end


for i = 1:nExtraP
    extraPIDVbin = num2str(extraPIDV(i, :));
    k(i + 1) = bin2dec(extraPIDVbin) + 1;
end

end

if (sum(sum(AMM)) == 0)
    IMM = zeros(2 ^ length2, 2 ^ length1);
else
        
        for i = 1:(2 ^ length1)
           
            pIDV(:, i) = AMM * (IDV1(i, :))';
            pIDVbin(i, :) = num2str(pIDV(:, i));

        end
    
        for l = 1:length(k)
        for i = 1:(2 ^ length1)        
           
            index = bin2dec(pIDVbin(i, :)) + k(l);
            IMM(index, i) = 1;
            
        end
        end

end

t1 = atom(met2, 3:2 + met(metP));
t2 = flip(atom(met2 + 1, 3:2 + met(metP)));
if ((str2double(atom(met2 + 1, 1:2)) == metP) && (strcmp(t1, t2) == 1))
    met2 = met2 + 1;
    AMM2 = AMMgen(met, stoi, atom, metR, metP, met1, met2);
    
    if (sum(sum(AMM2)) == 0)
    IMM2 = zeros(2 ^ length2, 2 ^ length1);
    else
        
        for i = 1:(2 ^ length1)
           
            pIDV(:, i) = AMM2 * (IDV1(i, :))';
            pIDVbin(i, :)=num2str(pIDV(:, i));

        end
    
        for l = 1:length(k)
        for i = 1:(2 ^ length1)        
           
            index = bin2dec(pIDVbin(i, :)) + k(l);
            IMM2(index, i) = 1;
            
        end
        end

    end
    IMM = (IMM + IMM2) / 2;
end
end

function AMM = AMMgen(met, stoi, atom, label1, label2, flag1, flag2)
% This program generates the Atomic Mapping Matrix for a single reactant-
% product pair of metabolites met1(r) and met2(p)
% Inputs: Metabolite vector, stoichiometric vector for reaction, atom 
% transition information, labels of the reactant and product corresponding 
% to their index 
% This program, written by Abhinav Achreja based on the AMM method 
% by Zupke et al (1994)
% Version 1.1, modified on 03-21-2013

[m,~] = size(atom);

% Convert global index of metabolite to index corresponding to atom
% transition matrix
for i = 1:m % Reactant   
    if (label1 == str2double(atom(i, 1:2)))
        met1 = i;
        break;
    end
end
for i = 1:m % Product
    if (label2 == str2double(atom(i, 1:2)))
        met2 = i;
        break;
    end
end

met1 = met1 + flag1;

if (flag2 - met2 > 0)
    met2 = flag2;
end

length1 = met(label1);
length2 = met(label2);

% Initialize AMM
AMM = zeros(length2, length1);
% Generate AMM
for r = 3:(length1 + 2)
    for p = 3:(length2 + 2)
        
            if (atom(met1, r) == atom(met2, p))
                AMM(p - 2, r - 2) = 1;
            end

   end
end
end

function IMM = referIMM(rxnInd, prodInd, reacInd, IMMfull, IMMindex, metSize, varFlag)
% Extract the IMM for corresponding reactant-product pair of the reaction
% from the parent matrix containing all IMMs for the model
% Version 1.2 created by Abhinav Achreja (10-2013)

indexA = [1 0];
if (prodInd ~= 1)
    indexA(1) = indexA(1) + sum((1 - varFlag(1:prodInd - 1)) .* (2 .^ metSize(1:prodInd - 1)));
end

indexA(2) = sum((1 - varFlag(1:prodInd)) .* (2 .^ metSize(1:prodInd)));

indexB = [1 1];

if (rxnInd ~= 1)
    indexB(1) = indexB(1) + sum(IMMindex(prodInd, 1:rxnInd - 1));
end
indexB(2) = indexB(1) + IMMindex(prodInd, rxnInd) - 1;

% indexA
% indexB
IMM = IMMfull(indexA(1):indexA(2), indexB(1):indexB(2));

end