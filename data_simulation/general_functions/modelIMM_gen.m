function [IMMfull,IMMindex] = modelIMM_gen(metSize,varFlag,stoimatrix,atommatrix)
% Function to store all combinations of product-reactant IMMs
% Each metabolite's IMM for all the reactions it is involved in 
% Reactions will be listed column-wise
% Pair-IMM for each reaction will be stored column-wise

[nM,~] = size(stoimatrix);
% nM is the number of mass-balanced metabolites

B = zeros(nM,1);
for k = 1:nM
    stoimet = stoimatrix(k,:);
    rxnind1 = find(stoimet>0);
    B(k,1) = length(rxnind1);
end

maxRxn = max(B);
A = zeros(nM,maxRxn);
IMMindex = zeros(nM,maxRxn);
indexA = [0 0];
rowInd = 0;
% Apply isotope balance for metabolite k
for k = 1:nM

    if (varFlag(k)==0)
        % track first and last indices of isotopers of a metabolite in a
        % vector of all metabolite isotopomers
        indexA(1) = 1 + rowInd;
        indexA(2) = indexA(2) + 2^metSize(k);
        rowInd = rowInd + 2^metSize(k);
        % stoichiometry vector of metabolite k in all reactions
        stoimet = stoimatrix(k,:);
        % in which reaction metabolite k is produced
        rxnind1 = find(stoimet>0);
        colIndA = 0;
        colIndC = 0;
        colIndIMM = 0;
        indexB = [0 0];
        % Iterate over the reactions that produce metabolite k
        for i = rxnind1
           
            colIndA = colIndA + 1;
            % stoichiometry vector of metabolites in reaction i in which
            % metabolite k is produced
            stoirxn = stoimatrix(:,i)';
            
            atom = atommatrix(:,:,i); % AMM of reaction i
            % whcih metabolite is consumed in reaction i
            rxnind2 = find(stoirxn<0);
            flag = 0;            
            [atomrow,~] = size(atom);
            tempIMMsize = 0;
            % find metabolite k (product) in AMM
            for l = 1:atomrow
                if (k==str2double(atom(l,1:2)))
                    met2 = l; % row number of metabolite k in AMM of reaction i
                    break;
                end
            end
            while (flag == 0)
                % Iterate over the meatobiltes j that produce metabolite k
                % through reaction i
                for j = rxnind2

                    stoi_j = abs(stoirxn(j)); % stoi coeff of metabolite j (reactant of reaction i)
                    j_counter = 1; % counts number of reactants in reaction i
                    while(stoi_j-j_counter+1>0)
                        % track first and last indices of isotopers of
                        % metabolite j (reactants)
                        if (colIndIMM==0)
                            indexB(1) = indexB(1)+1;
                        else
                            indexB(1) = 1 + colIndIMM;
                        end
                        indexB(2) = indexB(1) + 2^metSize(j)-1; 
                        tempIMMsize = tempIMMsize + 2^metSize(j);
                        met1 = stoi_j-j_counter;
                        
                        IMM = IMMgen(metSize,i,stoimatrix,atommatrix,j,k,met1,met2);
                        IMMfull(indexA(1):indexA(2),indexB(1):indexB(2)) = IMM;
                        colIndC = colIndC + 1;
                        colIndIMM = colIndIMM + 2^metSize(j);
                        IMMindex(k,colIndC) = 2^metSize(j);
                        j_counter = j_counter+1;
                    end
                end

                A(k,colIndA) = length(rxnind2);
                t1 = atom(met2,3:2+metSize(k)); % product k carbons
                t2 = flip(atom(met2+1,3:2+metSize(k)));
                if (str2double(atom(met2+1,1:2))==k&&(strcmp(t1,t2)==0))
                    colIndA = colIndA + 1;
                    met2=met2+1;
                else
                    flag = 1;
                end
               
            end
        end
    end 
end

fprintf('\nIMMs generated\n')
end