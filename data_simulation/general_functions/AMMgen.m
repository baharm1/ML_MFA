function AMM = AMMgen(met,stoi,atom,label1,label2,flag1,flag2)
% This program generates the Atomic Mapping Matrix for a single reactant-
% product pair of metabolites met1(r) and met2(p)
% Inputs: Metabolite vector, stoichiometric vector for reaction, atom 
% transition information, labels of the reactant and product corresponding 
% to their index 
% This function is written based on the AMM method by Zupke et al (1994)

[m,~] = size(atom);

% Convert global index of metabolite to index corresponding to atom
% transition matrix
for i = 1:m % Reactant   
    if (label1==str2double(atom(i,1:2)))
        met1 = i;
        break;
    end
end
for i = 1:m % Product
    if (label2==str2double(atom(i,1:2)))
        met2 = i;
        break;
    end
end

met1 = met1 + flag1;

if (flag2-met2>0)
    met2 = flag2;
end

length1 = met(label1);
length2 = met(label2);

% Initialize AMM
AMM = zeros(length2,length1);
% Generate AMM
for r = 3:(length1+2)
    for p = 3:(length2+2)
        
            if (atom(met1,r)==atom(met2,p))
                AMM(p-2,r-2)=1;
            end

   end
end
end

function IMM = referIMM(rxnInd,prodInd,reacInd,IMMfull,IMMindex,metSize,varFlag)
% Extract the IMM for corresponding reactant-product pair of the reaction
% from the parent matrix containing all IMMs for the model

indexA = [1 0];
if (prodInd~=1)
    indexA(1) = indexA(1)+sum((1-varFlag(1:prodInd-1)).*(2.^metSize(1:prodInd-1)));
end

indexA(2) = sum((1-varFlag(1:prodInd)).*(2.^metSize(1:prodInd)));

indexB = [1 1];

if (rxnInd~=1)
    indexB(1) = indexB(1)+sum(IMMindex(prodInd,1:rxnInd-1));
end
indexB(2) = indexB(1) + IMMindex(prodInd,rxnInd)-1;

% indexA
% indexB
IMM = IMMfull(indexA(1):indexA(2),indexB(1):indexB(2));

end