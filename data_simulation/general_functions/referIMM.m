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