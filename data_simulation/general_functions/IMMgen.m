function IMM = IMMgen(met,rxn,stoimatrix,atommatrix,metR,metP,met1,met2)
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
% This function is a subfunction for the IMM method based on 
% Schmidt et al. (1997)

stoi = stoimatrix(:,rxn)';
atom = atommatrix(:,:,rxn);

length1 = met(metR); % Size of reactant
length2 = met(metP); % Size of product

% Vectors containing all possible labeling patterns of reactant(IDV1)
% and product(IDV2)
IDV1 = zeros(2^length1,1);

for i = 0:(2^length1-1)
    bin = dec2bin(i,length1);
    for j = 1:length1
        IDV1(i+1,j) = str2double(bin(j));
    end
end

IMM = zeros(2^length2,2^length1);
AMM = AMMgen(met,stoi,atom,metR,metP,met1,met2);
m = ones(1,length2);


% Vector identifying the atoms in the product molecule which do not come
% from chosen reactant molecule
for i = 1:length2
    if (sum(AMM(i,:)~=0))
        m(i) = 0;
    end
end

% Generating indices of all possible product IDVs which can come from the
% labeled/unlabeled reactant
extraP = sum(m);
k = zeros(extraP+1,1);
k(1) = 1;


if (extraP>0)
nExtraP = 2^extraP-1;
extraIDV = zeros(nExtraP,extraP);

for i = 1:nExtraP
    bin = dec2bin(i,extraP);
    for j = 1:extraP
        extraIDV(i,j) = str2double(bin(j));
    end
end

extraPIDV = zeros(nExtraP,length2);

for i = 1:nExtraP
    u = 1;
    for j = 1:length2
        if (m(j)==1)
            extraPIDV(i,j) = extraIDV(i,u);
            u = u+1;
        end
    end
end

for i = 1:nExtraP
    extraPIDVbin = num2str(extraPIDV(i,:));
    k(i+1) = bin2dec(extraPIDVbin)+1;
end

end

if (sum(sum(AMM))==0)
    IMM = zeros(2^length2,2^length1);
else
        
        for i = 1:(2^length1)
           
            pIDV(:,i) = AMM*(IDV1(i,:))';
            pIDVbin(i,:)=num2str(pIDV(:,i));

        end

        for l = 1:length(k)
        for i = 1:(2^length1)        
           
            index = bin2dec(pIDVbin(i,:)) + k(l);
            IMM(index,i) = 1;
            
        end
        end

end


t1 = atom(met2,3:2+met(metP));
t2 = flip(atom(met2+1,3:2+met(metP)));
if ((str2double(atom(met2+1,1:2))==metP)&&(strcmp(t1,t2)==1))
    met2 = met2 + 1;
    AMM2 = AMMgen(met,stoi,atom,metR,metP,met1,met2);

    if (sum(sum(AMM2))==0)
    IMM2 = zeros(2^length2,2^length1);

    elseif(length1 < length2 && (abs(length1-length2) == 1) || length1 == 1) % Special case 
        
        extraPIDV2 = extraPIDV;
        extraPIDV2(:,1) = extraPIDV(:,end);
        extraPIDV2(:,end) = extraPIDV(:,1);
        k2 = zeros(numel(k),1);
        k2(1) = 1;
        for i = 1:nExtraP
            extraPIDVbin = num2str(extraPIDV2(i,:));
            k2(i+1) = bin2dec(extraPIDVbin)+1;
        end

        
        for i = 1:(2^length1)
           
            pIDV(:,i) = AMM2*(IDV1(i,:))';
            pIDVbin(i,:)=num2str(pIDV(:,i));

        end
        
        for i = 1:(2^length1)        
        for l = 1:length(k2)   
            index = bin2dec(pIDVbin(i,:)) + k2(l);
            IMM2(index,i) = 1;
            
        end
        end

    else
        
        for i = 1:(2^length1)
           
            pIDV(:,i) = AMM2*(IDV1(i,:))';
            pIDVbin(i,:)=num2str(pIDV(:,i));

        end
        
        for l = 1:length(k)
        for i = 1:(2^length1)        
           
            index = bin2dec(pIDVbin(i,:)) + k(l);
            IMM2(index,i) = 1;
            
        end
        end

    end

    IMM = (IMM + IMM2)/2;
end
end
