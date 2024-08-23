function [remove_idx,nterms] = ODE_fn_gen...
    (Sfull, metab_char, AM_full, metab_size, unlabeled_metabs,input_metabs,...
    balance_metabs, n_flux_param,n_c_param,...
    fwd_rxn_idx,bkd_rxn_idx,F,metab_eq,AMeq,Seq,eq_unbalance)


%% Setup amd IMM generation

% Create vector to identify unlabeled metabolites
varFlag = zeros(1,length(metab_size));
[~,unlabeled_idx] = ismember(unlabeled_metabs, metab_char); 
varFlag(unlabeled_idx) = 1;
kIDVsize = sum((varFlag).*(2.^metab_size));
mx = zeros(kIDVsize,1);
for i = unlabeled_idx'
    index = kIDVindex(i,metab_size,varFlag);
   mx(index(1),1) = 1;
end

% Create vector to identify input metabs
ipFlag = zeros(1,length(metab_size));
[~,input_idx] = ismember(input_metabs, metab_char); 
ipFlag(input_idx) = 1;

tnM = size(Sfull,1);

uIDVsize = sum((1-varFlag).*(2.^metab_size));

% Generating major matrix with all IMMs of the model
[IMMfull,IMMindex] = modelIMM_gen(metab_size,varFlag,Sfull,AM_full);

%% Symbolic variable generation

% Symbolic variable generation of model parameters and IDV vectors
flux = sym('v', [n_flux_param 1]); % Flux vector is of the form {vFwd,vBkd}
conc = sym('c',[numel(balance_metabs) 1]);

% Generate symbolic IDV vector
x = [];
for i = 1:tnM
    if (varFlag(i)==0)
        tempVar = sym([metab_char{i},'_%d'],[2^metab_size(i) 1]);
        x = [x;tempVar];
        clear tempVar
    end
end

% Function vector for isotopomer balance equations
A = zeros(uIDVsize,1);
A = sym(A);


%% Generate the expression for equilibrium metabolites as a function of other metabolites

% Generate IMMs for EQ reactions
if(isempty(metab_eq))
    vareq = [];
    eq_unbalance = [];
else
    vareq = zeros(1,numel(metab_eq));
    vareq(ismember(metab_eq,unlabeled_metabs)) = 1;
    sizeeq = metab_size(find(ismember(metab_char,metab_eq)));
    [IMMfulleq,IMMindexeq] = modelIMM_gen(sizeeq,vareq,Seq,AMeq);
end

% Loop over equilibrium metabolites and develop their expressions
for p = 1:numel(eq_unbalance)
     
        
     k = find(ismember(metab_eq,eq_unbalance(p)));
     stoimet = Seq(k,:);
     rxnind1 = find(stoimet>0);
     rxnInd = 0;     
     
     % Vector to store the expression for IDVs
     B = zeros(2^sizeeq(k),1); 
     B = sym(B);

     m = find(ismember(metab_char,eq_unbalance(p)));
     indexP = uIDVindex(m,metab_size,varFlag,0);
    
     for i = rxnind1   
    
            stoirxn = Seq(:,i)';
            atom = AMeq(:,:,i);
            rxnind2 = find(stoirxn<0);
            flag = 0;            
            [atomrow,~]=size(atom);
    
            % Finding the index of the product metabolite in the 'AMeq'
            for l = 1:atomrow
                if (k==str2double(atom(l,1:2)))
                    met2 = l;
                    break;
                end
            end
            
            % Loop to generate IMM for reactant 'j'->'k'
            while (flag==0)
            
            % Initializing the rxnIDV vector
            rxnIDV = ones(2^vareq(k),1);  
            clear IDVvec 
            %term = 1;
            for j = rxnind2
    
                stoi_j = abs(stoirxn(j));
                j_counter = 1; % Required for metabolites producing symmetric molecules
                while(stoi_j-j_counter+1>0)
                    rxnInd = rxnInd + 1;
    
                    IMM = referIMM(rxnInd,k,j,IMMfulleq,IMMindexeq,sizeeq,vareq);
    
                    % Product of flux j->k and IDVj
                    if vareq(j)==0
                        indexR = uIDVindex(find(ismember(metab_char,metab_eq(j))), metab_size,varFlag,0);
                        vIDV = x(indexR(1):indexR(2));
                    else
                        indexR = kIDVindex(find(ismember(metab_char,metab_eq(j))),metab_size,varFlag);
                        vIDV = mx(indexR(1):indexR(2));
                    end
    
                    if sum(sum(IMM))~=0
                        rxnIDV = rxnIDV.*(IMM*vIDV);
                    end
                    j_counter=j_counter+1;
                end
            end
            t1 = atom(met2,3:2+sizeeq(k));
            t2 = flip(atom(met2+1,3:2+sizeeq(k)));
            if (str2double(atom(met2+1,1:2))==k&&(strcmp(t1,t2)==0))
                met2 = met2 + 1;
            else
                flag = 1;
            end
              B = B + rxnIDV/Seq(k,i);
            end
            
     end

    % Set the expression 
    x(indexP(1):indexP(2)) = B;

end


%% IMM equation generation
 if(isempty(metab_eq))
    eq_unbalance = 'PLACEHOLDER';
 end
for k = 1:tnM
    
    if (varFlag(k)==0 && ipFlag(k)==0 && ~sum(ismember(eq_unbalance,metab_char(k)) ))
        
          % Divide by the concentration
        [~ , c_idx] = ismember(metab_char(k), balance_metabs); 
        k_conc = conc(c_idx);
        
        indexP = uIDVindex(k,metab_size,varFlag,0);

        stoimet = Sfull(k,:);
        rxnind1 = find(stoimet>0);
        rxnInd = 0;    
        
        % Loop for reactions that produce metabolite 'k'
        for i = rxnind1   
    
            stoirxn = Sfull(:,i)';
            atom = AM_full(:,:,i);
            rxnind2 = find(stoirxn<0);
            flag = 0;            
            [atomrow,~]=size(atom);
    
            % Finding the index of the product metabolite in the 'AM_full'
            for l = 1:atomrow
                if (k==str2double(atom(l,1:2)))
                    met2 = l;
                    break;
                end
            end
            
            % Loop to generate IMM for reactant 'j'->'k'
            while (flag==0)
            
            % Initializing the rxnIDV vector
            rxnIDV = ones(2^metab_size(k),1);    
            clear IDVvec dIDVvec
            term = 1;
            for j = rxnind2
    
                stoi_j = abs(stoirxn(j));
                j_counter = 1; % Required for metabolites producing symmetric molecules
                while(stoi_j-j_counter+1>0)
                    rxnInd = rxnInd + 1;
    
                    IMM = referIMM(rxnInd,k,j,IMMfull,IMMindex,metab_size,varFlag);
    
                    % Product of flux j->k and IDVj
                    if varFlag(j)==0
                        indexR = uIDVindex(j, metab_size,varFlag,0);
                        vIDV = x(indexR(1):indexR(2));
                    else
                        indexR = kIDVindex(j,metab_size,varFlag);
                        vIDV = mx(indexR(1):indexR(2));
                    end
    
                    if sum(sum(IMM))~=0
                        rxnIDV = rxnIDV.*(IMM*vIDV);
                        if varFlag(j)==0
                            IDVvec{term} = IMM*vIDV;

                            term = term+1;
                        end
                    end
                    j_counter=j_counter+1;
                end
            end
            t1 = atom(met2,3:2+metab_size(k));
            t2 = flip(atom(met2+1,3:2+metab_size(k)));
            if (str2double(atom(met2+1,1:2))==k&&(strcmp(t1,t2)==0))
                met2 = met2 + 1;
            else
                flag = 1;
            end
                
                if(sum(i == bkd_rxn_idx))
                    flux_val = (F*flux(i)/(1-flux(i)));
                elseif(sum(i == fwd_rxn_idx))
                    flux_val = (flux(i) + F*flux(i+1)/(1-flux(i+1)));
                else
                    flux_val = flux(i);
                end
                A(indexP(1):indexP(2)) = A(indexP(1):indexP(2)) + flux_val*rxnIDV/k_conc;
               
            
            end
        end
        
        % Loop for reactions that consume metabolite 'k'
        stoimet2 = Sfull(k,:);
        rxnind3 = find(stoimet2<0);
        
        for i = rxnind3
            % fprintf([num2str(i),': ',dispRxn(i,Sfull,metNAMe),'\n'])
            if(sum(i == bkd_rxn_idx))
                rxnIDV = stoimet2(i)*(F*flux(i)/(1-flux(i)))*x(indexP(1):indexP(2))/k_conc;
            elseif(sum(i == fwd_rxn_idx))
                rxnIDV = stoimet2(i)*(flux(i) + F*flux(i+1)/(1-flux(i+1)))*x(indexP(1):indexP(2))/k_conc;
            else
                rxnIDV = stoimet2(i)*flux(i)*x(indexP(1):indexP(2))/k_conc;
            end
            A(indexP(1):indexP(2)) = A(indexP(1):indexP(2)) + rxnIDV;
        end
        
        
    end
end

%% Remove the terms for equilibrium and input metabolites from the expressions

remove_idx = [];

% Remove equilibrium metabolites
if(isempty(metab_eq))
    eq_unbalance = [];
end
for i = 1:numel(eq_unbalance)

    m = find(ismember(metab_char,eq_unbalance(i)));
    indexP = uIDVindex(m,metab_size,varFlag,0);
    remove_idx = [remove_idx,[indexP(1):indexP(2)]];

end
x(remove_idx) = [];
A(remove_idx) = [];

% Remove input metabolites
remove_idx = [];
varFlag2 = varFlag;
if(~isempty(metab_eq))
    varFlag2(ismember(metab_char,eq_unbalance)) = 1;
end
labeled_inputs = setdiff(input_metabs, unlabeled_metabs);
for i = 1:numel(labeled_inputs)

    m = find(ismember(metab_char,labeled_inputs(i)));
    indexP = uIDVindex(m,metab_size,varFlag2,0);
    remove_idx = [remove_idx,[indexP(1):indexP(2)]];

end
Afull = A;
A(remove_idx) = [];
nterms = numel(A);


%% Generate M-file for IMM equations

filepath = 'serine_ODE.m';
matlabFunction(A,'file',filepath,'vars',{x,flux,conc});
disp('Successfully created serine gradient file.');


% Store IMM equations in text file
filepath_text = 'serine_ODE.txt';
fid = fopen(filepath_text,'w');
for k = 1:tnM
    if(varFlag2(k)==0)
        fprintf(fid,'%s\n',metab_char{k});
        A_idx = uIDVindex(k,metab_size,varFlag2,0);
        for fw = A_idx(1):A_idx(2)
            fw_flag = char(Afull(fw));
            fprintf(fid,'%s\n',fw_flag);
        end
    end
end
fclose(fid);
end