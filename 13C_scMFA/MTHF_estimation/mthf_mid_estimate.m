

%% Read the input data
input_data_file = 'patient_ser_gly.xlsx';
min_sd = sqrt(10^(-5));
[MID,SD,time_points,data_metab_list,data_i_list,~] = import_expt_data_mthf(input_data_file,min_sd);

%% Calculate the MID and MC distribution of Me-THF from serine and glycine MID

for t = 2:length(time_points)
    
    ser_MID = MID(ismember(data_metab_list,'SERINE'),t);
    gly_MID = MID(ismember(data_metab_list,'GLYCINE'),t);
    ser_sd =  SD(ismember(data_metab_list,'SERINE'),t);
    gly_sd = SD(ismember(data_metab_list,'GLYCINE'),t);
    
    % Find the optimum solution vector
    [ cthf_mid(:,t), x_opt] = calc_methf_mid(ser_MID, gly_MID, ser_sd, gly_sd,[]);

    % Create multiple data points for MC simulations
    nmc = 200;
    sim_mid = mc_sampler2([ser_MID;gly_MID],[ser_sd;gly_sd],["SERINE";"SERINE";"SERINE";"SERINE";"GLY";"GLY";"GLY"],nmc);   

    % Perform the optimization for the simulated MC data
    for n = 1:nmc
      [cthf_mc(:,n),~] = calc_methf_mid(sim_mid(1:4,n), sim_mid(5:end,n), ...
          ser_sd, gly_sd,x_opt);
    end
    
    % Calaculate the mean and SD of the simulated distrubution
    cthf_mid(:,t) = mean(cthf_mc,2);
    cthf_sd(:,t) = std(cthf_mc,0,2);

end

% Add the data for time = 0 (unlabeled case) 
 cthf_mid(:,1) = [1;0];
 cthf_sd(:,1) = [0;0];

 % Save the data as percentage MID
 cthf_mid = 100*cthf_mid;
 cthf_sd = 100*cthf_sd;
 results = table(["M0";"M1"],cthf_mid,cthf_sd);
 writetable(results ,[input_data_file,'_cthf_calc.xlsx']);


function [ cthf_mid, x_opt] = calc_methf_mid(ser_MID, gly_MID, ser_sd, gly_sd,x0_opt)
   
nvar = 6 + 3;

% define the lower and upper bounds
lb = zeros(nvar,1);
ub = ones(nvar,1);

% define the objective function
obj_fn = @(x)methf_objective...
    (x,ser_MID, gly_MID, ser_sd, gly_sd); 

% define the linear contraint matrix
b = zeros(3,1);
A = zeros(3,nvar);
% terms from glycine
A(1,7) = -1;
A(2,8) = -1;
A(3,9) = -1;
% terms from serine
A(1,[1,3]) = 1;
A(2,[2,4]) = 1;
A(3,[5,6]) = 1;

% create the initial guess vectors
if(numel(x0_opt) > 0)
    nstart = 1;
    x0 = x0_opt;
else
    nstart = 10;
    x0 = rand(nvar, 10);
end

% perform the optimization
kn_opt_E = knitro_options('UseParallel','true');
for iter = 1:nstart
 [ x(:,iter), fval(iter), exitflag(iter)] = knitro_nlp(obj_fn, x0(:,iter), [],[],A,b,lb,ub,[],[],kn_opt_E,'kn_methf.opt');
end

% select the solution with the minimum error
[~,min_index] = min(fval);
x_opt = x(:,min_index(1));

% calculate the MID of MeTHF
cthf_mid0 = x_opt(1) + x_opt(2) + x_opt(5);
cthf_mid1 = x_opt(3) + x_opt(4) +x_opt(6);
cthf_mid = [cthf_mid0; cthf_mid1];


end

function obj = methf_objective(x, ser_MID, gly_MID, ser_sd, gly_sd)
    
    ser00 = x(1);
    ser10 = x(2);
    ser01 = x(3);
    ser11 = x(4);
    ser20 = x(5);
    ser21 = x(6);

    gly0 = x(7);
    gly1 = x(8);
    gly2 = x(9);

    ser0 = ser00;
    ser1 = ser10 + ser01;
    ser2 = ser11 + ser20;
    ser3 = ser21;

    ser_sim = [ser0; ser1; ser2; ser3];
    gly_sim = [gly0; gly1; gly2];

    error = [ser_sim; gly_sim] - [ser_MID; gly_MID];
    sd = [ser_sd; gly_sd];
    obj = sumsqr(error./sd);    

end