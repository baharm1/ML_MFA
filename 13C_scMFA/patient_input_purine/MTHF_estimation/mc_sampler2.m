function mc_sim = mc_sampler2(MID, SD, metab_list, nsim)

% Create random variables sampled from a truncated normal distribution
% This is the inverse trabsform method
% http://web.michaelchughes.com/research/sampling-from-truncated-normal +
% wikipedia

alpha = -MID./SD;
beta = (1-MID)./SD;

a = normcdf(alpha);
b = normcdf(beta);
u = a + (b-a).*unifrnd(0,1,numel(MID),nsim);
x = norminv(u);

% Sample from the distribution
mc_sim = repmat(MID,1,nsim) + repmat(SD,1,nsim).*x;

options = optimoptions('fmincon','Display','off');
% Make sure that the MIDs sum to 1
% Work on one metabolite at a time
metabs = unique(metab_list);
for m = 1:length(metabs)

   metab_values = mc_sim(ismember(metab_list,metabs(m)),:);
   
   for p = 1:nsim

       x0 = metab_values(:,p);
       x0(1) = 1-sum(x0(2:end)); % Subtract the M+i from 1 to get M+0
       if(x0(1) < 0)
        xopt = fmincon(@(x)sumsqr(x-x0), x0,[],[],ones(1,numel(x0)),[1],zeros(numel(x0),1),ones(numel(x0),1),[],options );
       else
           xopt = x0;
       end
       metab_values(:,p) = xopt;
   end
    mc_sim(ismember(metab_list,metabs(m)),:) = metab_values;

end