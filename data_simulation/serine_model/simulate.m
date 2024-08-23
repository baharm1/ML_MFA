function M = simulate(v, c, n_flux_param, n_c_param, serp1, serp2, serp3, ...
    pgc1, pgc2, pgc3, pg1, pg2, pg3, M0, J, tspan, metab_char, metab_size, varFlag)

vg = v(1:n_flux_param); % serine model fluxes

cg = c(1:n_c_param); % serine model pools

% Set ODE options
starttime = datetime('now');
stoptime = starttime+minutes(5);
timer = @(tp,y)EventFcn(tp,y,stoptime);
reltol = 0.001;
abstol = 0.001;
odeOptions = odeset('RelTol', reltol, 'AbsTol', abstol, 'Events', timer);

% Define gradient function
odefunc = @(t,x)ODE_grad(t, x, J, metab_char, metab_size, varFlag, serp1, ...
    serp2, serp3, pgc1, pgc2, pgc3, pg1, pg2, pg3, vg, cg);

[t,M] = ode23s(odefunc, tspan, M0, odeOptions);

end

function [value, isterminal, direction] = EventFcn(t, y, stoptime)

    isterminal = 1;
    direction = 0;
    if(datetime('now') > stoptime)
        value = 0;
        disp('ODE solver terminated');
    else
        value = 1;
    end

end

function G = ODE_grad(t, x, J, metab_char, metab_size, varFlag, serp1, ...
    serp2, serp3, pgc1, pgc2, pgc3, pg1, pg2, pg3, vg, cg)

G = zeros(numel(x), 1);

% Get the serine IDVs
nSer = numel(x);
ser_idv = x(1:nSer);
ser_idv = J * ser_idv;

% Add the enrichment of circulating serine
m = find(ismember(metab_char, 'SERp'));
index = uIDVindex(m, metab_size, varFlag, 0);
ser_idv(index(1)) = 1 - serp1 - serp2 - serp3;
ser_idv(index(1) + 1) = serp1; % isotopomer 001
ser_idv(index(2) - 1) = serp2; % isotopomer 110
ser_idv(index(2)) = serp3;

% Add the enrichment of PGc
m = find(ismember(metab_char,'PGic'));
index = uIDVindex(m, metab_size, varFlag, 0);
ser_idv(index(1)) = 1 - pgc1 - pgc2 - pgc3;
ser_idv(index(1) + 1) = pgc1;
ser_idv(index(2) - 1) = pgc2;
ser_idv(index(2)) = pgc3;

% Add the enrichment of PGc
m = find(ismember(metab_char,'PGig'));
index = uIDVindex(m, metab_size, varFlag, 0);
ser_idv(index(1)) = 1 - pg1 - pg2 - pg3;
ser_idv(index(1) + 1) = pg1;
ser_idv(index(2) - 1) = pg2;
ser_idv(index(2)) = pg3;

% Estimate gradient of serine model metabolites
G(1:nSer) = serine_ODE(ser_idv, vg, cg);

end

function map = idv_to_mid(n)

map = zeros(n+1,2^n);
  for r = 1:size(map, 2)
     map(sum(dec2bin(r - 1) == '1') + 1, r) = 1;
  end

end

function sum = sum_diag(N)

[m, n] = size(N);
sum = zeros(m + n - 1, 1);
for i = 1:m
    for j = 1:n
        sum(i + j - 1) = sum(i + j - 1) + N(i, j);
    end
end

end
