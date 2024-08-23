function M = simulate(v, c, n_flux_param, n_c_param, ser1, ser2, ser3, ...
    r5p1, r5p2, r5p3, r5p4, r5p5, ...
    M0, J, tspan, metab_char, metab_size, varFlag)

vg = v(1:n_flux_param); % serine fluxes
vp = v(n_flux_param + 1:end); % purine fluxes

cg = c(1:n_c_param); % serine pools
cp = c(n_c_param + 1:end); % purine pool sizes

% Set ODE options
starttime = datetime('now');
stoptime = starttime + minutes(5);
timer = @(tp, y)EventFcn(tp, y, stoptime);
reltol = 0.001;
abstol = 0.001;
odeOptions = odeset('RelTol', reltol, 'AbsTol', abstol, 'Events', timer);

% Define gradient function
odefunc = @(t,x)ODE_grad(t, x, J, metab_char, metab_size, varFlag, ...
    ser1, ser2, ser3, r5p1, r5p2, r5p3, r5p4, r5p5, ...
    vg, cg, vp, cp);

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

function G = ODE_grad(t, x, J, metab_char, metab_size, varFlag, ...
    ser1, ser2, ser3, r5p1, r5p2, r5p3, r5p4, r5p5, ...
    vg, cg, vp, cp)

G = zeros(numel(x),1);

% Get the serine IDVs
nSer = numel(x) - 36;

ser_idv = x(1:nSer);
ser_idv = J * ser_idv;

% Add the enrichment of SERg
clear m index map
m = find(ismember(metab_char,'SERg'));
index = uIDVindex(m,metab_size,varFlag,0);
ser_idv(index(1)) = 1 - ser1 - ser2 - ser3;
ser_idv(index(1) + 1) = ser1; % isotopomer 001
ser_idv(index(2) - 1) = ser2; % isotopomer 110
ser_idv(index(2)) = ser3;

% Estimate gradient of glycine and MTHF 
G(1:nSer) = serine_ODE(ser_idv, vg, cg);

clear m index map
m = find(ismember(metab_char,'MTHFg'));
index = uIDVindex(m, metab_size, varFlag, 0);
M_mthf = ser_idv(index(1):index(2));

clear m index map
m = find(ismember(metab_char,'GLYg'));
index = uIDVindex(m, metab_size, varFlag, 0);
map = idv_to_mid(2); % no of carbons
ser_idv = ser_idv(index(1):index(2));
M_gly = map * ser_idv;

clear m index map


% Get the MIDs of inputs to the purine pathway
r5p0 = 1 - r5p1 - r5p2 - r5p3 - r5p4 - r5p5;
M_r5p = [r5p0; r5p1; r5p2; r5p3; r5p4; r5p5];


% Gradients for purine metabolites
M_amp = x(nSer + 1:nSer + 6);
M_gdp = x(nSer + 7:nSer + 12);
M_gmp = x(nSer + 13:nSer + 18);
M_guo = x(nSer + 19:nSer + 24);
M_imp = x(nSer + 25:nSer + 30);
M_ino = x(nSer + 31:nSer + 36);

c_amp = cp(1);
c_gdp = cp(2);
c_gmp = cp(3);
c_guo = cp(4);
c_imp = cp(5);
c_ino = cp(6);

imp_in = vp(1);
hpx_r5p_imp = vp(2);
imp_gmp = vp(3);
gua_r5p_gmp = vp(4);
gmp_gdp = vp(5);
gdp_out = vp(6);
gmp_guo = vp(7);
gua_r5p_guo = vp(8);
guo_r5p_gua = vp(9);
imp_ino = vp(10);
imp_amp = vp(11);
amp_imp = vp(12);
hpx_r5p_ino = vp(13);
ino_r5p_hpx = vp(14);
ade_amp = vp(15);
ade_ino = vp(16);
amp_out = vp(17);

% AMP MID gradient
G(nSer + 1:nSer + 6) = (-(amp_out + amp_imp) * M_amp + imp_amp * M_imp + ...
    ade_amp * [1; 0; 0; 0; 0; 0]) / c_amp;
% GDP MID gradient
G(nSer + 7:nSer + 12) = (-gdp_out * M_gdp + gmp_gdp * M_gmp) / c_gdp;
% GMP MID gradient
G(nSer + 13:nSer + 18) = (-(gmp_gdp + gmp_guo) * M_gmp + imp_gmp * M_imp + ...
    gua_r5p_gmp * M_r5p) / c_gmp;
% Guanosine MID gradient
G(nSer + 19:nSer + 24) = (-guo_r5p_gua * M_guo + gua_r5p_guo * M_r5p + ...
    gmp_guo * M_gmp) / c_guo;
% IMP MID gradient
imp_de_novo = sum_diag(M_mthf * M_mthf');
imp_de_novo = sum_diag(imp_de_novo * M_gly');
imp_de_novo = sum_diag(imp_de_novo * M_r5p');
G(nSer + 25:nSer + 30) = (-(imp_ino + imp_gmp + imp_amp) * M_imp + ...
    imp_in * imp_de_novo(1:6) + hpx_r5p_imp * M_r5p + amp_imp * M_amp) / c_imp;
% INOSINE MID gradient
G(nSer + 31:nSer + 36) = (-ino_r5p_hpx * M_ino + imp_ino * M_imp + ...
    hpx_r5p_ino * M_r5p + ade_ino * [1; 0; 0; 0; 0; 0]) / c_ino;

end

function map = idv_to_mid(n)

map = zeros(n + 1, 2 ^ n);
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
