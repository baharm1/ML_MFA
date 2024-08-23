function [C, Ceq, diffC, diffCeq] = NL_const_wrapper(x, nflux, P, nRev, file1)
% Function to generate all nonlinear equations for balancing isotopomers
% of all metabolites
% Returns Ceq matrix 

C = [];
diffC = [];
diffCeq = [];

% Transformed flux vector (from {vNet,vExch} to {vFwd,vBkd})
flux = x(1:nflux);
fluxTrans = fluxTransform1(flux, P, nRev);
flux(1:nRev) = fluxTrans;

IDV = x(nflux + 1:end);
Ceq = f_NL_const(flux, IDV);
