function v = fluxTransform1(x, P, nRev)

vNet = x(1:2:nRev-1);
vXchTemp = x(2:2:nRev);
vXch = P .* vXchTemp ./ (1 - vXchTemp);
zer = zeros(length(vNet), 1);
vFwd =  vXch - min(-vNet, zer);
vRev = vXch - min(vNet, zer);
v = zeros(nRev, 1);
v(1:2:nRev-1) = vFwd;
v(2:2:nRev) = vRev;