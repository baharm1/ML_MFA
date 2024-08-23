% Function that returns the first and the last index of the IDV vector of
% specified metabolite 'metindex'
% Version 1.3, modified from 1.0 by Abhinav Achreja on 03-31-2013
% Updates: Includes correction for unknown fluxes and IDVs in the same 'x'
% vector

function index = uIDVindex(metindex,met,varFlag,nFlux)

index = zeros(1,2);
IDVsizes = 0;
nM = length(met);

if (metindex>nM||varFlag(metindex)==1)
    index = [0 0];
else
    if (metindex==1&&varFlag(1)==0)
    index(1,1) = 1 + nFlux;
    else
        for i = 1:metindex-1
            IDVsizes = IDVsizes + (1-varFlag(i))*(2^met(i));
        end 
        index(1,1) = 1 + IDVsizes + nFlux;   
    end
    index(1,2)=index(1,1) + 2^met(metindex) - 1;
end