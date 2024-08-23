% Returns first and last index of the metabolite from the vector of known
% IDVs

function index = kIDVindex(metindex,met,varFlag)

index = zeros(1,2);
IDVsizes = 0;
nM = length(met);
if (metindex>nM||varFlag(metindex)==0)
    index = [0 0];
else
    if (metindex==1&&varFlag(1)==1)
    index(1) = 1;
    else
        for i = 1:metindex-1
            IDVsizes = IDVsizes + varFlag(i)*(2^met(i));
            index(1) = 1 + IDVsizes;
        end    
    end
    index(2)=index(1) + 2^met(metindex) - 1;
end