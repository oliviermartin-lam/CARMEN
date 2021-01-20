function [Cn2eq altEq] = eqLayers(Cn2, altitudes, nEqLayers, power)
%{
            Cn2         ::  The input Cn2 profile (vector)
            altitudes   ::  The input altitudes (vector)
            nEqLayers   ::  The number of output equivalent layers (scalar)
            power       ::  the exponent of the turbulence (default 5/3)
            
            See: Saxenhuber17: Comparison of methods for the reduction of
            reconstructed layers in atmospheric tomography, App Op, Vol. 56, No. 10 / April 1 2017
%}
nCn2 = numel(Cn2);
nAltitudes = numel(altitudes);
if nargin ~= 4
    power = 5/3;
end

% if nargin ~= 5
nSlab = floor(round(nCn2)/fix(nEqLayers));
ppp = 1;
posSlab =  round((linspace(0, nEqLayers-1, nEqLayers))*nSlab)+1;
for iii = 1:nEqLayers-1
    if posSlab(iii) >= posSlab(iii+1)
        posSlab(iii+1) = posSlab(iii)+1;
    end
end
posSlab = [posSlab, nAltitudes+1];
Cn2eq = zeros(1,nEqLayers);
altEq = zeros(1,nEqLayers);
for ii = 1:nEqLayers
    Cn2eq(ii) = sum(Cn2(posSlab(ii):posSlab(ii+1)-1));
    altEq(ii) =  (sum(altitudes(posSlab(ii):posSlab(ii+1)-1) .^ (power) .* Cn2(posSlab(ii):posSlab(ii+1)-1))/Cn2eq(ii)) .^ (1./power);
end
end
