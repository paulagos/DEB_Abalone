% acidcorr
% computes corections factor for ocean acidificatcion
function [p_Mc, C_X] = acidcorr(pCO2, p, V, TC)
% Created 16/02/2023 by Paulo F. Lagos

%% Syntax
% [p_Mc, C_X] = acidcorr.m*acidcorr*(pCO2, q, C_XH, C_XL, C_MH, C_ML, V, TC)

%% Description
% Calculated two corerction factors for increasing parcial pressure of pCO_2 and multiplies the factors to: 
% (1) C_X apply depression costs associated to decrease in feedings rates with increasing pH
% (2) p_Mc calculates a new manintenance costs (p_Mc) associated to partial increses in partial pressure of pCO_2
%
% Input:
% 
% * p: 4-vector with parameters [q C_XH C_XL C_MH C_ML]
% * pCO2: partiaL pressure of CO_2
% * V: Volumetric lenghts
% * TC: Temperature correction factor

% Output:
%
% * p_Mc vector with factor(s) that affect(s) all rates

if pCO2 <= p(3)
    C_X = 1;
    
elseif (p(3) < pCO2 && pCO2 < p(2))
    C_X = (p(2) - pCO2) / (p(2) - p(3));
    
elseif pCO2 >= p(2)
    C_X = 0;     
end

p_Mc = p(1) * min(max(pCO2 - p(5), 0), p(4) - p(5)) * V^(2/3) * TC;
end
