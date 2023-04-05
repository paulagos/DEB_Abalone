% dELHR
function dELHR = dget_ELHR_j(t, ELHR, Lb, Lj, Lm, p_Am, v, g, kJ, kap, kap_R, Hb, Hj, Hp, T_ref, T_A, tT, f)
% Created 21/02/2023 by Paulo F. Lagos

%% Syntax

%% Description
  %  Calculates change in state variables during juvenile and adult stage at given temperature for species with bj acceleration
  %
  % Input:
  %   dELHR = dget_ELHR_pj(t, ELH)
  %   ELHR: 4-vector
  %
  % output:
  %   E: reserve E in j
  %   L: structural length in cm
  %   H: maturity E_H in j
  %   R: reproduction buffer E_R in j
  %   dELHR: change in reserve, length, maturity and reproduction buffer
  %
  %  unpack initial conditions variables
  E = ELHR(1); 
  L = ELHR(2); 
  H = ELHR(3); 
  R = ELHR(4);
  
  % interpolation for TC
  TC    = tempcorr(interp1(tT(:,1),tT(:,2), t), T_ref, T_A);%
  vT    = v * TC;
  kT_J  = kJ * TC; 
  pT_Am = p_Am * TC;
  
  if H < Hb 
    s = 1; % -, multiplication factor for v and {p_Am}
  elseif H < Hj
    s = L / Lb;
  else
    s = Lj / Lb;
  end
  
  % e = vT * E/ L^3/ pT_Am;                   % -, scaled reserve density; [E]/[Em]
  e    = E / pT_Am / vT;                      % -, scaled reserve density; [E]/[Em]
  rT   = s * vT * (e/ L - 1/ Lm/ s)/ (e + g); % 1/d, spec growth rate
  pT_C = E * (s * vT/ L - rT);                % cm^2, scaled mobilisation
  
  % generate dH/dt, dE/dt, dL/dt, dR/dt
  if H < Hp
    dH = (1 - kap) * pT_C - kT_J * H;  
    dR = 0;
  else
    dH = 0;
    dR = kap_R * (1 - kap) * pT_C - kT_J * Hp;
  end
  
  dE = (L > Lb) * s * pT_Am * f * L^2 - pT_C;
  dL = rT * L/3;

  % pack derivatives
  dELHR = [dE; dL; dH; dR];
end        