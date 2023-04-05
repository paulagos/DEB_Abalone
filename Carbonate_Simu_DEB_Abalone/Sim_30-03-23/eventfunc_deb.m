function [value, isterminal, direction] = eventfunc_deb(t, EVHRSO, simu)

%------------ Value -------------------------------
%
% EVENT 1: Egg to Larva transition
par   = simu.par;
E_Hb  = par(13);  % J, E_H^b, maturity at birth
E_Hj  = par(28);  % J, E_H^J, maturity at metamorphosis   
E_Hp  = par(14);  % J, E_H^p, maturity at puberty

value(1) = (EVHRSO(3) >= E_Hb);  % egg to larva transition

value(2) = (EVHRSO(3) >= E_Hj); % larva to juvenile transition

value(3) = (EVHRSO(4) >= E_Hp); % Juvenile to adult transition

%---------------- Terminal--------------------------
isterminal = ones(size(value));

% -------------- Direction -------------------------
direction = zeros(size(value)); 

end