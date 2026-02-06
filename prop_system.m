%% power equations
% iterating through to find battery power required
P_req_mech = W_int/W_P;
P_req_bat = P_req_mech/(etaM*etaESC);

% this is how these would change if our weight increased by 50 percent
P_req_mech_50 = (W_int*1.5)/(W_P);
P_req_bat_50 = P_req_mech_50/(etaM*etaESC);

fprintf("Battery Power Required: %.4f [W]\n",P_req_bat);

%% Thrust Equations
% calculating static thrust CLr = CLTo
L_TO =  0.7*V_TO; % need to calculate this, this is an estimate
D_TO = 0.7*V_TO; % need to calculate this
T_0 = ((W_int*(V_TO^2))/(2*g*S_TO)) + D_TO + mu*(W_int - L_TO);

fprintf("Estimated Static Thrust: %.4f [N]\n",T_0);

% calculating cruise thrust

% T_cr = D_cr;

