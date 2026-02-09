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
L_TO =  0.7*V_TO; % estimate for static thrust
D_TO = 0.7*V_TO; 
L_TO =  0.7*V_TO; % need to calculate this, this is an estimate
D_TO = 0.7*V_TO; % need to calculate this
T_0 = ((W_int*(V_TO^2))/(2*g*S_TO)) + D_TO + mu_ground*(W_int - L_TO);

fprintf("Estimated Static Thrust: %.4f [N]\n",T_0);

CL_CR = W_int/(0.5*rho*(V_CR^2)*wing_area_total);
CD_i = (CL_CR^2)/(pi()*AR_wing*e);
% calculating cruise thrust
C_D = CD_0 + CD_i; % add aditional C_Ds after calcualtion, cdi and cdwave
D_cr = 0.5*C_D*rho*(V_CR^2)*wing_area_total;
T_cr = D_cr;

% plotting thrust and drag vs airspeed
velocity_array = linspace(15,60,1000);
cl_array = (2*W_int)./(rho.*(velocity_array.^2).*wing_area_total);
cdi_array = (cl_array.^2)/(pi()*AR_wing*e);
cd_array = CD_0 + cdi_array;
Thrust_available = P_req_mech./velocity_array;
D_array = 0.5*cd_array*rho.*((velocity_array).^2).*wing_area_total;

figure(6);
plot(velocity_array, Thrust_available, "m-", linewidth = 1.5);
hold on
plot(velocity_array, D_array, "g-", linewidth = 1.5);
title("Drag and Thrust vs Airspeed");
ylabel("Forces (N)");
xlabel('Airspeed (m/s)');
legend('Thrust', 'Drag');
grid on;

% plotting minimum power required vs airspeed
P_mech_min = D_array.*velocity_array;
P_elec_min = P_mech_min./(etaESC*etaM);

figure(7);
plot(velocity_array, P_elec_min, "b-", LineWidth=1.5);
xlabel('Airspeed (m/s)');
ylabel("Power (W)");
grid on
title("Electrical Power vs Airspeed");

% plotting propeller efficiency vs airspeed
n_p_array = (T_cr.*velocity_array)./P_elec_min;

figure(8);

plot(velocity_array, n_p_array, "c-", linewidth = 1.5);
ylabel("Propeller Efficiency");
xlabel('Airspeed (m/s)');
title("Propeller Efficiency vs Airspeed");
grid on;
CD_i = CL_CR^2 / (pi() * AR_wing * e);

CD = CD_i + CD_0;
D_cr = 0.5* CD*rho*(V_CR^2)*wing_area_total;
T_cr = D_cr;

fprintf("Thrust: %f [N]\n", T_cr)
