%% power equations
% iterating through to find battery power required
P_req_mech = weight_TO/W_P;
P_req_bat = P_req_mech/(etaM*etaESC);

% this is how these would change if our weight increased by 50 percent
P_req_mech_50 = (weight_TO*1.5)/(W_P);
P_req_bat_50 = P_req_mech_50/(etaM*etaESC);
fprintf("\nPOWER REQUIREMENTS ------------------------\n")

fprintf("Battery Power Required: %.4f [W]\n",P_req_bat);

%% Thrust Equations
% calculating static thrust CLr = CLTo
L_TO =  0.7*V_TO; % estimate for static thrust
D_TO = 0.7*V_TO; 
T_0 = ((weight_TO*(V_TO^2))/(2*g*S_TO)) + D_TO + mu_ground*(weight_TO - L_TO);

fprintf("Estimated Static Thrust: %.4f [N]\n",T_0);

CL_CR = weight_TO/(0.5*rho*(V_CR^2)*wing_area_total);
CD_i = (CL_CR^2)/(pi()*AR_wing*e);
% calculating cruise thrust
C_D = CD_0 + CD_i; % add aditional C_Ds after calcualtion, cdi and cdwave
D_cr = 0.5*C_D*rho*(V_CR^2)*wing_area_total;
T_cr = D_cr;

% plotting thrust and drag vs airspeed
velocity_array = linspace(0,30,1000);
% for loop for calculating the cl array at different velocities
cl_array = zeros(size(velocity_array)); 
for i = 1:length(velocity_array)
    if velocity_array(i) <= V_TO % for before v takeoff make the value constant
        cl_array(i) = (2*weight_TO)/(rho*wing_area_total*(V_TO^2));
    else
        cl_array(i) = (2 * weight_TO) / ...
                      ((rho * (velocity_array(i)^2)) * wing_area_total); % calculating cl
    end
end

cdi_array = (cl_array.^2)/(pi()*AR_wing*e);
D_i = cdi_array.*0.5*rho.*(velocity_array.^2);
D_0 = 0.5*CD_0*rho.*((velocity_array).^2).*wing_area_total;
D_array = D_0 + D_i;
% values from chosen prop system - 12X6E at 9,989 RPM
velocity_fromprop_10 = [0.000
1.127
2.253
3.384
4.510
5.637
6.763
7.897
9.022
10.148
11.274
12.410
13.533
14.659
15.785
16.910
18.036
19.162
20.295
21.421
22.543
23.673
24.800
25.927
27.053
28.180
29.307
30.435
31.571];

thrust_fromprop_10 = [27.785
27.328
26.837
26.311
25.747
25.145
24.503
23.819
23.092
22.322
21.508
20.649
19.747
18.801
17.813
16.781
15.708
14.598
13.457
12.287
11.093
9.879
8.648
7.404
6.152
4.898
3.647
2.405
1.181];

V_75per = velocity_fromprop_10.*0.75;
V_50per = velocity_fromprop_10.*0.50;
V_25per = velocity_fromprop_10.*0.25;

T_75per = thrust_fromprop_10.*0.75;
T_50per = thrust_fromprop_10.*0.50;
T_25per = thrust_fromprop_10.*0.25;

% Find intersection of 100% throttle thrust and drag

T_interp = interp1(velocity_fromprop_10, thrust_fromprop_10, ...
                   velocity_array, 'linear', 'extrap');

diff_T_D = T_interp - D_array;

idx = find(diff_T_D(1:end-1).*diff_T_D(2:end) <= 0, 1, 'first');

V_max = interp1(diff_T_D(idx:idx+1), velocity_array(idx:idx+1), 0);
fprintf("Maximum Velocity (Thrust = Drag): %.4f [m/s]\n", V_max);

figure();
plot(velocity_fromprop_10, thrust_fromprop_10, "m-", linewidth = 1.5);
hold on
plot(V_75per, T_75per, "m-.", linewidth = 1.5);
plot(V_50per, T_50per, "m:", linewidth = 1.5);
plot(V_25per, T_25per, "m--", linewidth = 1.5);
plot(velocity_array, D_array, "g-", linewidth = 1.5);
plot(V_max, interp1(velocity_array, D_array, V_max), ...
     'bo', 'MarkerSize', 8, linewidth = 2);
title("Drag and Thrust vs Airspeed");
ylabel("Forces (N)");
xlabel('Airspeed (m/s)');
legend("100 % Throttle", "75 % Throttle", "50 % Throttle", "25 % Throttle", "Drag");
xlim([0 28]);
grid on;

% plotting propeller efficiency vs airspeed
n_p_fromprop_10 = [0.0000
0.0646
0.1257
0.1834
0.2377
0.2888
0.3367
0.3817
0.4236
0.4627
0.4991
0.5327
0.5637
0.5919
0.6175
0.6402
0.6601
0.6768
0.6901
0.6996
0.7047
0.7046
0.6980
0.6833
0.6576
0.6168
0.5533
0.4536
0.2898];


np_75per = n_p_fromprop_10.*0.75;
np_50per = n_p_fromprop_10.*0.50;
np_25per = n_p_fromprop_10.*0.25;

figure();
plot(velocity_fromprop_10, n_p_fromprop_10, "c-", linewidth = 1.5);
hold on
plot(V_75per, np_75per, "c-.", linewidth = 1.5);
plot(V_50per, np_50per, "c:", linewidth = 1.5);
plot(V_25per, np_25per, "c--", linewidth = 1.5);
ylabel("Propeller Efficiency");
xlabel('Airspeed (m/s)');
title("Propeller Efficiency vs Airspeed");
legend("100 % Throttle", "75 % Throttle", "50 % Throttle", "25 % Throttle");
grid on;

fprintf("Thrust: %f [N]\n", T_cr)

% plotting minimum power required vs airspeed

P_mech_min_100 = (thrust_fromprop_10.*velocity_fromprop_10)./(n_p_fromprop_10);
P_elec_min_100 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

P_mech_min_75 = (T_75per.*V_75per)./(np_75per);
P_elec_min_75 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

P_mech_min_50 = (T_50per.*V_50per)./(np_50per);
P_elec_min_50 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

P_mech_min_25 = (T_25per.*V_25per)./(np_25per);
P_elec_min_25 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

figure();
plot(velocity_fromprop_10, P_elec_min_100, "b-", LineWidth=1.5);
hold on
plot(V_75per, P_elec_min_75, "b-.", LineWidth=1.5);
plot(V_50per, P_elec_min_50, "b:", LineWidth=1.5);
plot(V_25per, P_elec_min_25, "b--", LineWidth=1.5);
xlabel('Airspeed (m/s)');
ylabel("Power (W)");
legend("100 % Throttle", "75 % Throttle", "50 % Throttle", "25 % Throttle");
grid on
title("Electrical Power vs Airspeed");