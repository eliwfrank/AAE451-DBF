%% power equations
% iterating through to find battery power required
P_req_mech = W_int/W_P;
P_req_bat = P_req_mech/(etaM*etaESC);

% this is how these would change if our weight increased by 50 percent
P_req_mech_50 = (W_int*1.5)/(W_P);
P_req_bat_50 = P_req_mech_50/(etaM*etaESC);
fprintf("\nPOWER REQUIREMENTS ------------------------\n")

fprintf("Battery Power Required: %.4f [W]\n",P_req_bat);

%% Thrust Equations
% calculating static thrust CLr = CLTo
L_TO =  0.7*V_TO; % estimate for static thrust
D_TO = 0.7*V_TO; 
T_0 = ((W_int*(V_TO^2))/(2*g*S_TO)) + D_TO + mu_ground*(W_int - L_TO);

fprintf("Estimated Static Thrust: %.4f [N]\n",T_0);

CL_CR = W_int/(0.5*rho*(V_CR^2)*wing_area_total);
CD_i = (CL_CR^2)/(pi()*AR_wing*e);
% calculating cruise thrust
C_D = CD_0 + CD_i; % add aditional C_Ds after calcualtion, cdi and cdwave
D_cr = 0.5*C_D*rho*(V_CR^2)*wing_area_total;
T_cr = D_cr;

% plotting thrust and drag vs airspeed
velocity_array = linspace(0,40,1000);
% for loop for calculating the cl array at different velocities
cl_array = zeros(size(velocity_array)); 
for i = 1:length(velocity_array)
    if velocity_array(i) <= V_TO % for before v takeoff make the value constant
        cl_array(i) = CL_R;
    else
        cl_array(i) = (2 * W_int) / ...
                      (rho * (velocity_array(i)^2) * wing_area_total); % calculating cl
    end
end

cdi_array = (cl_array.^2)/(pi()*AR_wing*e);
cd_array = CD_0 + cdi_array;
% values from chosen prop system - 11X55E at 14,000 RPM
velocity_fromprop = [0.000
1.462
2.928
4.390
5.856
7.318
8.780
10.247
11.709
13.171
14.638
16.100
17.566
19.028
20.490
21.957
23.419
24.881
26.348
27.810
29.276
30.738
32.200
33.667
35.129
36.591
38.058
39.520
40.987];

thrust_fromprop = [40.981
40.291
39.549
38.754
37.903
36.993
36.023
34.990
33.894
32.733
31.507
30.216
28.859
27.439
25.955
24.405
22.797
21.137
19.431
17.685
15.904
14.093
12.258
10.405
8.541
6.673
4.811
2.964
1.143];

D_array = 0.5*cd_array*rho.*((velocity_array).^2).*wing_area_total;

figure();
plot(velocity_fromprop, thrust_fromprop, "m-", linewidth = 1.5);
hold on
plot(velocity_array, D_array, "g-", linewidth = 1.5);
title("Drag and Thrust vs Airspeed");
ylabel("Forces (N)");
xlabel('Airspeed (m/s)');
legend('Thrust', 'Drag');
xlim([0 40]);
grid on;

% plotting minimum power required vs airspeed
P_mech_min = D_array.*velocity_array;
P_elec_min = P_mech_min./(etaESC*etaM);

figure();
plot(velocity_array, P_elec_min, "b-", LineWidth=1.5);
xlabel('Airspeed (m/s)');
ylabel("Power (W)");
grid on
title("Electrical Power vs Airspeed");

% plotting propeller efficiency vs airspeed
n_p_array = (T_cr.*velocity_array)./P_elec_min;

figure();

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
