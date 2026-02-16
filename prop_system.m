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
velocity_array = linspace(0,28,1000);
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
% values from chosen prop system - 11X55E at 12,000 RPM
velocity_fromprop_10 = [0.000 
1.033
2.065
3.098
4.131
5.163
6.196
7.230
8.263
9.296
10.328
11.363
12.396
13.429
14.462
15.494
16.528
17.561
18.594
19.627
20.659
21.693
22.726
23.759
24.792
25.824
26.858
27.891
28.924];

velocity_fromprop_11 = [0.000
1.136
2.275
3.411
4.547
5.682
6.821
7.957
9.093
10.232
11.368
12.504
13.640
14.779
15.915
17.051
18.190
19.326
20.462
21.598
22.737
23.873
25.009
26.148
27.284
28.420
29.559
30.695
31.831];

thrust_fromprop_10 = [20.372
20.026
19.655
19.259
18.836
18.385
17.906
17.397
16.858
16.288
15.688
15.056
14.394
13.701
12.977
12.223
11.440
10.632
9.802
8.952
8.085
7.204
6.311
5.408
4.499
3.588
2.678
1.775
0.882];

thrust_fromprop_11 = [24.786
24.366
23.917
23.436
22.923
22.375
21.793
21.175
20.519
19.826
19.096
18.327
17.521
16.677
15.796
14.878
13.924
12.940
11.928
10.893
9.837
8.763
7.674
6.575
5.467
4.357
3.248
2.147
1.060];

% interpolation for rpm
rpm_low  = 10000;
rpm_high = 11000;   % example

rpm_target = 10721;

interp = (rpm_target - rpm_low) / (rpm_high - rpm_low);

V_10721 = velocity_fromprop_10 + interp * (velocity_fromprop_11 - velocity_fromprop_10);
V_75per = V_10721.*0.75;
V_50per = V_10721.*0.50;
V_25per = V_10721.*0.25;

T_10721  = thrust_fromprop_10  + interp * (thrust_fromprop_11  - thrust_fromprop_10);
T_75per = T_10721.*0.75;
T_50per = T_10721.*0.50;
T_25per = T_10721.*0.25;

figure();
plot(V_10721, T_10721, "m-", linewidth = 1.5);
hold on
plot(V_75per, T_75per, "m-.", linewidth = 1.5);
plot(V_50per, T_50per, "m:", linewidth = 1.5);
plot(V_25per, T_25per, "m--", linewidth = 1.5);
plot(velocity_array, D_array, "g-", linewidth = 1.5);
title("Drag and Thrust vs Airspeed");
ylabel("Forces (N)");
xlabel('Airspeed (m/s)');
legend("100 % Throttle", "75 % Throttle", "50 % Throttle", "25 % Throttle", "Drag");
xlim([0 28]);
grid on;

% plotting propeller efficiency vs airspeed
n_p_fromprop_10 = [0.0000
0.0630
0.1226
0.1790
0.2322
0.2823
0.3294
0.3736
0.4149
0.4535
0.4894
0.5226
0.5531
0.5810
0.6062
0.6286
0.6480
0.6643
0.6772
0.6861
0.6906
0.6898
0.6824
0.6668
0.6403
0.5987
0.5350
0.4366
0.2785];

n_p_fromprop_11 = [0.0000
0.0633
0.1232
0.1799
0.2333
0.2837
0.3310
0.3755
0.4170
0.4559
0.4920
0.5255
0.5563
0.5845
0.6101
0.6328
0.6526
0.6694
0.6827
0.6922
0.6973
0.6972
0.6906
0.6759
0.6502
0.6094
0.5462
0.4474
0.2859];

np_10721 = n_p_fromprop_10 + interp * (n_p_fromprop_11 - n_p_fromprop_10); % 100% throttle
np_75per = np_10721.*0.75;
np_50per = np_10721.*0.50;
np_25per = np_10721.*0.25;

figure();

plot(V_10721, np_10721, "c-", linewidth = 1.5);
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

P_mech_min_100 = (T_10721.*V_10721)./(np_10721);
P_elec_min_100 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

P_mech_min_75 = (T_10721.*V_10721)./(np_10721);
P_elec_min_75 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

P_mech_min_50 = (T_10721.*V_10721)./(np_10721);
P_elec_min_50 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

P_mech_min_25 = (T_10721.*V_10721)./(np_10721);
P_elec_min_25 = (P_mech_min_100)./(etaM*etaESC); % converting mechanical power to electrical

figure();
plot(V_10721, P_elec_min_100, "b-", LineWidth=1.5);
hold on
plot(V_75per, P_elec_min_75, "b-.", LineWidth=1.5);
plot(V_50per, P_elec_min_50, "b:", LineWidth=1.5);
plot(V_25per, P_elec_min_25, "b--", LineWidth=1.5);
xlabel('Airspeed (m/s)');
ylabel("Power (W)");
legend("100 % Throttle", "75 % Throttle", "50 % Throttle", "25 % Throttle");
grid on
title("Electrical Power vs Airspeed");