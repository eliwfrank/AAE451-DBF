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

CL_R = (2*weight_TO)/(rho*wing_area_total*(V_TO^2));

cdi_array = (cl_array.^2)/(pi()*AR_wing*e);
D_i = cdi_array.*0.5*rho.*(velocity_array.^2);
D_0 = 0.5*CD_0*rho.*((velocity_array).^2).*wing_area_total;
D_array = D_0 + D_i;
% values from chosen prop system - 12X6E at 9,989 RPM
rpm_100per = 9989;
rpm_75per = 0.75 * rpm_100per;
rpm_50per = 0.50 * rpm_100per;
rpm_25per = 0.25 * rpm_100per;

for i = 1:length(velocity_array)
thrust_100per(i) = Thrust_func(velocity_array(i),rpm_100per);
thrust_75per(i) = Thrust_func(velocity_array(i),rpm_75per);
thrust_50per(i) = Thrust_func(velocity_array(i),rpm_50per);
thrust_25per(i) = Thrust_func(velocity_array(i),rpm_25per);
end

% Find intersection of 100% throttle thrust and drag

T_interp = interp1(velocity_array, thrust_100per, ...
                   velocity_array, 'linear', 'extrap');

diff_T_D = T_interp - D_array;

idx = find(diff_T_D(1:end-1).*diff_T_D(2:end) <= 0, 1, 'first');

V_max = interp1(diff_T_D(idx:idx+1), velocity_array(idx:idx+1), 0);
fprintf("Maximum Velocity (Thrust = Drag): %.4f [m/s]\n", V_max);

figure();
plot(velocity_array, thrust_100per, "m-", linewidth = 1.5);
hold on
plot(velocity_array, thrust_75per, "m-.", linewidth = 1.5);
plot(velocity_array, thrust_50per, "m:", linewidth = 1.5);
plot(velocity_array, thrust_25per, "m--", linewidth = 1.5);
plot(velocity_array, D_array, "g-", linewidth = 1.5);
plot(V_max, interp1(velocity_array, D_array, V_max), ...
     'bo', 'MarkerSize', 8, linewidth = 2);
ylim([0 30]);
title("Drag and Thrust vs Airspeed");
ylabel("Forces (N)");
xlabel('Airspeed (m/s)');
legend("100 % Throttle", "75 % Throttle", "50 % Throttle", "25 % Throttle", "Drag");
xlim([0 28]);
grid on;

% plotting propeller efficiency vs airspeed


for i = 1:length(velocity_array)
np_100per(i) = etaP_func(velocity_array(i),rpm_100per);
np_75per(i) = etaP_func(velocity_array(i),rpm_75per);
np_50per(i) = etaP_func(velocity_array(i),rpm_50per);
np_25per(i) = etaP_func(velocity_array(i),rpm_25per);
end

figure();
plot(velocity_array, np_100per, "c-", linewidth = 1.5);
hold on
plot(velocity_array, np_75per, "c-.", linewidth = 1.5);
plot(velocity_array, np_50per, "c:", linewidth = 1.5);
plot(velocity_array, np_25per, "c--", linewidth = 1.5);
ylim([0 1]);
ylabel("Propeller Efficiency");
xlabel('Airspeed (m/s)');
title("Propeller Efficiency vs Airspeed");
legend("100 % Throttle", "75 % Throttle", "50 % Throttle", "25 % Throttle");
grid on;

fprintf("Thrust: %f [N]\n", T_cr)

for i = 1:length(velocity_array)
pow_100per(i) = (pow_func(velocity_array(i),rpm_100per))/(etaM*etaESC);
pow_75per(i) = (pow_func(velocity_array(i),rpm_75per))/(etaM*etaESC);
pow_50per(i) = (pow_func(velocity_array(i),rpm_50per))/(etaM*etaESC);
pow_25per(i) = (pow_func(velocity_array(i),rpm_25per))/(etaM*etaESC);
end


figure();
plot(velocity_array, pow_100per, 'b-', 'LineWidth', 1.5);
hold on
plot(velocity_array, pow_75per,  'b-.', 'LineWidth', 1.5);
plot(velocity_array, pow_50per,  'b:', 'LineWidth', 1.5);
plot(velocity_array, pow_25per,  'b--', 'LineWidth', 1.5);
ylim([0 700]);
xlabel('Airspeed (m/s)');
ylabel('Electrical Power (W)');
title('Electrical Power Required');
legend('100% Throttle','75% Throttle','50% Throttle','25% Throttle');
grid on;

etaP_CR = etaP_func(V_CR, rpm_100per);

etaP_TO = etaP_func(V_TO, rpm_100per);

etaP_CL = etaP_func(V_CL, rpm_100per);

etaP_M = etaP_func(V_M, rpm_100per);