clc
clear
close all

rho = 1.20343; % Air density [kg/m^3]
CD_0 = 0.03; % Parasitic Drag Coefficient
g = 9.81; % [m/s^2]

max_val = 100;

%% Stall Constraint
V_s = 10; % Stall Velocity [m/s]
CL_max = 1.5; % Max Lift Coefficient

wingloading_stall = 0.5 * rho * V_s^2 * CL_max;
powerloading_stall = linspace(0,max_val,1000);

%% Cruise Constraint
eta_CR = 0.7; % Cruise Propellor Efficiency
phi = 0.75; % Cruise Throttle
V_CR = 18; % Cruise Velocity

wingloading_CR = linspace(0,max_val,1000);
powerloading_CR = ((eta_CR * phi) / (0.5 * rho * 1.1 * CD_0 * V_CR^3)) * wingloading_CR;

%% Climb Constraint
eta_CL = 0.55; % Climb Propellor efficiency
V_CL = 13; % Climb Velocity
LD_max = 12;
gamma = 20; % climb angle [deg]

powerloading_CL = eta_CL / (V_CL * (1 / (0.866 * LD_max) + sind(gamma)));

%% Maneuver Constraint
eta_M = 0.7;
V_M = 18; % Maneuver Velocity [m/s]
AR = 7;
e = 0.85; % Oswald Efficiency Factor
n = 2; % Maneuvering load factor
q = 0.5 * rho * V_M^2;

wingloading_M = linspace(0,max_val,1000);
powerloading_M = eta_M ./ (q * V_M * (CD_0 ./ wingloading_M + 1/(pi * AR * e) * (n/q)^2 * wingloading_M));

%% Takeoff Constraint
CD_G = .01; % change!!!!!!!!
S_TO = 18; % TO run [m]
mu = 0.04; % Ground friction coefficient
CL_R = 2.1075; % CL at rotation
eta_TO = 0.5; % TO prop efficiency
V_TO = 12; % TO velocity [m/s]

wingloading_TO = linspace(0,max_val,1000);
exp_term = exp(0.6 .* rho .* g .* CD_G .* S_TO .* 1 ./ wingloading_TO);
powerloading_TO = (eta_TO ./ V_TO) .* (1 - exp_term) ./ (mu - (mu + CD_G ./ CL_R) .* exp_term);

%% Plot
figure()
xline(wingloading_stall,"b")
hold on
grid on
plot(wingloading_CR, powerloading_CR,"g");
yline(powerloading_CL,"r");
plot(wingloading_M, powerloading_M,"m")
plot(wingloading_TO,powerloading_TO,"k")
legend("stall","cruise","climb","maneuver","takeoff",Location="northwest")
title("Aircraft Constraint Diagram")
ylabel("Power Loading (W/P) [N/W]")
xlabel("Wing Loading (W/S) [N/m^2]")


