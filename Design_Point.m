fprintf("WING AND POWER LOADING ------------------------\n")
num_vals = 100;

%% Stall Constraint
wingloading_stall = 0.5 * rho * V_si^2 * CL_max;

%% Cruise Constraint
wingloading_CR = linspace(0,num_vals,1000);
powerloading_CR = ((etaP_CRi * phi) / (0.5 * rho * 1.1 * CD_0i * V_CR^3)) * wingloading_CR;

%% Climb Constraint
powerloading_CL = etaP_CLi / (V_CL * ( (1 / (0.866 * LD_max)) + sind(gamma)));
 
%% Maneuver Constraint
wingloading_M = linspace(1e-8,num_vals,1000);
powerloading_M = etaP_Mi ./ (q * V_M * (CD_0i ./ wingloading_M + 1/(pi * AR_wing * e) * (n/q)^2 * wingloading_M));

%% Takeoff Constraint
wingloading_TO = linspace(1e-8,num_vals,1000);
exp_term = exp(0.6 .* rho .* g .* CD_G .* S_TO .* 1 ./ wingloading_TO);
powerloading_TO = (etaP_TOi ./ V_TOi) .* (1 - exp_term) ./ (mu_ground - (mu_ground + CD_G ./ CL_Ri) .* exp_term);

%% Plot
figure()
xline(wingloading_stall,"b")
hold on
grid on
plot(wingloading_CR, powerloading_CR,"g");
yline(powerloading_CL,"r");
plot(wingloading_M, powerloading_M,"m")
plot(wingloading_TO,powerloading_TO,"k")
ylim([0 0.6])

margin = 0.90;
W_S = wingloading_stall * margin;
W_P = interp1(wingloading_TO, powerloading_TO, wingloading_stall) * margin;
plot(W_S,W_P,'ko', 'MarkerFaceColor', "yellow", 'MarkerSize', 10)

fprintf("Wingloading Design Point (W/S): %.4f [N/m^2]\nPowerloading Design Point (W/P): %.4f [N/W]\n", W_S,W_P)

legend("Stall","Cruise","Climb","Maneuver","Takeoff","Design Point",Location="northwest")
title("Aircraft Constraint Diagram")
ylabel("Power Loading (W/P) [N/W]")
xlabel("Wing Loading (W/S) [N/m^2]")

V_s = sqrt(W_S / (.5 * rho * CL_max));
V_TO = n * V_s; % TO velocity [m/s]
