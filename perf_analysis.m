%% Performance analysis calcs
V0 = 0;
x0 = 0;

y0 = [x0; V0];

% Time span
tspan = [0 15];

rpm_100 = 9989;
rpm_80 = 0.80 * rpm_100;
rpm_68 = 0.68 * rpm_100;
rpm_60 = 0.6 * rpm_100;

%% Takeoff Performance
[t,y] = ode45(@(t,y) takeoff_ode(t,y,weight_TO, rho, wing_area_total, mu_ground, CL_0, CD_0, K_w,Thrust_func,K_t,Sh_S,CL0_t,CLa_t,CL_trim,CL_clean,CD_clean,CD_trim) , tspan, y0);

dist = y(:,1);
Vel = y(:,2);
t_TO = interp1(Vel, t, V_TO); % Interpolate to find exact takeoff time
x_TO = interp1(t, dist, t_TO); % Corresponding takeoff distance
    
figure();
yyaxis left
plot(t, Vel, "b-",'LineWidth', 1.5)
hold on
ylabel('Velocity (m/s)','Color',"b")
yline(V_TO, "k--")
ax = gca;
ax.YColor = 'b'; % Changes y-axis line and tick marks (RGB)
yyaxis right
plot(t, dist, "m-", 'LineWidth', 1.5);
hold on
xline(t_TO, 'k--', 'LineWidth', 1.5)
plot(t, dist, "m-", 'LineWidth', 1.5);
ylabel('Distance (m)')
ax = gca;
ax.YColor = 'm';
title('Aircraft Takeoff Performance')
xlabel('Time (s)')
grid on
legend("Velocity","TO Velocity","Distance","Takeoff Location",Location="southeast")

% Force Sanity Check
% % Thrust
% T = zeros(size(Vel));
% for i = 1:length(Vel)
%         T(i) = Thrust_func(Vel(i),9989);
% end   
% Lift
% L = 0.5*rho*Vel.^2*wing_area_total*CL_0;
% 
% % Drag
% D = 0.5*rho*Vel.^2*wing_area_total * ...
%         (interp1(CL_clean, CD_clean, CL_0)); % ask on Friday about this
% 
% % Friction
% Ff = mu_ground*(weight_TO - L);
% figure() 
% plot(Vel, T)
% hold on
% plot(Vel,D)
% plot(Vel,Ff)
% plot(Vel, D+Ff)
% plot(Vel,L)
% plot(Vel,ones(size(L)) * weight_TO)
% grid on
% legend("T","D","Ff","D+Ff","L","W")
% title("Force sanity check")
% ylabel("N")
% xlabel("airspeed")

%% Climb & Cruise Performance Calculation
Vel_array = linspace(V_s,35,1000);

for i = 1:length(Vel_array)
    CL_perf(i) = weight_TO ./ ( 1/2 * rho * Vel_array(i).^2 * wing_area_total);
    CD_perf_clean(i) = interp1(CL_clean, CD_clean, CL_perf(i));
    CD_perf_trim(i) = interp1(CL_clean,CD_trim, CL_perf(i));

    thrust_100(i) = Thrust_func(Vel_array(i),rpm_100);
    thrust_80(i) = Thrust_func(Vel_array(i),rpm_80);
    thrust_68(i) = Thrust_func(Vel_array(i),rpm_68);
    thrust_60(i) = Thrust_func(Vel_array(i),rpm_60);
end

D_perf = 1/2  * rho * Vel_array.^2 * wing_area_total .* CD_perf_trim;

Thrust_req_trim = 1/2 * rho * Vel_array.^2 * wing_area_total .* (CD_perf_trim);
Thrust_req_clean = 1/2 * rho * Vel_array.^2 * wing_area_total .* (CD_perf_clean);

ROC_100 = (thrust_100 - D_perf) .* Vel_array / weight_TO;
ROC_80 = (thrust_80 - D_perf) .* Vel_array / weight_TO;
ROC_60 = (thrust_60 - D_perf) .* Vel_array / weight_TO;

gamma_climb_100 = asind((thrust_100 - D_perf) / weight_TO);
gamma_climb_80 = asind((thrust_80 - D_perf) / weight_TO);
gamma_climb_60 = asind((thrust_60 - D_perf) / weight_TO);

%% Climb Performance Plot
figure();
yyaxis left
plot(Vel_array, ROC_100,"b-",'LineWidth', 1.5)
hold on
plot(Vel_array, ROC_80,"b--",'LineWidth', 1.5)
plot(Vel_array,ROC_60,"b:",'LineWidth', 1.5)
ylabel('Rate of Climb (m/s)')
grid on
ax = gca;
ax.YColor = 'b'; % Changes y-axis line and tick marks (RGB)

yyaxis right
plot(Vel_array, gamma_climb_100, "m-", 'LineWidth', 1.5);
hold on
plot(Vel_array,gamma_climb_80,"m--", 'LineWidth', 1.5)
plot(Vel_array,gamma_climb_60, "m:", 'LineWidth', 1.5)
xlabel('Airspeed (m/s)')
ylabel('Climb Angle (degrees)')
title('Aircraft Climb Performance')
xlim([5 25])
xline(V_s,"r-.",LineWidth=1.5)
legend("RoC at 100%","ROC at 80%","RoC at 60%","Angle at 100%","Angle at 80%","Angle at 60%","Stall Velocity")
ax = gca;
ax.YColor = 'm';

delta_h_climb = 30; % height for cruising altitude in m
% delta t climb at different velocities
gamma_climb_100(gamma_climb_100 <= 3) = 3;
gamma_climb_80(gamma_climb_80 <= 3) = 3;
gamma_climb_60(gamma_climb_60 <= 3) = 3;

delta_t_climb_100 = delta_h_climb./gamma_climb_100;
delta_t_climb_80 = delta_h_climb./gamma_climb_80;
delta_t_climb_60 = delta_h_climb./gamma_climb_60;

figure();
plot(Vel_array, delta_t_climb_100, "b-", linewidth = 1.5);
hold on
plot(Vel_array, delta_t_climb_80, "m-", linewidth = 1.5);
plot(Vel_array, delta_t_climb_60, "g-", linewidth = 1.5);
grid on
xline(V_CL);
xlabel("Velocity (m/s)");
ylabel("delta t climb (s)");
legend("t climb - 100 % throttle", "t climb - 80 %", "t climb - 60 %", "Takeoff Velocity");
title("t climb vs velocity");
xlim([5 25]);



%% Cruise performace plot
figure()
plot(Vel_array,Thrust_req_trim,"k-",'LineWidth', 1.5)
hold on
plot(Vel_array,Thrust_req_clean,"k--",'LineWidth', 1.5)
plot(Vel_array, thrust_100,"b-",'LineWidth', 1.5)
plot(Vel_array, thrust_80,"b--",'LineWidth', 1.5)
plot(Vel_array,thrust_68,"b-.",'LineWidth', 1.5)
title("Aircraft Cruise Performance - Thrust")
ylabel("Thrust/Drag (N)")
xlabel("Airspeed (m/s)")
grid on
xline(V_s,"r-.",LineWidth=1.5)
legend("Thrust Required Trimmed (Drag)","Thrust Required Clean (Drag)","Thrust at 100%","Thrust at 80%","Thrust at 68%","Stall Velocity")

for i = 1:length(Vel_array)
np_100perf(i) = etaP_func(Vel_array(i),rpm_100);
np_80perf(i) = etaP_func(Vel_array(i),rpm_80);
np_68perf(i) = etaP_func(Vel_array(i),rpm_68);
end

np_100perf(np_100perf <= 0.1) = 0.1;
np_80perf(np_80perf <= 0.1)   = 0.1;
np_68perf(np_68perf <= 0.1)   = 0.1;

figure()
plot(Vel_array,(Thrust_req_trim .* Vel_array)/(etaESC*etaM*etaP_CR),"k-",'LineWidth', 1.5)
hold on
plot(Vel_array,(Thrust_req_clean .* Vel_array)/(etaESC*etaM*etaP_CR),"k--",'LineWidth', 1.5)
plot(Vel_array, (thrust_100 .* Vel_array)./(etaESC*etaM*np_100perf),"b-",'LineWidth', 1.5)
plot(Vel_array, (thrust_80 .* Vel_array)./(etaESC*etaM*np_80perf),"b--",'LineWidth', 1.5)
plot(Vel_array,(thrust_68 .* Vel_array)./(etaESC*etaM*np_68perf),"b-.",'LineWidth', 1.5)
xlim([5 28])
ylim([0 750])
title("Aircraft Cruise Performance - Power")
ylabel("Power (W)")
xlabel("Airspeed (m/s)")
grid on
xline(V_s,"r-.",LineWidth=1.5)
legend("Power Required Trimmed (Drag)","Power Required Clean (Drag)","Power at 100%","Power at 80%","Power at 68%","Stall Velocity")


function dydt = takeoff_ode(t,y,weight_TO, rho, wing_area_total, mu_ground, CL_0, CD_0, K_w, ...
    Thrust_func,K_t,Sh_S,CL0_t,CLa_t,CL_trim,CL_clean,CD_clean,CD_trim) 

    % velocity
    dist = y(1);     % position
    Vel = y(2); 
    
    %% Thrust
    for i = 1:length(Vel)
        T(i) = Thrust_func(Vel(i),9989);
    end
    
    %% Lift
    L = 0.5*rho*Vel.^2*wing_area_total*CL_0;
    
    %% Drag
    D = 0.5*rho*Vel.^2*wing_area_total * ...
        (interp1(CL_clean, CD_clean, CL_0)); % ask on Friday about this
    
    %% Friction
    Ff = mu_ground*(weight_TO - L);
    
    %% Acceleration
    m = weight_TO/9.81;
    a_long = (T - D - Ff)/m;
    
    %% ODE system
    dydt = zeros(2,1);
    dydt(1) = Vel;        % dx/dt
    dydt(2) = a_long;   % dV/dt
end

%% Range and Endurance calculations
E_batt_wh = 4*18.5; % Amp hours 8 voltage = Wh
E_batt_J  = E_batt_wh * 3600; % Convert to Joules

% Throttle sweep
throttle_array = linspace(0.4,1.0,1000);   % adjust lower limit if needed

% Preallocate
V_cruise   = zeros(size(throttle_array));
P_elec = zeros(size(throttle_array));
Range  = zeros(size(throttle_array));

for i = 1:length(throttle_array)

    throttle = throttle_array(i);
    rpm = throttle * rpm_100;

    % Create velocity sweep to check feasibility
    V_test = linspace(V_s+0.5, 35, 200);

    Drag_test = 0.5*rho*V_test.^2*wing_area_total .* ...
        interp1(CL_clean, CD_trim, ...
        weight_TO./(0.5*rho*V_test.^2*wing_area_total), ...
        'linear','extrap');

    Thrust_test = arrayfun(@(Velo) Thrust_func(Velo,rpm), V_test);

    diff_test = Thrust_test - Drag_test;

    % Check if sign change exists
    if max(diff_test) < 0
        % Not enough thrust for cruise
        V_cruise(i) = NaN;
        P_elec(i) = NaN;
        Range(i) = NaN;
        continue
    end

    % Now solve safely
    cruise_equation = @(Velo) Thrust_func(Velo,rpm) - ...
        (0.5*rho*Velo.^2*wing_area_total .* ...
        interp1(CL_clean, CD_trim, ...
        weight_TO./(0.5*rho*Velo.^2*wing_area_total), ...
        'linear','extrap'));

    V_cruise(i) = fzero(cruise_equation, [V_s+0.5 35]);

    if V_cruise(i) <= V_s
        V_cruise(i) = NaN;
        continue
    end

    P_shaft = pow_func(V_cruise(i), rpm);
    P_elec(i) = P_shaft / (etaESC * etaM);

    Range(i) = V_cruise(i) * (E_batt_J / P_elec(i));

end
Endurance_sec = E_batt_J ./ P_elec;   % seconds
Endurance_min = Endurance_sec ./ 60;  % minutes

% Convert range to km
Range_km = Range ./ 1000;

%% Remove invalid points
valid = ~isnan(Range_km) & Range_km > 0;

%% Plot Range vs Throttle
figure
plot(100*throttle_array(valid), Range_km(valid), 'b-', 'LineWidth', 1.5)
grid on
xlabel('Throttle (%)')
ylabel('Range (km)')
title('Range vs Throttle Setting')

figure()
plot(Vel_array, Range_km, "m-", linewidth = 1.5);
grid on
hold on
yline(max(Range_km));
title("Range vs Velocity");
legend("Range vs Velocity", "Max Range");


%% Endurance vs Throttle
figure
plot(100*throttle_array(valid), Endurance_min(valid), 'r-', 'LineWidth', 2)
grid on
xlabel('Throttle (%)')
ylabel('Endurance (minutes)')
title('Endurance vs Throttle Setting')

figure()
plot(Vel_array, Endurance_min, "g-", linewidth = 1.5);
grid on
hold on
yline(max(Endurance_min));
title("Endurance vs Velocity");
legend("Endurance vs Velocity", "Max Endurance");

% Throttles you want displayed
display_throttles = [1.0, 0.80, 0.68, 0.60];

fprintf('\nAircraft Range for Sustainable Cruise Conditions:\n')

for k = 1:length(display_throttles)

    throttle = display_throttles(k);
    rpm = throttle * rpm_100;

    % --- Solve cruise condition (Thrust = Drag) ---
    cruise_eq = @(Velo) Thrust_func(Velo,rpm) - ...
        (0.5*rho*Velo.^2*wing_area_total .* ...
        interp1(CL_clean, CD_trim, ...
        weight_TO./(0.5*rho*Velo.^2*wing_area_total), ...
        'linear','extrap'));

    V_test = linspace(V_s+0.5, 35, 200);
    Drag_test = 0.5*rho*V_test.^2*wing_area_total .* ...
        interp1(CL_clean, CD_trim, ...
        weight_TO./(0.5*rho*V_test.^2*wing_area_total), ...
        'linear','extrap');

    Thrust_test = arrayfun(@(Velo) Thrust_func(Velo,rpm), V_test);

    if max(Thrust_test - Drag_test) < 0
        fprintf('Throttle %.0f%%: Not enough thrust for cruise\n\n', throttle*100)
        continue
    end

    V_cruise = fzero(cruise_eq, [V_s+0.5 35]);

    % --- Power ---
    P_shaft = pow_func(V_cruise, rpm);
    P_elec = P_shaft / (etaESC * etaM);

    % --- Endurance ---
    Endurance_sec = E_batt_J / P_elec;
    Endurance_hr  = Endurance_sec / 3600;

    % --- Range ---
    Range_km = (V_cruise * Endurance_sec) / 1000;

    % --- Display ---
    fprintf('Throttle %.0f%%:\n', throttle*100)
    fprintf('  - Speed: %.2f m/s\n', V_cruise)
    fprintf('  - Endurance: %.2f hours\n', Endurance_hr)
    fprintf('  - Range: %.2f km\n\n', Range_km)

end