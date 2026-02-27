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

figure()
plot(Vel_array,Thrust_req_trim .* Vel_array,"k-",'LineWidth', 1.5)
hold on
plot(Vel_array,Thrust_req_clean .* Vel_array,"k--",'LineWidth', 1.5)
plot(Vel_array, thrust_100 .* Vel_array,"b-",'LineWidth', 1.5)
plot(Vel_array, thrust_80 .* Vel_array,"b--",'LineWidth', 1.5)
plot(Vel_array,thrust_68 .* Vel_array,"b-.",'LineWidth', 1.5)
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

