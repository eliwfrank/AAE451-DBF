% performance analysis calcs

% Initial conditions
V0 = 0;
x0 = 0;

y0 = [V0; x0];

% Time span
tspan = [0 15];

% Solve
[t,y] = ode45(@(t,y) takeoff_ode(t,y,weight_TO, rho, wing_area_total, mu_ground, CL_0, CD_0, K_w,Thrust_func,K_t,Sh_S,CL0_t,CLa_t), tspan, y0);

dist = y(:,1);
Vel = y(:,2);

% Interpolate to find exact takeoff time
t_TO = interp1(Vel, t, V_TO);

% Corresponding takeoff distance
x_TO = interp1(t, dist, t_TO);

figure();
yyaxis left
plot(t, Vel, "c-",'LineWidth', 1.5)
ylabel('Velocity (m/s)')
yline(V_TO, "k--")
grid on

yyaxis right
plot(t, dist, "m-", 'LineWidth', 1.5);
ylabel('Distance (m)')

xlabel('Time (s)')
title('Takeoff Performance')

xline(t_TO, 'k--', 'LineWidth', 1.5)

yyaxis right
plot(t, dist, "m-", 'LineWidth', 1.5);
ylabel('Distance (m)')

legend('Distance','Velocity','Location','best')

% cl_g is cl_0 and de is 0 at ground conditions
function dydt = takeoff_ode(t,y,weight_TO, rho, wing_area_total, mu_ground, CL_0, CD_0, K_w,Thrust_func,K_t,Sh_S,CL0_t,CLa_t) 

    % velocity
    dist = y(1);     % position
    Vel = y(2); 
    
    %% --- Thrust ---
    for i = 1:length(Vel)
        T(i) = Thrust_func(Vel(i),9989);
    end
    
    %% --- Lift ---
    L = 0.5*rho*Vel^2*wing_area_total*CL_0;
    
    %% --- Drag ---
    alpha_t = 0;
    CL_t = CL0_t + CLa_t * alpha_t;
    D = 0.5*rho*Vel.^2*wing_area_total * ...
        (CD_0 + K_w*CL_0^2 + Sh_S * K_t * CL_t.^2); % ask on Friday about this
    
    %% --- Friction ---
    Ff = mu_ground*(weight_TO - L);
    
    %% --- Acceleration ---
    m = weight_TO/9.81;
    a_long = (T - D - Ff)/m;
    
    %% --- ODE system ---
    dydt = zeros(2,1);
    dydt(1) = Vel;        % dx/dt
    dydt(2) = a_long;   % dV/dt


end

