% def calculate_shrenk_lift_distribution(y, semi_span, chord, lift_coefficient, air_density, assumed_velocity, load_factor):
% """
% Calculates the lifting loads along the wing semi-span using a simplified approximation.
% Args:
% y: Points along the semi-span in meters.
% semi_span: Wing semi-span in meters.
% chord: Wing chord in meters.
% lift_coefficient: Lift coefficient (dimensionless).
% air_density: Air density in kg/m^3.
% assumed_velocity: Assumed flight velocity for load calculation (m/s).
% Returns:
% An array of lifting loads along the semi-span.
% """
% # Calculate dynamic pressure (q)
% q = 0.5 * air_density * (assumed_velocity)**2 * load_factor
% # Calculate lift distribution using a simplified Shrenk-like approximation
% # This formula gives a semi-elliptical like distribution, which is typical for
% # wings and often used as a simplification.
% lift_distribution = lift_coefficient * q * chord * np.sqrt(1 - (y / semi_span)**2

%% Wing

y = linspace(0,b/2,1000);
q = (1/2)*rho*(V_CR^2)*n_plus; % Max Load for Pitch at Cruise
L_dist = CL_max_w*q*c*sqrt(1-(y/(b/2)).^2); % CL_max for wing at AoA Stall

% --------- SHEAR + BENDING (tip -> root integration) ----------
% Tip boundary conditions: V(L)=0, M(L)=0
yR = flip(y);
wR = flip(L_dist);

V_R = cumtrapz(yR, wR);      % shear distribution from tip inward (N)
M_R = cumtrapz(yR, V_R);     % bending moment distribution (N*m)

V = -flip(V_R);
My = flip(M_R);

Mz = zeros(size(My));            % set if you have lateral loads too

% Optional sanity checks (root values)
V_root_check = trapz(y, L_dist);
check_root = 0.5*weight_TO*n_plus;
fprintf('V_root = %.3f N (check %.3f)\n', V(1), V_root_check);
fprintf('V_root check = %.3f N \n', check_root);

%% Material Properties

% ---- Section properties: hollow carbon pipe ----
Do = 0.015;                          % outer diameter (m)
Di = 0.014;                        % inner diameter (m)
I_carbon  = (pi/64)*(Do^4 - Di^4);      % Iy = Iz = I
E_carbon  = 150*10^9; % Carbon Fiber Young Modulus (Pa);           

% ---- Exact equations from your sheet, simplified for pipe ----
vpp_carbon = -(Mz)/(E_carbon*I_carbon);               % v''(x)
wpp_carbon =  (My)/(E_carbon*I_carbon);               % w''(x)

% ---- Integrate twice with clamped root BC: v(0)=v'(0)=0, w(0)=w'(0)=0 ----
vp_carbon = cumtrapz(y, vpp_carbon);           % v'(0)=0 automatically
v_carbon  = cumtrapz(y, vp_carbon);            % v(0)=0 automatically

wp_carbon = cumtrapz(y, wpp_carbon);           % w'(0)=0 automatically
w_carbon  = cumtrapz(y, wp_carbon);            % w(0)=0 automatically

stress_wc = (My*(Do/2))/I_carbon;

% Aluminum Hollow Tube

Do = 0.015;                      % outer diameter (m)
Di = 0.014;                      % inner diameter (m)
I_aluminum  = (pi/64)*(Do^4 - Di^4);      % Iy = Iz = I
E_aluminum  = 71.7*10^9; % Aluminum 6061 T-651 Young Modulus (Pa)         

% ---- Exact equations from your sheet, simplified for pipe ----
vpp_aluminum = -(Mz)/(E_aluminum*I_aluminum);               % v''(x)
wpp_aluminum =  (My)/(E_aluminum*I_aluminum);               % w''(x)

% ---- Integrate twice with clamped root BC: v(0)=v'(0)=0, w(0)=w'(0)=0 ----
vp_aluminum = cumtrapz(y, vpp_aluminum);           % v'(0)=0 automatically
v_aluminum  = cumtrapz(y, vp_aluminum);            % v(0)=0 automatically

wp_aluminum = cumtrapz(y, wpp_aluminum);           % w'(0)=0 automatically
w_aluminum  = cumtrapz(y, wp_aluminum);            % w(0)=0 automatically

% Balsa Solid Tube

Do = 0.015;                      % outer diameter (m)
Di = 0;                      % inner diameter (m)
I_balsa  = (pi/64)*(Do^4 - Di^4);      % Iy = Iz = I
E_balsa  = 3.12*10^9; % Balsa Wood Young Modulus (Pa)           

% ---- Exact equations from your sheet, simplified for pipe ----
vpp_balsa = -(Mz)/(E_balsa*I_balsa);               % v''(x)
wpp_balsa =  (My)/(E_balsa*I_balsa);               % w''(x)

% ---- Integrate twice with clamped root BC: v(0)=v'(0)=0, w(0)=w'(0)=0 ----
vp_balsa = cumtrapz(y, vpp_balsa);           % v'(0)=0 automatically
v_balsa  = cumtrapz(y, vp_balsa);            % v(0)=0 automatically

wp_balsa = cumtrapz(y, wpp_balsa);           % w'(0)=0 automatically
w_balsa  = cumtrapz(y, wp_balsa);            % w(0)=0 automatically

stress_wb = (My*(Do/2))/I_balsa;


%% Horizontal Stabilizer 

y_t = linspace(0,b_tail_h/2,1000);
q_t = (1/2)*rho*(V_CR^2)*n_plus; % Max Load for Pitch at Cruise

CL_t_el_max = CL0_t + CLa_t * deg2rad(aoa_tail_stall) + CL_de * deg2rad(-15); % -15 (deg) is max elevator deflection (pitch up)

L_dist_t = CL_t_el_max*q_t*c_tail_h*sqrt(1-(y_t/(b_tail_h/2)).^2); % CL_max for horizontal tail at AoA

% --------- SHEAR + BENDING (tip -> root integration) ----------
% Tip boundary conditions: V(L)=0, M(L)=0
yR_t = flip(y_t);
wR_t = flip(L_dist_t);

V_R_t = cumtrapz(yR_t, wR_t);      % shear distribution from tip inward (N)
M_R_t = cumtrapz(yR_t, V_R_t);     % bending moment distribution (N*m)

V_t = -flip(V_R_t);
My_t = flip(M_R_t);

Mz_t = zeros(size(My_t));            % set if you have lateral loads too

%% Material Properties

% ---- Section properties: hollow carbon pipe ----
Do_t = 0.005;                          % outer diameter (m)
Di_t = 0.004;                        % inner diameter (m)
I_carbon_t  = (pi/64)*(Do_t^4 - Di_t^4);      % Iy = Iz = I
E_carbon  = 150*10^9; % Carbon Fiber Young Modulus (Pa);           

% ---- Exact equations from your sheet, simplified for pipe ----
vpp_carbon_t = -(Mz_t)/(E_carbon*I_carbon_t);               % v''(x)
wpp_carbon_t =  (My_t)/(E_carbon*I_carbon_t);               % w''(x)

% ---- Integrate twice with clamped root BC: v(0)=v'(0)=0, w(0)=w'(0)=0 ----
vp_carbon_t = cumtrapz(y_t, vpp_carbon_t);           % v'(0)=0 automatically
v_carbon_t  = cumtrapz(y_t, vp_carbon_t);            % v(0)=0 automatically

wp_carbon_t = cumtrapz(y_t, wpp_carbon_t);           % w'(0)=0 automatically
w_carbon_t  = cumtrapz(y_t, wp_carbon_t);            % w(0)=0 automatically

stress_tc = (My_t*(Do/2))/I_carbon_t;

% Aluminum Hollow Tube

Do_t = 0.005;                      % outer diameter (m)
Di_t = 0.004;                      % inner diameter (m)
I_aluminum_t  = (pi/64)*(Do_t^4 - Di_t^4);      % Iy = Iz = I
E_aluminum  = 71.7*10^9; % Aluminum 6061 T-651 Young Modulus (Pa)         

% ---- Exact equations from your sheet, simplified for pipe ----
vpp_aluminum_t = -(Mz_t)/(E_aluminum*I_aluminum_t);               % v''(x)
wpp_aluminum_t =  (My_t)/(E_aluminum*I_aluminum_t);               % w''(x)

% ---- Integrate twice with clamped root BC: v(0)=v'(0)=0, w(0)=w'(0)=0 ----
vp_aluminum_t = cumtrapz(y_t, vpp_aluminum_t);           % v'(0)=0 automatically
v_aluminum_t  = cumtrapz(y_t, vp_aluminum_t);            % v(0)=0 automatically

wp_aluminum_t = cumtrapz(y_t, wpp_aluminum_t);           % w'(0)=0 automatically
w_aluminum_t  = cumtrapz(y_t, wp_aluminum_t);            % w(0)=0 automatically

% Balsa Solid Tube

Do_t = 0.005;                      % outer diameter (m)
Di_t = 0;                      % inner diameter (m)
I_balsa_t  = (pi/64)*(Do_t^4 - Di_t^4);      % Iy = Iz = I
E_balsa  = 3.12*10^9; % Balsa Wood Young Modulus (Pa)           

% ---- Exact equations from your sheet, simplified for pipe ----
vpp_balsa_t = -(Mz_t)/(E_balsa*I_balsa_t);               % v''(x)
wpp_balsa_t =  (My_t)/(E_balsa*I_balsa_t);               % w''(x)

% ---- Integrate twice with clamped root BC: v(0)=v'(0)=0, w(0)=w'(0)=0 ----
vp_balsa_t = cumtrapz(y_t, vpp_balsa_t);           % v'(0)=0 automatically
v_balsa_t  = cumtrapz(y_t, vp_balsa_t);            % v(0)=0 automatically

wp_balsa_t = cumtrapz(y_t, wpp_balsa_t);           % w'(0)=0 automatically
w_balsa_t  = cumtrapz(y_t, wp_balsa_t);            % w(0)=0 automatically

stress_tb = (My_t*(Do/2))/I_balsa_t;

%% Tail Boom (full elevator deflection/Symmetric)

L = 1.0; % tail boom length
x = linspace(0,L,1000).';

Fz = 0.5*rho*(V_CR^2)*S_wet_ht*CL_t_el_max;  % N (tail load at tip)
E  = E_carbon;  % Pa
Do_bcarbon = 0.03;  % m
Di_bcarbon = 0.02;  % m

I_boom = (pi/64)*(Do_bcarbon^4 - Di_bcarbon^4);

V_boom = Fz*ones(size(x));
M_boom = Fz*(L-x);

wpp_boom = M_boom/(E*I_boom);
wp_boom  = cumtrapz(x, wpp_boom);
w_boom   = cumtrapz(x, wp_boom);

stress_boom = (M_boom*(Do/2))/I_boom;

%% Tail Boom (full rudder deflection/ Assymetric)

CL_vt_rudder_max = CL0_t + CL_de * deg2rad(30); % yaw to the right

Fy = 0.5*rho*(V_CR^2)*S_wet_vt*CL_vt_rudder_max;  % N (tail load at tip)
E  = E_carbon;  % Pa
Do_bcarbon = 0.03;  % m
Di_bcarbon = 0.02;  % m

I_boom = (pi/64)*(Do_bcarbon^4 - Di_bcarbon^4);

V_yboom = Fy*ones(size(x));
M_yboom = Fy*(L-x);

vpp_boom = M_yboom/(E*I_boom);
vp_boom  = cumtrapz(x, vpp_boom);
v_boom   = cumtrapz(x, vp_boom);

stress_booma = (M_yboom*(Do/2))/I_boom;

T_boom = V_yboom*(0.42*b_tail_v); % Torsion on Boom generated by the distance between the tube center and Vert. Stab. Centroid Force Dist.
