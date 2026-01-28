%% Design Point
rho = 1.20343; % Air density [kg/m^3]
CD_0 = 0.03; % Parasitic Drag Coefficient
g = 9.81; % [m/s^2]

V_s = 10; % Stall Velocity [m/s]
CL_max = 1.5; % Max Lift Coefficient

eta_CR = 0.7; % Cruise Propellor Efficiency
phi = 0.75; % Cruise Throttle
V_CR = 18; % Cruise Velocity

eta_CL = 0.5; % Climb Propellor efficiency
V_CL = 13; % Climb Velocity
LD_max = 12;
gamma = 12; % climb angle [deg]

eta_M = 0.65;
V_M = 17; % Maneuver Velocity [m/s]
AR = 4; % the same for individual wing as both wings combined since using the same size wing
e = 0.85; % Oswald Efficiency Factor
n = 1.2; % Maneuvering load factor
q = 0.5 * rho * V_M^2;

CD_G = .0875; % Ground Drag Coefficient
S_TO = 18; % TO run [m]
mu = 0.04; % Ground friction coefficient
CL_R = 0.8876; % CL at rotation
eta_TO = 0.5; % TO prop efficiency
V_TO = 12; % TO velocity [m/s]

%% Sizing
t_LF = 135; % level flight time (s)
LD_LF = 12.3; % lift to drag ratio in level flight
etaP_LF = 0.7; % level flight propeller efficiency
etaM = 0.8; % motor efficiency
etaESC = 0.95; % ESC efficiency
rho_B = 527000/g; % batter energy density (J/N)

V_M = 17; % Maneuver velocity (m/s)
t_TF = 61; % turning flight time (s)
n = 1.3; % turning flight load factor
LD_TF = 11; % lift to drag ratio in turning flight
etaP_M = 0.65; % maneuvering propeller efficiency

t_CL = 13.29; % climb time
etaP_CL = 0.55; % climb propeller efficiency
LD_CL = 10; % lift to drag ratio in climb

t_TO = 3; % take off time (s)
N = 10; % multiple of takeoff time, warmup

etaUse = 0.8; % use efficiency
etaTemp = 0.8; % temperature efficiency

w_p = 5.1012; % weight of payload in newtons