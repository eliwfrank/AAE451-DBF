%% environment
rho = 1.20343; % Air density [kg/m^3]
g = 9.81; % [m/s^2]
mu = 0.04; % Ground friction coefficient

%% Velocities
V_s = 10; % Stall Velocity [m/s]
V_M = 17; % Maneuver Velocity [m/s]
V_CR = 18; % Cruise Velocity
V_CL = 13; % Climb Velocity
V_TO = 12; % TO velocity [m/s]

%% Coefficients
CL_max = 1.5; % Max Lift Coefficient
CD_0 = 0.03; % Parasitic Drag Coefficient
CD_G = .0875; % Ground Drag Coefficient
CL_R = 0.8876; % CL at rotation

%% Efficiencies
etaP_CR = 0.7; % Cruise/level flight Propellor Efficiency
etaP_TO = 0.5; % TO prop efficiency
etaP_CL = 0.5; % Climb Propellor efficiency
etaP_M = 0.65; % maneuvering propeller efficiency
etaM = 0.8; % motor efficiency
etaESC = 0.95; % ESC efficiency
etaUse = 0.8; % use efficiency
etaTemp = 0.8; % temperature efficiency

%% L/D
LD_max = 12;
LD_CR = 12.3; % lift to drag ratio in level flight
LD_M = 11; % lift to drag ratio in turning flight
LD_CL = 10; % lift to drag ratio in climb

%% phase times
t_CR = 135; % level flight time (s)
t_M = 61; % turning flight time (s)
t_CL = 13.29; % climb time
t_TO = 3; % take off time (s)

%% random other
gamma = 12; % climb angle [deg]
phi = 0.75; % Cruise Throttle
AR = 4; % the same for individual wing as both wings combined since using the same size wing
e = 0.85; % Oswald Efficiency Factor
q = 0.5 * rho * V_M^2;
S_TO = 18; % TO run [m]
rho_B = 527000/g; % batter energy density (J/N)
n = 1.3; % turning flight load factor
N = 10; % multiple of takeoff time, warmup
w_p = 5.1012; % weight of payload in newtons

%% Tail Sizing
SM = 0.2; % Static Margin
x_ac = 0.25; % Aerodynamic Center of the Wing
etaH = 1; % tail efficiency, "if it is not aligned with the wing the number is close to 1.0"
CLa_w = 5.7; % CLA of wing
AR_tail = 4;
de_da = CLa_w / (pi*AR);
selected_Sh_S = 0.12;
Cm0_w = -0.05; % from clark y
CLt_nose_up = -0.8; % [-0.5 to -0.8] max negative tail lift for nose up control
Cm_req_rot = 0.2; % [0.1 to 0.2] dependant on your configuration and the position of the main landing gear
CLa_t = 6.86 ; % tail CLa per radian
aoa_tail_stall = deg2rad(12.25); % Stall angle of attack
