%% environment
rho = 1.20343; % Air density [kg/m^3]
g = 9.81; % [m/s^2]
mu_ground = 0.04; % Ground friction coefficient
mu_air = 1.73e-5;
n = 1.2; % turning flight load factor

%% Velocities
V_si = 10; % Stall Velocity [m/s], recalculated in Design Point.
V_TOi = n * V_si; % TO velocity [m/s]

V_M = 17; % Maneuver Velocity [m/s]
V_CR = 18; % Cruise Velocity
V_CL = 13; % Climb Velocity

%% Coefficients
CD_0i = 0.0383; % Parasitic Drag Coefficient
CL_Ri = 1.0541; % CL at rotation - recalculated, 
CD_G = .0875; % Ground Drag Coefficient

%% Efficiencies
etaP_CRi = 0.6595; % Cruise/level flight Propellor Efficiency, calculated in prop system periodically update
etaP_TOi = 0.5024; % TO prop efficiency, calculated in prop system periodically update
etaP_CLi = 0.5490; % Climb Propellor efficiency, calculated in prop system periodically update
etaP_Mi = 0.6418; % maneuvering propeller efficiency, calculated in prop system periodically update

etaM = 0.8; % motor efficiency
etaESC = 0.95; % ESC efficiency
etaUse = 0.8; % use efficiency
etaTemp = 0.8; % temperature efficiency
etaH = 1; % tail efficiency, "if it is not aligned with the wing the number is close to 1.0"

%% phase times
t_CR = 135; % level flight time (s)
t_M = 61; % turning flight time (s)
t_CL = 13.29; % climb time
t_TO = 3; % take off time (s)

%% random other
gamma = 12; % climb angle [deg]
phi = 0.75; % Cruise Throttle
e = 0.85; % Oswald Efficiency Factor
q = 0.5 * rho * V_M^2;
S_TO = 22.5; % TO run [m]
rho_B = 527000/g; % batter energy density (J/N)
N = 10; % multiple of takeoff time, warmup
w_p = 5.1012; % weight of payload in newtons

%% Tail Sizing
SMi = 0.2; % Static Margin
x_ac = 0.25; % Aerodynamic Center of the Wing
Sh_S = 0.15;
CLt_nose_up = -0.8; % [-0.5 to -0.8] max negative tail lift for nose up control
Cm_req_rot = 0.25; % [0.1 to 0.2] dependant on your configuration and the position of the main landing gear
lt = 1; % Tail Arm
Lambda_m_vtail = .405; %sweep angle of vertical tail

%% Airfoil Characteristics
% Clark Y
x_c_max = 0.3;
t_c = 0.117;

cl0_w = 0.36; 
cm0_w = -0.05; % from clark y
cl_max_w = 1.52;
cla_w = 5.7; % CLA of wing

% Naca 0012
t_c_tail = 0.12;
aoa_tail_stall = deg2rad(12.25); % Stall angle of attack

cl0_t = 0; % Tail CL0
cl_max_t = 1.1101;
cla_t = 6.86 ; % tail CLa per radian

%% Wing Characteristics
AR_wing = 7; % the same for individual wing as both wings combined since using the same size wing
CL0_w = cl0_w; % Wing CL0
Cm0_w = cm0_w; % from clark y
CM0 = Cm0_w;
CLa_w = cla_w / (1 + cla_w / (pi() * AR_wing));
CL_max_w = 0.9 * cl_max_w;

K_w = 1./(pi() * e * AR_wing);

%% Tail Characteristics
AR_tail = 4;
CL0_t = cl0_t; % Tail CL0
CLa_t = cla_t / (1 + cla_t / (pi() * AR_wing)); % found in lecture notes
CL_max_t = 0.9 * cl_max_t;

de_da = CLa_w / (pi*AR_wing);
K_t = 1./(pi() * e * AR_tail);

%% Aircraft Coefficients
CL_max = CL_max_w + Sh_S * CL_max_t;
CL_a = CLa_w + Sh_S * CLa_t * (1 - de_da);
CL_0 = CL0_w + Sh_S * CL0_t;

%% Fuselage Design
D_f = 0.1; % Diameter of the fuselage
lambda_f = 5.1; %fineness/slenderness ratio, 5.1 is optimal per Sadraey

%% Control Surface
SE_Sh = 0.4; % Elevator surface Area / lifting surface (horzontal tail) area [0.2 - 0.4]
SA_S = 0.12; % Aileron surface Area / lifting surface (wing) area [0.03 - 0.12]
SR_Sv = 0.35; % Rudder surface Area / lifting surface (vertical tail) area [0.15 - 0.35]

CE_Ch = 0.4; % Elevator chord / lifting surface chord [0.2 - 0.4]
CA_C = 0.3; % Aileron chord / lifting surface chord [0.15 - 0.3]

%% L/D
LD_max = 1 / (2 * sqrt(K_w * CD_0i));
LD_CR = 12.3; % lift to drag ratio in level flight
LD_M = 11; % lift to drag ratio in turning flight
LD_CL = 10; % lift to drag ratio in climb
