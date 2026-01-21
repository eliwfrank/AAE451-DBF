% A3 assignment 

g = 9.81;
% constants - level flight
V_CR = 18; 
t_LF = 135;
LD_LF = 12.3; 
etaP_LF = 0.7;
etaM = 0.8;
etaESC = 0.95;
rho_B = 527000/g;

WB_W_LF = (V_CR*t_LF)/(LD_LF*etaP_LF*etaM*etaESC*rho_B);

fprintf("Level flight battery weight fraction: \n");
disp(WB_W_LF);

% constants - turning flight
V_M = 17;
t_TF = 61;
n = 2;
LD_TF = 11;
etaP_M = 0.65;

WB_W_TF = (V_M*t_TF*n)/(LD_TF*etaP_M*etaM*etaESC*rho_B);

fprintf("Turning Flight battery weight fraction: \n");
disp(WB_W_TF);

% constants - climb
V_CL = 13;
t_CL = 13.29; % needs to be adjusted based on calculations
etaP_CL = 0.55;
gamma = 10;
LD_CL = 10;

WB_W_CL = ((V_CL*t_CL)/(etaP_CL*etaM*etaESC*rho_B))*(((cosd(gamma))/(LD_CL)) + sind(gamma));

fprintf("Climb battery weight fraction: \n");
disp(WB_W_CL);

% constants - takeoff
t_TO = 15; % needs to be adjusted based on calculations
W_P = 0.1049;

WB_W_TO = t_TO/(etaM*etaESC*rho_B*W_P);

fprintf("Takeoff battery weight fraction: \n");
disp(WB_W_TO);

% constants - warmup 
N = 10;

WB_W_WU = N*(WB_W_TO);

fprintf("Warmup battery weight fraction: \n");
disp(WB_W_WU);

% total weight fraction calculation

WB_W_total = WB_W_TO + WB_W_WU + WB_W_CL + WB_W_LF + WB_W_TF;

fprintf("Total battery weight fraction: \n");
disp(WB_W_total);
