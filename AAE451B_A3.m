% A3 assignment 
clc
clear
close all

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
t_CL = 13.29; % calculated in A3
etaP_CL = 0.55;
gamma = 10;
LD_CL = 10;

WB_W_CL = ((V_CL*t_CL)/(etaP_CL*etaM*etaESC*rho_B))*(((cosd(gamma))/(LD_CL)) + sind(gamma));

fprintf("Climb battery weight fraction: \n");
disp(WB_W_CL);

% constants - takeoff
t_TO = 3; % needs to be adjusted based on calculations
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


% accounting for efficiencies 

etaUse = 0.8; 
etaTemp = 0.7;

WB_W_eff= WB_W_total/(etaP_LF*etaM*etaESC*etaUse*etaTemp);

fprintf("Total battery weight fraction accounting for efficiency: \n");
disp(WB_W_eff);

% finding the empty weight
W = linspace(0, 100, 1000);

Wb_plus__wp = WB_W_eff.*W + 5.1012;
W_minus_we = 0.1582.*W + 2.9575;

W_int = (2.9575 - 5.1012) / (WB_W_eff - 0.1582);
Y_int = WB_W_eff*W_int + 5.1012;

figure(1);
plot(W, Wb_plus__wp, "m-", LineWidth=2);
hold on
plot(W, W_minus_we, "b-", LineWidth=2);
plot(W_int, Y_int, 'go', 'MarkerSize', 8, 'LineWidth', 2)
text(W_int, Y_int, sprintf('  (%.2f, %.2f)', W_int, Y_int))
xlabel('Weight (Newtons)');
ylabel('We-W and Wb+Wp (Newtons)');
title('Weight estimation');
legend('Wb + wp', 'W - we');
grid on;


% empty weight fraction
w_e = W_int - 0.1582*W_int - 2.9575;
we_W = w_e/W_int;
fprintf("Empty weight fraction: \n");
disp(we_W);

w_b = WB_W_eff*W_int;
w_p = 5.1012; % weight of payload in newtons

figure(2);
weights = [w_e w_b w_p];
names = ["Empty Weight", "Battery Weight", "Payload Weight"];
piechart(weights, names);
title("Weight Pie Chart");

% energy calculations
E_LF = (WB_W_LF*W_int) * rho_B;
E_TF = (WB_W_TF*W_int) * rho_B;
E_CL = (WB_W_CL*W_int) * rho_B;
E_TO = (WB_W_TO*W_int) * rho_B;
E_WU = (WB_W_WU*W_int) * rho_B;
E_total = (E_LF + E_TF + E_CL + E_TO + E_WU)/(etaTemp*etaUse);

fprintf("Total Energy = \n");
disp(E_total);



