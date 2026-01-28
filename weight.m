% A3 assignment, weight and weight fraction

% level flight
WB_W_LF = (V_CR*t_LF)/(LD_LF*etaP_LF*etaM*etaESC*rho_B);

% fprintf("Level flight battery weight fraction: \n");
% disp(WB_W_LF);

% constants - turning flight
WB_W_TF = (V_M*t_TF*n)/(LD_TF*etaP_M*etaM*etaESC*rho_B);

% fprintf("Turning Flight battery weight fraction: \n");
% disp(WB_W_TF);

% climb
WB_W_CL = ((V_CL*t_CL)/(etaP_CL*etaM*etaESC*rho_B))*(((cosd(gamma))/(LD_CL)) + sind(gamma));

% fprintf("Climb battery weight fraction: \n");
% disp(WB_W_CL);

% constants - takeoff
WB_W_TO = t_TO/(etaM*etaESC*rho_B*W_P);

% fprintf("Takeoff battery weight fraction: \n");
% disp(WB_W_TO);

% warmup 
WB_W_WU = N*(WB_W_TO); % takeoff battery weight fraction

% fprintf("Warmup battery weight fraction: \n");
% disp(WB_W_WU);

% total weight fraction calculation
WB_W_total = WB_W_TO + WB_W_WU + WB_W_CL + WB_W_LF + WB_W_TF;

% accounting for efficiencies 
WB_W_eff= WB_W_total/(etaP_LF*etaM*etaESC*etaUse*etaTemp);

% fprintf("Total battery weight fraction accounting for efficiency: \n");
% disp(WB_W_eff);

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
% fprintf("Empty weight fraction: \n");
% disp(we_W);

w_b = WB_W_eff*W_int;

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

% fprintf("Total Energy = \n");
% disp(E_total);

fprintf("Total Weight is: %f [N] \n%27.4d [kg]\n\n",W_int,W_int/g)


