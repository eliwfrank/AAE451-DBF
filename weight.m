% A3 assignment, weight and weight fraction
fprintf("\nWEIGHT ------------------------\n")

% level flight/cruise
WB_W_LF = (V_CR*t_CR)/(LD_CR*etaP_CR*etaM*etaESC*rho_B);

% turning flight/maneuver
WB_W_TF = (V_M*t_M*n)/(LD_M*etaP_M*etaM*etaESC*rho_B);

% climb
WB_W_CL = ((V_CL*t_CL)/(etaP_CL*etaM*etaESC*rho_B))*(((cosd(gamma))/(LD_CL)) + sind(gamma));

% takeoff
WB_W_TO = t_TO/(etaM*etaESC*rho_B*W_P);

% warmup 
WB_W_WU = N*(WB_W_TO); % takeoff battery weight fraction

% total weight fraction calculation
WB_W_total = WB_W_TO + WB_W_WU + WB_W_CL + WB_W_LF + WB_W_TF;

% accounting for efficiencies 
WB_W_eff= WB_W_total/(etaP_CR*etaM*etaESC*etaUse*etaTemp);

% finding the empty weight
W = linspace(0, 100, 1000);
Wb_plus__wp = WB_W_eff.*W + 5.1012;
W_minus_we = 0.1582.*W + 2.9575;

weight_TO = (2.9575 - 5.1012) / (WB_W_eff - 0.1582);
Y_int = WB_W_eff*weight_TO + 5.1012;

figure();
plot(W, Wb_plus__wp, "m-", LineWidth=2);
hold on
plot(W, W_minus_we, "b-", LineWidth=2);
plot(weight_TO, Y_int, 'go', 'MarkerSize', 8, 'LineWidth', 2)
text(weight_TO, Y_int, sprintf('  (%.2f, %.2f)', weight_TO, Y_int))
xlabel('Weight (Newtons)');
ylabel('We-W and Wb+Wp (Newtons)');
title('Weight estimation');
legend('Wb + wp', 'W - we');
grid on;


% empty weight fraction
w_e = weight_TO - 0.1582*weight_TO - 2.9575;
we_W = w_e/weight_TO;

w_b = WB_W_eff*weight_TO;

figure();
weights = [w_e w_b w_p];
names = ["Empty Weight", "Battery Weight", "Payload Weight"];
piechart(weights, names);
title("Weight Pie Chart");

% energy calculations
E_LF = (WB_W_LF*weight_TO) * rho_B;
E_TF = (WB_W_TF*weight_TO) * rho_B;
E_CL = (WB_W_CL*weight_TO) * rho_B;
E_TO = (WB_W_TO*weight_TO) * rho_B;
E_WU = (WB_W_WU*weight_TO) * rho_B;
E_total = (E_LF + E_TF + E_CL + E_TO + E_WU)/(etaTemp*etaUse);

fprintf("Total Weight is: %f [N] or %.4f [kg]\n",weight_TO,weight_TO/g)


