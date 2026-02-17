% v-n diagram
velocity_array_knots = (velocity_array).*1.944;
n_plus = 2.1 + (24000/(weight_TO + 10000));
n_minus = (-0.4)*n_plus;

n_v_plus = (rho/2)*(1.688^2).*(((velocity_array_knots.^2)*(CL_max))/(W_S)); 
% allegedly cl_max + and - should be different? idk we need to look at it
n_v_minus = (rho/2)*(1.688^2).*(((velocity_array_knots.^2)*(-1*CL_max))/(W_S));
V_D_knots = (V_CR/0.8)*1.944; % estimated dive speed in knots

V_int_plus = sqrt((2*W_S*n_plus)/(rho*(1.688^2)*CL_max));
V_int_minus = sqrt((2*W_S*abs(n_minus))/(rho*(1.688^2)*CL_max));

n_v_plus_trim = n_v_plus;
n_v_plus_trim(velocity_array_knots > V_int_plus) = NaN;

n_v_minus_trim = n_v_minus;
n_v_minus_trim(velocity_array_knots > V_int_minus) = NaN;

V_1g = sqrt((2*W_S*1)/(rho*(1.688^2)*CL_max));

% Envelope construction

% Positive stall curve
idx_pos = velocity_array_knots >= V_1g & velocity_array_knots <= V_int_plus;
V_pos_curve = velocity_array_knots(idx_pos);
n_pos_curve = n_v_plus(idx_pos);

% Negative stall curve
idx_neg = velocity_array_knots >= V_1g & velocity_array_knots <= V_int_minus;
V_neg_curve = velocity_array_knots(idx_neg);
n_neg_curve = n_v_minus(idx_neg);

% Intermediate segments
V_top = linspace(V_int_plus, V_D_knots, 50);
n_top = n_plus * ones(size(V_top));

V_vert = V_D_knots * ones(1,50);
n_vert = linspace(n_plus, n_minus, 50);

V_bot = linspace(V_D_knots, V_int_minus, 50);
n_bot = n_minus * ones(size(V_bot));

V_envelope = [V_1g, V_pos_curve, V_top, V_vert,V_bot, fliplr(V_neg_curve), ...
    V_1g];

n_envelope = [-1, n_pos_curve, n_top, n_vert, n_bot, fliplr(n_neg_curve),1];

% Gust lines
mu_gust = (2*(W_S))/(rho*c*CL_a*g);
K_gust = (0.88*mu_gust)/(5.3 + mu_gust);
Ude_max = 50*0.305; % converting ft/s to m/s
Ude_25 = 25*0.305; % converting ft/s to m/s
n_gust_50 = 1 + (K_gust.*Ude_max.*velocity_array.*CL_a)./(2*(W_S));
n_gust_25 = 1 + (K_gust.*Ude_25.*velocity_array.*CL_a)./(2*(W_S));
n_gust_neg50 = 1 - (K_gust.*Ude_max.*velocity_array.*CL_a)./(2*(W_S));
n_gust_neg25 = 1 - (K_gust.*Ude_25.*velocity_array.*CL_a)./(2*(W_S));

figure();
yline(n_plus, "b--", linewidth = 1.5);
hold on
yline(n_minus, "b--", linewidth = 1.5);
plot([0 V_1g], [1 1], 'b--', 'LineWidth', 1);
plot([0 V_1g], [-1 -1], 'b--', 'LineWidth', 1);
xline(V_CR*1.944, "b--", LineWidth=1.5);
xline(V_D_knots, "b--", LineWidth=1.5);
plot(velocity_array_knots, n_v_plus_trim, 'g-', linewidth = 1.5);
plot(velocity_array_knots, n_v_minus_trim, 'g-', linewidth = 1.5);
plot(V_int_plus, n_plus, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
plot(V_int_minus, n_minus, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
plot(V_envelope, n_envelope,'Color', [0 0.9 0.9],'LineWidth', 6);
plot(velocity_array_knots, n_gust_50, "m--", linewidth = 1.5);
plot(velocity_array_knots, n_gust_25, "m--", linewidth = 1.5);
plot(velocity_array_knots, n_gust_neg25, "m--", linewidth = 1.5);
plot(velocity_array_knots, n_gust_neg50, "m--", linewidth = 1.5);
title("V-N diagram");
xlabel("Velocity (Knots)");
ylabel('Load Factor (n)');
grid on