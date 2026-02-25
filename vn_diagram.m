% v-n diagram
velocity_array_knots = (velocity_array).*1.944;
weight_TO_pounds = weight_TO/4.44822;
n_plus = 3.5;
n_minus = -1.5;

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

%% Key Speeds
V_A_plus  = V_int_plus;      % maneuvering speed (+)
V_A_minus = V_int_minus;     % maneuvering speed (-)

% Gust lines
mu_gust = (2*(W_S))/(rho*c*CL_a*g);
K_gust = (0.88*mu_gust)/(5.3 + mu_gust);
Ude_max = 50*0.305; % converting ft/s to m/s
Ude_25 = 25*0.305; % converting ft/s to m/s
n_gust_50 = 1 + (K_gust.*Ude_max.*velocity_array.*CL_a)./(2*(W_S));
n_gust_25 = 1 + (K_gust.*Ude_25.*velocity_array.*CL_a)./(2*(W_S));
n_gust_neg50 = 1 - (K_gust.*Ude_max.*velocity_array.*CL_a)./(2*(W_S));
n_gust_neg25 = 1 - (K_gust.*Ude_25.*velocity_array.*CL_a)./(2*(W_S));

% building the v-n envelope

V_gust_top = ((n_plus - 1) * 2 * (W_S)) / ...
    (K_gust * Ude_max * CL_a);

V_gust_bot = ((1 - n_minus) * 2 * (W_S)) / ...
    (K_gust * Ude_max * CL_a);

V_gust_top_knots = V_gust_top * 1.944;
V_gust_bot_knots = V_gust_bot * 1.944;

V_corner_top = min([V_A_plus, V_gust_top_knots, V_D_knots]);
V_corner_bot = min([V_A_minus, V_gust_bot_knots, V_D_knots]);

idx_pos = velocity_array_knots <= V_corner_top;
V_pos_curve = velocity_array_knots(idx_pos);
n_pos_curve = n_v_plus(idx_pos);

idx_neg = velocity_array_knots <= V_corner_bot;
V_neg_curve = velocity_array_knots(idx_neg);
n_neg_curve = n_v_minus(idx_neg);

V_top = linspace(V_corner_top, V_D_knots, 50);
n_top = n_plus * ones(size(V_top));

V_vert = V_D_knots * ones(1,50);
n_vert = linspace(n_plus, n_minus, 50);

V_bot = linspace(V_D_knots, V_corner_bot, 50);
n_bot = n_minus * ones(size(V_bot));

V_envelope = [0,V_pos_curve,V_top,V_vert,V_bot,fliplr(V_neg_curve), 0];
n_envelope = [0,n_pos_curve,n_top,n_vert,n_bot,fliplr(n_neg_curve),0];


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
% V_A (positive side)
text(V_A_plus, n_plus, '  V_A', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold');

% V_C (Cruise speed)
text(V_CR*1.944, n_plus*0.9, '  V_C', ...
    'Rotation',90, ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold');

% V_D (Dive speed)
text(V_D_knots, n_plus*0.9, '  V_D', ...
    'Rotation',90, ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold');

% N_max
text(V_D_knots*0.02, n_plus, '  n_{max}', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold');

% n_min
text(V_D_knots*0.02, n_minus, '  n_{min}', ...
    'VerticalAlignment','top', ...
    'FontWeight','bold');
title("V-N diagram");
xlabel("Velocity (Knots)");
ylabel('Load Factor (n)');
grid on