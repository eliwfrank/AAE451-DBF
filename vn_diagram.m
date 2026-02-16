% v-n diagram
n_plus = 2.1 + (24000/(W_int + 10000));
n_minus = (-0.4)*n_plus;

n_v_plus = (rho/2)*(1.688^2).*(((velocity_array.^2)*(CL_max))/(W_S));
n_v_minus = (rho/2)*(1.688^2).*(((velocity_array.^2)*(-1*CL_max))/(W_S));
V_D = V_CR/0.8; % estimated dive speed

V_int_plus = sqrt((2*W_S*n_plus)/(rho*(1.688^2)*CL_max));
V_int_minus = sqrt((2*W_S*abs(n_minus))/(rho*(1.688^2)*CL_max));

n_v_plus_trim = n_v_plus;
n_v_plus_trim(velocity_array > V_int_plus) = NaN;

n_v_minus_trim = n_v_minus;
n_v_minus_trim(velocity_array > V_int_minus) = NaN;

figure();
yline(n_plus, "b--", linewidth = 1.5);
hold on
yline(n_minus, "b--", linewidth = 1.5);
xline(V_CR, "b--", LineWidth=1.5);
xline(V_D, "b--", LineWidth=1.5);
plot(velocity_array, n_v_plus_trim, 'g-', linewidth = 1.5);
plot(velocity_array, n_v_minus_trim, 'g-', linewidth = 1.5);
plot(V_int_plus, n_plus, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
plot(V_int_minus, n_minus, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
grid on