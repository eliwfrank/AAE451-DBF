prop_data = readmatrix("12x6E.xlsx");

RPM = prop_data(:,1);
prop_vel = prop_data(:,2) /  2.237; % in m/s
prop_thrust = prop_data(:,12);

prop_vel = prop_vel(:);
prop_RPM = RPM(:);
prop_thrust = prop_thrust(:);

valid = isfinite(prop_vel) & isfinite(prop_RPM) & isfinite(prop_thrust);

prop_vel = prop_vel(valid);
prop_RPM = prop_RPM(valid);
prop_thrust = prop_thrust(valid);

Thrust_func = scatteredInterpolant(prop_vel, prop_RPM, prop_thrust);

[vel_q, RPM_q] = meshgrid( ...
    linspace(min(prop_vel), max(prop_vel), 50), ...
    linspace(min(prop_RPM), max(prop_RPM), 50));

zq = Thrust_func(vel_q, RPM_q);

figure()
scatter3(prop_vel, prop_RPM, prop_thrust, 15, prop_thrust, 'filled')
hold on
xlabel("Velocity (m/s)")
ylabel("RPM")
set(gca, 'YDir', 'reverse')
zlabel("Thrust (N)")
colormap('parula')
colorbar

surf(vel_q, RPM_q, zq, 'FaceAlpha', 0.5)
shading interp
colorbar
title("Multivariable Fit for Thrust vs. RPM and Velocity for APC 12x6E")