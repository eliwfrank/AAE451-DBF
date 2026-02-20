% def calculate_shrenk_lift_distribution(y, semi_span, chord, lift_coefficient, air_density, assumed_velocity, load_factor):
% """
% Calculates the lifting loads along the wing semi-span using a simplified approximation.
% Args:
% y: Points along the semi-span in meters.
% semi_span: Wing semi-span in meters.
% chord: Wing chord in meters.
% lift_coefficient: Lift coefficient (dimensionless).
% air_density: Air density in kg/m^3.
% assumed_velocity: Assumed flight velocity for load calculation (m/s).
% Returns:
% An array of lifting loads along the semi-span.
% """
% # Calculate dynamic pressure (q)
% q = 0.5 * air_density * (assumed_velocity)**2 * load_factor
% # Calculate lift distribution using a simplified Shrenk-like approximation
% # This formula gives a semi-elliptical like distribution, which is typical for
% # wings and often used as a simplification.
% lift_distribution = lift_coefficient * q * chord * np.sqrt(1 - (y / semi_span)**2

y = linspace(0,b/2,1000);
q = (1/2)*rho*(V_CR^2)*n_top(1); % Max Load for Pitch at Cruise
L_dist = CL_max_w*q*c*sqrt(1-(y/(b/2)).^2); % CL_max for wing at AoA Stall

% --------- SHEAR + BENDING (tip -> root integration) ----------
% Tip boundary conditions: V(L)=0, M(L)=0
yR = flip(y);
wR = flip(L_dist);

V_R = cumtrapz(yR, wR);      % shear distribution from tip inward (N)
M_R = cumtrapz(yR, V_R);     % bending moment distribution (N*m)

V = -flip(V_R);
M_bend = flip(M_R);
T = V*(0.25*c);

% Optional sanity checks (root values)
V_root_check = trapz(y, L_dist);
M_root_check = trapz(y, y.*L_dist);
check_root = 0.5*weight_TO*n_top(1);
fprintf('V_root = %.3f N (check %.3f)\n', V(1), V_root_check);
fprintf('V_root check = %.3f N \n', check_root);
fprintf('M_root = %.3f N*m (check %.3f)\n', M_bend(1), M_root_check);

% --------- PLOTS ----------
figure();
plot(y, L_dist, 'r-')
xlabel('Wing Semi-span (m)')
ylabel('Lift Distribution (N/m)')
title('Schrenk Lift Approximation')
grid on

figure();
plot(y, V, 'b-')
xlabel('Wing Semi-span (m)')
ylabel('V (N)')
title('Shear Force Distribution')
grid on

figure();
plot(y, M_bend, 'k-')
xlabel('Wing Semi-span (m)')
ylabel('M_{bending} (N*m)')
title('Bending Moment Distribution')
grid on

figure()
plot(y, T, "g-");
xlabel("wing semi-span (m)");
ylabel("M_{torsion} (N*m)");
title("Torsion Moment Distribution at 0.25*c");
grid on;