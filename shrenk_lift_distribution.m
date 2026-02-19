function lift_distribution = shrenk_lift_distribution( ...
            y, semi_span, chord, lift_coefficient, air_density, assumed_velocity, load_factor)

% CALCULATE_SHRENK_LIFT_DISTRIBUTION
% Calculates the lifting loads along the wing semi-span using a simplified approximation.
%
% Inputs:
% y                 : Points along the semi-span (m)
% semi_span         : Wing semi-span (m)
% chord             : Wing chord (m)
% lift_coefficient : Lift coefficient (dimensionless)
% air_density      : Air density (kg/m^3)
% assumed_velocity : Flight velocity (m/s)
% load_factor      : Load factor (n)
%
% Output:
% lift_distribution : Lifting load distribution along semi-span

% ---- Dynamic pressure (q) ----
q = 0.5 * air_density * assumed_velocity^2 * load_factor;

% ---- Shrenk-like semi-elliptical lift distribution ----
lift_distribution = lift_coefficient .* q .* chord .* ...
                    sqrt(1 - (y ./ semi_span).^2);

end
