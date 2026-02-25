%% ============================================================
%           FULL AIRCRAFT STRUCTURAL MARGIN OF SAFETY
% ============================================================

fprintf('\n====================================================\n')
fprintf('           STRUCTURAL MARGIN OF SAFETY\n')
fprintf('====================================================\n\n')

%% ------------------------------------------------------------
% MATERIAL ALLOWABLES  (EDIT TO MATCH YOUR MATERIAL DATA)
%% ------------------------------------------------------------

sigma_allow_carbon = 600e6;     % Pa
sigma_allow_al     = 276e6;     % Pa (6061-T6)
sigma_allow_balsa  = 20e6;      % Pa (conservative)

tau_allow_carbon   = 80e6;      % Pa
tau_allow_balsa    = 8e6;       % Pa

%% ============================================================
%% WING SPAR BENDING
%% ============================================================

sigma_wing_carbon = max(abs(stress_wc));
sigma_wing_balsa  = max(abs(stress_wb));

MS_wing_carbon = sigma_allow_carbon / sigma_wing_carbon - 1;
MS_wing_balsa  = sigma_allow_balsa  / sigma_wing_balsa  - 1;

fprintf('Wing Spar (Carbon) MS: %.3f\n', MS_wing_carbon)
fprintf('Wing Spar (Balsa)  MS: %.3f\n\n', MS_wing_balsa)

%% ============================================================
%% WING SHEAR (MAX ROOT SHEAR)
%% ============================================================

V_root = abs(V(1));
A_spar = pi*(Do^2 - Di^2)/4;

tau_wing = V_root / A_spar;

MS_wing_shear = tau_allow_carbon / tau_wing - 1;

fprintf('Wing Spar Shear MS: %.3f\n\n', MS_wing_shear)

%% ============================================================
%% WING TORSION (FROM SCRIPT 2)
%% ============================================================

if exist('T','var')

    J = (pi/32)*(Do^4 - Di^4);
    tau_torsion = T(1)*(Do/2) / J;

    MS_wing_torsion = tau_allow_carbon / tau_torsion - 1;

    fprintf('Wing Torsion MS: %.3f\n\n', MS_wing_torsion)

else
    fprintf('Wing Torsion MS: (No torsion input found)\n\n')
end

%% ============================================================
%% HORIZONTAL TAIL SPAR
%% ============================================================

sigma_tail_carbon = max(abs(stress_tc));
sigma_tail_balsa  = max(abs(stress_tb));

MS_tail_carbon = sigma_allow_carbon / sigma_tail_carbon - 1;
MS_tail_balsa  = sigma_allow_balsa  / sigma_tail_balsa  - 1;

fprintf('Horizontal Tail Spar (Carbon) MS: %.3f\n', MS_tail_carbon)
fprintf('Horizontal Tail Spar (Balsa)  MS: %.3f\n\n', MS_tail_balsa)

%% ============================================================
%% TAIL BOOM (ELEVATOR CASE)
%% ============================================================

sigma_boom = max(abs(stress_boom));

MS_boom_bending = sigma_allow_carbon / sigma_boom - 1;

fprintf('Tail Boom Bending (Elevator) MS: %.3f\n', MS_boom_bending)

%% ------------------------------------------------------------
%% TAIL BOOM (RUDDER CASE)
%% ------------------------------------------------------------

sigma_boom_rudder = max(abs(stress_booma));

MS_boom_rudder = sigma_allow_carbon / sigma_boom_rudder - 1;

fprintf('Tail Boom Bending (Rudder)   MS: %.3f\n', MS_boom_rudder)

%% ------------------------------------------------------------
%% TAIL BOOM TORSION
%% ------------------------------------------------------------

if exist('T_boom','var')

    tau_boom = max(abs(T_boom))*(Do_bcarbon/2) / I_boom;

    MS_boom_torsion = tau_allow_carbon / tau_boom - 1;

    fprintf('Tail Boom Torsion MS: %.3f\n', MS_boom_torsion)
end

%% ------------------------------------------------------------
%% TAIL BOOM COLUMN BUCKLING
%% ------------------------------------------------------------

L_boom = L;   % length from your script
E = 150e9;

P_cr = (pi^2*E*I_boom)/(L_boom^2);
P_actual = max(abs(Fz));   % elevator force

MS_boom_buckling = P_cr/P_actual - 1;

fprintf('Tail Boom Buckling MS: %.3f\n\n', MS_boom_buckling)

%% ============================================================
%% SUMMARY
%% ============================================================

MS_all = [
    MS_wing_carbon
    MS_wing_balsa
    MS_wing_shear
    MS_tail_carbon
    MS_tail_balsa
    MS_boom_bending
    MS_boom_rudder
    MS_boom_buckling
];

MS_min = min(MS_all);

fprintf('----------------------------------------------------\n')
fprintf('MINIMUM AIRCRAFT MARGIN OF SAFETY: %.3f\n', MS_min)
fprintf('----------------------------------------------------\n\n')

fprintf("Wing bending stress = %.5f MPa\n", sigma_wing_carbon/1e6)
fprintf("Wing shear stress   = %.5f MPa\n", tau_wing/1e6)
fprintf("Wing torsion stress = %.5f MPa\n", tau_torsion/1e6)