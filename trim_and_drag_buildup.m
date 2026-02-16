CL_req = linspace(-1,2,1000);

x_cg_chosen = (x_cg_aft + x_cg_forward) * 0.5;
VH = Sh_S * lt_c;

CL_det = CLa_t/pi() * (acos(1-2*E) + 2*sqrt(E*(1-E)));
CM0 = Cm0_w;
CL_de = Sh_S * CL_det;

CM_de = CL_det * Sh_S * (x_cg_chosen - x_ac) - CL_det * VH;

CM_a = CL_a * (x_cg_chosen - x_ac) - CLa_t * (1-de_da) * VH;

a_trim = rad2deg((CM0 * CL_de + CM_de * (CL_req - CL_0)) ./ (CL_a * CM_de - CL_de * CM_a));
de_trim = rad2deg(-((CM0 * CL_a + CM_a*(CL_req - CL_0)) ./ (CL_a * CM_de - CL_de * CM_a)));

figure()
plot(CL_req,a_trim)
xlabel("CL_{trim}")
ylabel("\alpha_{trim} [deg]")
title("Trim Angle of Attack vs. Trim Lift Coefficient")
grid on

figure()
plot(CL_req, de_trim)
xlabel("CL_{trim}")
ylabel("\delta_{trim} [deg]")
title("Trim Elevator Deflection vs. Trim Lift Coefficient")
grid on

Re = rho * V_CR * c / mu_air; %Reynolds number with chord as ref length

%% FUSELAGE
l_f = D_f* lambda_f; % fuselage diameter
FF_f = 0.9 + (5/lambda_f^1.5) + lambda_f/400; %fuselage form factor
S_wet_f = 1.02*((2*c*D_f) + 4*(c*l_f));
Q_f = 1.2;

% Re_f = rho * V_CR * l_f/mu_air; %Reynolds number with fuselage length as ref length
C_f_f = 0.455 / (log10(Re))^2.58;

%% WING
Lambda_m = 0;
M = 0.0524781; %flight mach number at 18m/s
FF_w = (1 + 0.6/(x_c_max)*(t_c) + 100*(t_c)^4) * (1.34*M^0.18*(cos(Lambda_m))^0.28); %wing form factor
S_wet_w = wing_area_total * (1.02)*2; %wetted area of wing
Q_w = 1;
C_f_w = 0.455 / (log10(Re))^2.58;

%% TAIL
FF_ht = (1 + 0.6/(x_c_max)*(t_c_tail) + 100*(t_c_tail)^4) * (1.34*M^0.18*(cos(Lambda_m))^0.28); %horizontal tail form factor
S_wet_ht = tail_area_h * 2 * (1.02); %wetted area of horizontal tail

FF_vt = (1 + 0.6/(x_c_max)*(t_c_tail) + 100*(t_c_tail)^4) * (1.34*M^0.18*(cos(Lambda_m_vtail))^0.28); %vertical tail form factor
S_wet_vt = tail_area_v * 2 * (1.02); %wetted area of vertical tail

% Re_t = rho * V_CR * c_tail_h/mu_air; %Reynolds number with chord as ref length
C_f_t = 0.455 / (log10(Re))^2.58;

Q_t = 1.05;

%% NACELLES

%we don't have nacelles yet

%% LANDING GEAR DRAG
CD_S_main = 0.18; % drag coefficient based on type of tires
CD_S_nose = 0.18; % drag coefficient based on type of tires

wheel_diameter = 0.06;
wheel_width = 0.003;

%% TOTAL PARASITE DRAG
% Re = rho * V_CR * c;
C_f = 0.455 / (log10(Re))^2.58;
S_ref = wing_area_total;
C_D_0_f = (FF_f * Q_f * C_f_f * S_wet_f) / S_ref;
C_D_0_w = (FF_w * Q_w * C_f_w * S_wet_w) / S_ref;
C_D_0_ht = (FF_ht * Q_t * C_f_t * S_wet_ht) / S_ref;
C_D_0_lg_nose = 1/2 * (wheel_diameter*wheel_width) * CD_S_nose / S_ref; % Landing gear
C_D_0_lg_main = 1/2 * (wheel_diameter*wheel_width) * CD_S_nose / S_ref; % Landing gear
C_D_0_vt = (FF_vt * Q_t * C_f_t * S_wet_vt) / S_ref;

%TOTAL PARASITIC DRAG
CD_0 = ((C_D_0_f + C_D_0_w + C_D_0_ht + C_D_0_vt + C_D_0_lg_main + C_D_0_lg_nose + 0.0027)*1.1)*1.25; 
fprintf("\nDRAG BUILD-UP ------------------------\n")

fprintf("Fuselage CD0: %.4f", C_D_0_f)
fprintf("\nWing CD0: %.4f", C_D_0_w)
fprintf("\nTail V CD0: %.4f", C_D_0_vt)
fprintf("\nTail H CD0: %.4f", C_D_0_ht)

fprintf("\nTotal CD0: %.4f\n", CD_0)

CL_w = CL0_w + CLa_w *  deg2rad(a_trim);

CL_t_trim = CL0_t + CLa_t *  deg2rad(a_trim) + CL_de * deg2rad(de_trim);
CL_t_clean = CL0_t + CLa_t *  deg2rad(a_trim);

K_w = 1./(pi() * e * AR_wing);
K_t = 1./(pi() * e * AR_tail);

CD_clean = CD_0 + K_w .* CL_w.^2 + Sh_S .* CL_t_clean.^2 .* K_t;
CD_trim =  CD_0 + K_w .* CL_w.^2 +  CL_t_trim.^2 * Sh_S .* K_t;

CL_trim = CL_w + CL_t_trim * Sh_S;
CL_clean = CL_w + CL_t_clean * Sh_S;

figure()
plot(CL_trim,CD_trim);
hold on
plot(CL_clean,CD_clean);
grid on
ylabel("CD")
xlabel("CL")
title("Drag Polar")
legend("Trimmed", "Clean")