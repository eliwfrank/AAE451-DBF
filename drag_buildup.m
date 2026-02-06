%% FUSELAGE
x_c_max = 0.3;
t_c = 0.117;

l_f = 2 * c; %generic fuselage length based on wing chord
lambda_f = 5.1; %fineness/slenderness ratio, 5.1 is optimal per Sadraey
D_f = l_f / lambda_f; % fuselage diameter
FF_f = 0.9 + (5/lambda_f^1.5) + lambda_f/400; %fuselage form factor
S_wet_f = pi*D_f*l_f * (1 - 2/lambda_f)^(2/3) * (1 + 1/lambda_f^2);
Q_f = 1;

Re_f = rho * V_CR * l_f; %Reynolds number with fuselage length as ref length
C_f_f = 0.455 / (log10(Re_f))^2.58;

%% WING
Lambda_m = 0;
M = 0.0524781; %flight mach number at 18m/s
FF_w = (1 + 0.6/(x_c_max)*(t_c) + 100*(t_c)^4) * (1.34*M^0.18*(cos(Lambda_m))^0.28); %wing form factor
S_wet_w = wing_area_total * 2 * 1.02; %wetted area of wing
Q_w = 1;
mu_air = 1.73e-5;
Re_w = rho * V_CR * c / mu_air; %Reynolds number with chord as ref length
C_f_w = 0.455 / (log10(Re_w))^2.58;

%% TAIL
t_c_tail = 0.12;

Lambda_m_vtail = .405; %sweep angle of vertical tail
FF_ht = (1 + 0.6/(x_c_max)*(t_c_tail) + 100*(t_c_tail)^4) * (1.34*M^0.18*(cos(Lambda_m))^0.28); %horizontal tail form factor
S_wet_ht = tail_area_h * 2 * 1.02; %wetted area of horizontal tail

FF_vt = (1 + 0.6/(x_c_max)*(t_c_tail) + 100*(t_c_tail)^4) * (1.34*M^0.18*(cos(Lambda_m_vtail))^0.28); %vertical tail form factor
S_wet_vt = S_vt * 2 * 1.02; %wetted area of vertical tail

Re_t = rho * V_CR * c_tail_h; %Reynolds number with chord as ref length
C_f_t = 0.455 / (log10(Re_t))^2.58;


Q_t = 1.05;

%% NACELLES

%we don't have nacelles yet

%% TOTAL PARASITE DRAG
Re = rho * V_CR * c;
C_f = 0.455 / (log10(Re))^2.58;

C_D_0_f = (FF_f * Q_f * C_f_f * S_wet_f) / S_ref_f;
C_D_0_w = (FF_w * Q_w * C_f_w * S_wet_w) / wing_area_total;
C_D_0_ht = (FF_t * Q_t * C_f_t * S_wet_ht) / (tail_area_h);
C_D_0_vt = (FF_t * Q_t * C_f_t * S_wet_vt) / (tail_area_h);


%TOTAL PARASITIC DRAG
C_D_0 = C_D_0_f + C_D_0_w + C_D_0_ht + C_D_0_vt; 

fprintf("\nFuselage CD0: %.4f", C_D_0_f)
fprintf("\nWing CD0: %.4f", C_D_0_w)
fprintf("\nTail CD0: %.4f", C_D_0_t)
fprintf("\nTotal CD0: %.4f", C_D_0)

