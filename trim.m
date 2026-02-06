CL_trim = linspace(-1,2,1000);

E = 0.3;
CL0_t = 0;
CL0_w = 0.36;

x_ac_w = x_ac;
x_cg_chosen = (x_cg_aft + x_cg_forward) * 0.5;
VH = selected_Sh_S * lt_c;

CL_det = CLa_t/pi() * (acos(1-2*E) + 2*sqrt(E*(1-E)));
CM0 = Cm0_w;
CL_de = selected_Sh_S * CL_det;

CM_de = CL_de + selected_Sh_S * (x_cg_chosen - x_ac_w) - CL_det * VH;
CL_0 = CL0_w + selected_Sh_S * CL0_t;

CL_a = CLa_w + selected_Sh_S * CLa_t * (1 - de_da);
CM_a = CL_a * (x_cg_chosen - x_ac) - CLa_t * (1-de_da) * VH;

a_trim = rad2deg((CM0 * CL_de + CM_de * (CL_trim - CL_0)) / (CL_a * CM_de - CL_de * CM_a));
de_trim = rad2deg(-((CM0 * CL_a + CM_a*(CL_trim - CL_0)) / (CL_a * CM_de - CL_de * CM_a)));

figure()
plot(CL_trim,a_trim)
xlabel("CL_{trim}")
ylabel("\alpha_{trim} [deg]")
title("Trim Angle of Attack vs. Trim Lift Coefficient")
grid on

figure()
plot(CL_trim, de_trim)
xlabel("CL_{trim}")
ylabel("\delta_{trim} [deg]")
title("Trim Elevator Deflection vs. Trim Lift Coefficient")
grid on

