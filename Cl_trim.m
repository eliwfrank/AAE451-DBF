St_Sw = 0.089;
x_cg_def = 0.22;
Vh = St_Sw*lt_c;
CL_0 = 0.36;
CL_0t = 0;
Cma_w = (3.2*10^-3)*(180/pi); %Cm_alpha wing in rad
Cma_t = (-2*10^-3)*(180/pi); %Cm_alpha horz tail in rad

CL_trim = linspace(-1,2,1000);

E = (0.3*c_tail_h)/c_tail_h; % elevator to tail chord ratio
CLt_deltae = (CLa_t/pi)*(acos(1-2*E) + 2*sqrt(E*(1-E))); 
CLd_et = CLt_deltae;

Cmd_e = (CLt_deltae*St_Sw*(x_cg_def - x_ac)) - (CLd_et*Vh);

a_trim = ((Cm0_w*CLt_deltae) + Cmd_e*(CL_trim-CL_0))/((CLa_w*Cmd_e) - (CLt_deltae*Cma_w));
del_etrim = -((Cm0_w*CLa_w) + Cma_w*(CL_trim - CL_0))/((CLa_w*Cmd_e) - (CLt_deltae*Cma_w));

a_trim = rad2deg(a_trim);
del_etrim = rad2deg(del_etrim);

CL_target = 0.69;

alpha_trim_MTOW = interp1(CL_trim, a_trim, CL_target, 'linear');      % deg
del_trim_MTOW = interp1(CL_trim, del_etrim, CL_target, 'linear');   % deg

fprintf('At CL = %.2f: alpha_trim = %.3f deg, delta_e_trim = %.3f deg\n', ...
        CL_target, alpha_trim_MTOW, del_trim_MTOW);

CL_w1 = CL_0 + CLa_w*deg2rad(a_trim);
CL_t1 = CL_0t + CLa_t*deg2rad(a_trim) + CLt_deltae*deg2rad(del_etrim);

CL_t2 = CL_0t + CLa_t*deg2rad(a_trim); % CL Clean

K1 = 0.053;
Cd01 = 0.03238;

Cd1 = Cd01 + K1*CL_w1.^2 + St_Sw*((CL_t1.^2)/(pi*4*0.75));
Cd2 = Cd01 + K1*CL_w1.^2 + St_Sw*((CL_t2.^2)/(pi*4*0.75));


figure()
plot(CL_trim,a_trim,'r-',LineWidth=1.5)
grid on
xlabel('$C_{L_{trim}}$','Interpreter','latex')
ylabel('$\alpha_{trim} (^{\circ})$','Interpreter','latex')
title('$\alpha_{trim} vs. C_{L_{trim}}$','Interpreter','latex')

figure()
plot(CL_trim,del_etrim,'r-',LineWidth=1.5)
grid on
xlabel('$C_{L_{trim}}$','Interpreter','latex')
ylabel('$\delta_{e_{trim}} (^{\circ})$','Interpreter','latex')
title('$\delta_{e_{trim}} vs. C_{L_{trim}}$','Interpreter','latex')

figure()
plot(CL_trim,Cd1,'r-',LineWidth=1.5)
grid on
hold on
plot(CL_trim,Cd2,'b-',LineWidth=1.5)
xlabel('$C_{L_{trim}}$','Interpreter','latex')
ylabel('$C_D$','Interpreter','latex')
legend('Trimmed','Clean','Location','best')
title('Aircraft Drag Polar')

