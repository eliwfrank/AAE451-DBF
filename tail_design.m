x_cg = linspace(0,1,10000);
x_np = SM + x_cg;
lt_c = lt/c;

%% forward limit (stability)
forward_lim_STS = (x_cg - x_ac + SM) ./ (etaH * (1 - de_da) * lt_c - (x_cg - x_ac + SM));

%% stall recovery control
CLt_nose_down = CLa_t * aoa_tail_stall; % nose down
Cm_req_rec = -CL_max * SM + Cm0_w;
stall_recovery_STS = (Cm0_w + CL_max * (x_cg - x_ac) - Cm_req_rec) ./ (CLt_nose_down * etaH * (lt_c - x_cg + x_ac));

%% TO Rotation, nose up control
takeoff_rot_STS = (Cm0_w + CL_R * (x_cg - x_ac) - Cm_req_rot) ./ (CLt_nose_up * etaH * (lt_c - x_cg + x_ac));

%% Aft Limit (stability, SM = 0)
aft_lim_STS = (x_cg - x_ac - 0) ./ (etaH .* (1 - de_da) .* lt_c - (x_cg - x_ac + 0));

x_cg_forward = interp1(takeoff_rot_STS,x_cg, Sh_S);
x_cg_aft = interp1(forward_lim_STS, x_cg,Sh_S);
tail_area_h = Sh_S * wing_area_total;
b_tail_h = sqrt(tail_area_h * AR_tail);
c_tail_h = tail_area_h / b_tail_h;

%% Vertical Tail Area
VC_vt = 0.04; % chosen from values in lecture
lt = 1; % same as horizontal tail arm
lambda_vt = 0.4;
S_vt = (VC_vt*b*wing_area_total)/(lt);
AR_vt = 2; % chosen from values 
b_tail_v = sqrt(S_vt*AR_vt);
c_tail_v_root = ((2*b_tail_v) / (AR_vt*(1 + lambda_vt)));
c_tail_v_tip = lambda_vt * c_tail_v_root;

%% Display Results
figure()
plot(x_cg,forward_lim_STS,"k")
hold on
grid on
plot(x_cg,stall_recovery_STS,"k--",LineWidth=2)
plot(x_cg,takeoff_rot_STS,"k:",LineWidth=2)
plot(x_cg, aft_lim_STS,"g")
yline(Sh_S)
plot(x_cg_forward,Sh_S,"r.",markersize=20)
plot(x_cg_aft,Sh_S,"r.",markersize=20)
plot([x_cg_forward x_cg_aft], [Sh_S Sh_S], 'r-', 'LineWidth',2)
plot(0.2694,0.0409)
legend("Forward limit (stability)","Aft Limit (Stall recovery control)","Forward Limit (Nose-up Control)","Aft Limit (Stability, SM = 0)","selected S_h/S","Forward Limit","Aft Limit","Feasible cg range", "Intersection Point")
title("Scissor Plot")
ylabel("S_h/S (Horizontal Tail Area Ratio)")
xlabel("x_{cg}/c (Center of Gravity Position)")

fprintf("\nForward Limit: x_cg/c = %.4f, Sh/S = %.3f\n",x_cg_forward, Sh_S)
fprintf("Aft Limit: x_cg/c = %.4f, Sh/S = %.3f\n\n",x_cg_aft, Sh_S)
fprintf("Horizontal Tail Area: %f [m^2]\n",tail_area_h)
fprintf("Tail Chord length: %.4f [m]\nTail Span: %.4f [m]\nTail Half Span: %0.4f [m]\n\n",c_tail_h,b_tail_h,b_tail_h/2)
fprintf("Vertical Stabilizer Area:  %.4f [m^2]\n", S_vt);
fprintf("Vertical Stabilizer Span:  %.4f [m^2]\n", b_tail_v);
fprintf("Vertical Stabilizer Root Chord:  %.4f [m^2]\n", c_tail_v_root);
fprintf("Vertical Stabilizer Tip Chord:  %.4f [m^2]\n", c_tail_v_tip);

