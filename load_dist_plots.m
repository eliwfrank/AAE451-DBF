%% Initialization
y = linspace(0,b/2,1000);
x = linspace(0,L,1000).';

%% Wing
figure();
subplot(3,2,1)
plot(y, L_dist, 'r-')
xlabel('x (m)')
ylabel('Lift Distribution (N/m)')
title('Schrenk Lift Approximation')
grid on

subplot(3,2,2)
plot(y, V, 'b-')
xlabel('x (m)')
ylabel('V (N)')
title('Shear Force Distribution')
grid on

subplot(3,2,3)
plot(y, My, 'k-')
xlabel('x (m)')
ylabel('M_{bending} (N*m)')
title('Bending Moment Distribution')
grid on

subplot(3,2,4)
plot(y, T, 'g-')
xlabel('x (m)')
ylabel('M_{torsion} (N*m)')
title('Torsion Moment Distribution from c/4 to x_{sc}')
grid on

subplot(3,2,5)
plot(y,w_carbon); 
grid on
hold on
plot(y,w_aluminum); 
plot(y,w_balsa)
xlabel('x (m)'); 
ylabel('w(x) (m)');
legend('Carbon Fiber Hollow Tube','Aluminum Hollow Tube','Balsa Solid Tube','Location','best')
title('Vertical deflection (Euler–Bernoulli)');

% Bending Stress Calculation (Not dependent on material, is only dependent on Cross Sectional Area --> Diameter of Rods) 
subplot(3,2,6)
plot(y,stress_wc); 
grid on
hold on
plot(y,stress_wb); 
xlabel('x (m)'); 
ylabel('\sigma_{zz} (Pa)');
legend('Carbon Fiber Hollow Tube','Balsa Solid Tube','Location','best')
title('Bending Stress Distribution');

sgtitle('Wing Load Distribution')
%% Horizontal Stabilizer

figure();
subplot(3,2,1)
plot(y_t, L_dist_t, 'r-')
xlabel('x (m)')
ylabel('Lift Distribution (N/m)')
title('Schrenk Lift Approximation')
grid on

subplot(3,2,2)
plot(y_t, V_t, 'b-')
xlabel('x (m)')
ylabel('V (N)')
title('Shear Force Distribution')
grid on

subplot(3,2,3)
plot(y_t, My_t, 'k-')
xlabel('x (m)')
ylabel('M_{bending} (N*m)')
title('Bending Moment Distribution')
grid on

subplot(3,2,4)
plot(y_t, T_h, 'g-')
xlabel('x (m)')
ylabel('M_{torsion} (N*m)')
title('Torsion Moment Distribution from c/4 to x_{sc}')
grid on

subplot(3,2,5)
plot(y_t,w_carbon_t); 
grid on
hold on
plot(y_t,w_aluminum_t); 
plot(y_t,w_balsa_t)
xlabel('x (m)'); 
ylabel('w(x) (m)');
legend('Carbon Fiber Hollow Tube','Aluminum Hollow Tube','Balsa Solid Tube','Location','best')
title('Vertical deflection (Euler–Bernoulli)');

% Bending Stress Calculation (Not dependent on material, is only dependent on Cross Sectional Area --> Diameter of Rods)
subplot(3,2,6)
plot(y_t,stress_tc); 
grid on
hold on
plot(y_t,stress_tb); 
xlabel('x (m)'); 
ylabel('\sigma_{zz} (Pa)');
legend('Carbon Fiber Hollow Tube','Balsa Solid Tube','Location','best')
title('Bending Stress Distribution');

sgtitle('Horizontal Stabilizer Load Distribution')
%% Tail Boom - Symmetric
%be consistent
figure();
subplot(2,2,1)
plot(x,w_boom);
grid on;
ylabel('w (m)');
xlabel('x (m)');
title('Vertical Deflection');

subplot(2,2,2)
plot(x, V_boom);
grid on;
xlabel('x (m)');
ylabel('V (N)');
title('Vertical Shear');

subplot(2,2,3)
plot(x, M_boom);
grid on;
xlabel('x (m)');
ylabel('M_{bending} (N*m)');
title('Vertical Bending Moment');

% Bending Stress Calculation (Not dependent on material, is only dependent on Cross Sectional Area --> Diameter of Rods)
subplot(2,2,4)
plot(x,stress_boom); 
grid on
xlabel('x (m)'); 
ylabel('\sigma_{zz} (Pa)');
legend('Carbon Fiber Hollow Tube','Location','best')
title('Vertical Bending Stress Distribution');

sgtitle('Boom Load Distribution - Symmetric (T = 0)')

%% Tail Boom - Asymmetric
figure();
subplot(3,2,1)
plot(v_boom,x);
grid on;
xlabel('v (m)');
ylabel('x (m)');
title('Horizontal Deflection');

subplot(3,2,2)
plot(V_yboom, x);
grid on;
xlabel('V (N)');
ylabel('x (m)');
title('Horizontal Shear');

subplot(3,2,3)
plot(M_yboom, x);
grid on;
xlabel('M_{bending} (N*m)');
ylabel('x (m)');
title('Horizontal Bending Moment');

subplot(3,2,4)
plot(T_boom, x);
grid on;
xlabel('M_{torsion} (N*m)');
ylabel('x (m)');
title('Torsion Moment at Force Centroid'); % Ask Thiago if correct

% Bending Stress Calculation (Not dependent on material, is only dependent on Cross Sectional Area --> Diameter of Rods)
subplot(3,2,[5 6])
plot(stress_booma,x); 
grid on
xlabel('\sigma_{zz} (Pa)'); 
ylabel('x (m)');
legend('Carbon Fiber Hollow Tube','Location','best')
title('Horizontal Bending Stress Distribution');

sgtitle('Boom Load Distribution - Asymmetric')