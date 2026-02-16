
fprintf("\nCONTROL SURFACES AND SERVO TORQUE ------------------------\n")

v = 26.2 * 2.237; % max velocity found frum thrust and drag plot

elevator_area = SE_Sh * tail_area_h;
aileron_area = SA_S * wing_area_total;
rudder_area = SR_Sv * tail_area_v; 


elevator_c = CE_Ch * c_tail_h;
elevator_b = elevator_area / elevator_c;

aileron_c = CA_C * c;
aileron_b = 0.5 * aileron_area / aileron_c;

rudder_b = b_tail_v;
rudder_c = rudder_area / rudder_b;

elevator_demax = 25; % deg
aileron_demax = 30; % deg
rudder_demax = 30; % deg
servo_demax = 30;

elevator_torque = (8.5e-6) * ((elevator_c * 100)^2 * (elevator_b * 100)* v^2 * sind(elevator_demax) * tand(elevator_demax) / tand(servo_demax));
rudder_torque = (8.5e-6) * ((rudder_c* 100)^2 * (rudder_b * 100)* v^2 * sind(rudder_demax) * tand(rudder_demax) / tand(servo_demax));
aileron_torque = (8.5e-6) * ((aileron_c * 100)^2 * (aileron_b * 100) * v^2 * sind(aileron_demax) * tand(aileron_demax) / tand(servo_demax));

fprintf("Rudder Chord (percent of tail root chord) = %0.4f\n",rudder_c / c_tail_v_root)
fprintf("Control Surface Area:\n    Aileron: %.4f m^2\n    Elevator: %.4f m^2\n    Rudder: %.4f m^2\n",aileron_area,elevator_area,rudder_area)
fprintf("Torque Required:\n    elevator: %.4f oz-in or %.4f kg-cm\n", elevator_torque, elevator_torque / 13.887)
fprintf("    rudder: %.4f oz-in or %.4f kg-cm\n", rudder_torque, rudder_torque / 13.887)
fprintf("    aileron: %.4f oz-in or %.4f kg-cm\n", aileron_torque, aileron_torque / 13.887)