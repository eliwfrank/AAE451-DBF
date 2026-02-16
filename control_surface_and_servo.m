v = 26.2 * 2.237; % max velocity found frum thrust and drag plot

elevator_area = 0.4 * tail_area_h; % E = 0.3
aileron_area = 0.12 * wing_area_total;
rudder_area = 0.45 * tail_area_v; % 0.45 chosen for Sr/Sv from saedry textbook

elevator_c = 0.4 * c_tail_h;
elevator_b = elevator_area / elevator_c;

aileron_c = 0.3 * c;
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

fprintf("Torque Required:\nelevator: %.4f oz-in or %.4f kg-cm\n", elevator_torque, elevator_torque / 13.887)
fprintf("rudder: %.4f oz-in or %.4f kg-cm\n", rudder_torque, rudder_torque / 13.887)
fprintf("aileron: %.4f oz-in or %.4f kg-cm\n", aileron_torque, aileron_torque / 13.887)