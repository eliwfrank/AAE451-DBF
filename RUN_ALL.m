clear
clc
close all

run("param.m") 
run("Design_Point.m")
run("weight.m")
run("wing_design.m")
run("tail_design.m")
run("trim_and_drag_buildup.m")

run("prop_system.m")
run("battery_sizing.m")
run("control_surface_and_servo.m")
run("vn_diagram.m")

run("parameter_update_check.m")
run("wing_design.m")

run("component_load_dist.m")
run("wing_box_flexural.m")
run("wing_box_symb_D_box.m")
% close all