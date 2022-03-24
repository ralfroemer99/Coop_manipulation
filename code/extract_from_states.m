%%This script extracts the necessary data from the saved variables in the
%%respective dataset. Simply change the dataset ('.../xy.mat') to something
%%that is in the respective folder. Change the wanted dataset in the load
%%and clear command at the top and the bottom.

ddx_left = States(2:7,:);
dx_left = States(8:13,:);
x_left = States(14:20,:);

ddx_right = States(21:26,:);
dx_right = States(27:32,:);
x_right = States(33:39,:);

ddx_jd = States(40:51,:);
dx_jd = States(52:63,:);
x_jd = States(64:77,:);

dx_d = States(78:89,:);
flag_start_coop_des = States(90,:);
flag_start_coop_dyn = States(91,:);
flag_dyn_over = States(92,:);

flag_coop = length(dx_jd)-length(flag_start_coop_dyn(flag_start_coop_dyn>0));
disp('Starting of dynamic cooperative manipulation at')
disp(['Sample: ',num2str(flag_coop)])

flag_coop_done = length(dx_jd)-length(flag_dyn_over(flag_dyn_over>0));
disp('Ending dynamic cooperative manipulation at')
disp(['Sample: ',num2str(flag_coop_done)])

if size(States,1) > 93
    forces_left = States(93:98,:);
    forces_right = States(99:104,:);
end

clear States_25