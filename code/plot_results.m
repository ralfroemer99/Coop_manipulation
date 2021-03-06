%% Plot results
plot_errors = 0;
plot_uncertainty = 1;
plot_model_acc = 0;
plot_r = 1;
plot_m = 1;
plot_j = 1;

crop_size = 10;

if plot_r
%% r1 estimated by manipulator 1
num_sig_reg = 3;
estim = r1_dyn_man1_save;
estim_true = O_r_1;
unc_est = sqrt(Sigma_r1_dyn_man1_save);

diff_x = estim(:,1);
upper_x = diff_x + num_sig_reg*unc_est(:,1);
lower_x = diff_x - num_sig_reg*unc_est(:,1);

diff_y = estim(:,2);
upper_y = diff_y + num_sig_reg*unc_est(:,2);
lower_y = diff_y - num_sig_reg*unc_est(:,2);

diff_z = estim(:,3);
upper_z = diff_z + num_sig_reg*unc_est(:,3);
lower_z = diff_z - num_sig_reg*unc_est(:,3);

% Plot
figure()
plot(t,diff_x,'LineWidth',.7)
grid on
hold on
plot(t,diff_y,'LineWidth',.7)
plot(t,diff_z,'LineWidth',.7)
plot([0,t(end)],[estim_true(1),estim_true(1)], 'b--')
plot([0,t(end)],[estim_true(2),estim_true(2)], 'r--')
plot([0,t(end)],[estim_true(3),estim_true(3)], 'k--')
if plot_uncertainty
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_x(1:crop_size:length(estim))',fliplr(lower_x(1:crop_size:length(estim))')],'b','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_y(1:crop_size:length(estim))',fliplr(lower_y(1:crop_size:length(estim))')],'r','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_z(1:crop_size:length(estim))',fliplr(lower_z(1:crop_size:length(estim))')],'y','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
end
legend('r1_x','r1_y','r1_z','r1_x true','r1_y true','r1_z true','unc_x','unc_y','unc_z')   
xlabel('Time (s)')
ylabel('Estimation r1 (m)')
title('Estimation of r_1 from Model 1')

% Plot errors
if plot_errors
low_border = 1.5e-1;
estim_true = O_r_1;
unc = Sigma_r1_dyn_man1_save;
if length(estim_true) > 1
    diff = VecNorm((estim-estim_true')',2);
    unc_est = sum(sqrt(unc),2);
else
    diff = abs((estim-estim_true')');
    unc_est = sqrt(unc);
end
upper = diff + num_sig_reg*unc_est';
lower = low_border*ones(1,length(estim));
figure()
semilogy(t,diff)
grid on
hold on
% fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper(1:crop_size:length(estim)),fliplr(lower(1:crop_size:length(estim)))],'k','FaceAlpha',.15,'EdgeAlpha',.4,'LineStyle','--','EdgeColor','k')
% fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper(1:crop_size:length(estim)),fliplr(diff(1:crop_size:length(estim)))],'k','FaceAlpha',.15,'EdgeAlpha',.4,'LineStyle','--','EdgeColor','k')
xlabel('Time (s)')
ylabel('Estimation r1 (m)')
legend('diff','unc')
title('Error of r_1 from Model 1')
end

%% r2 estimated by manipulator 2
estim = r2_dyn_man2_save;
estim_true = O_r_2;
unc_est = sqrt(Sigma_r2_dyn_man1_save);

diff_x = estim(:,1);
upper_x = diff_x + num_sig_reg*unc_est(:,1);
lower_x = diff_x - num_sig_reg*unc_est(:,1);

diff_y = estim(:,2);
upper_y = diff_y + num_sig_reg*unc_est(:,2);
lower_y = diff_y - num_sig_reg*unc_est(:,2);

diff_z = estim(:,3);
upper_z = diff_z + num_sig_reg*unc_est(:,3);
lower_z = diff_z - num_sig_reg*unc_est(:,3);

% Plot
figure()
plot(t,diff_x,'LineWidth',.7)
grid on
hold on
plot(t,diff_y,'LineWidth',.7)
plot(t,diff_z,'LineWidth',.7)
plot([0,t(end)],[estim_true(1),estim_true(1)], 'b--')
plot([0,t(end)],[estim_true(2),estim_true(2)], 'r--')
plot([0,t(end)],[estim_true(3),estim_true(3)], 'k--')
if plot_uncertainty
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_x(1:crop_size:length(estim))',fliplr(lower_x(1:crop_size:length(estim))')],'b','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_y(1:crop_size:length(estim))',fliplr(lower_y(1:crop_size:length(estim))')],'r','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_z(1:crop_size:length(estim))',fliplr(lower_z(1:crop_size:length(estim))')],'y','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
end
legend('r2_x','r2_y','r2_z','r2_x true','r2_y true','r2_z true','unc_x','unc_y','unc_z')   
xlabel('time (s)')
ylabel('Estimation r2 (m)')
title('Estimation of r_2 from Model 2')

% Plot errors
if plot_errors
low_border = 1.5e-1;
unc = Sigma_r2_dyn_man1_save;
if length(estim_true) > 1
    diff = VecNorm((estim-estim_true')',2);
    unc_est = sum(sqrt(unc),2);
else
    diff = abs((estim-estim_true')');
    unc_est = sqrt(unc);
end
upper = diff + num_sig_reg*unc_est';
lower = low_border*ones(1,length(estim));
figure()
semilogy(t,diff)
grid on
hold on
% fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper(1:crop_size:length(estim)),fliplr(lower(1:crop_size:length(estim)))],'k','FaceAlpha',.15,'EdgeAlpha',.4,'LineStyle','--','EdgeColor','k')
% fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper(1:crop_size:length(estim)),fliplr(diff(1:crop_size:length(estim)))],'k','FaceAlpha',.15,'EdgeAlpha',.4,'LineStyle','--','EdgeColor','k')
xlabel('Time (s)')
ylabel('Estimation r2 (m)')
legend('diff','unc')
title('Error of r_2 from Model 2')
end


%% r1 consensus estimate
estim = r_cons_save(:,1:3);
estim_true = O_r_1;
unc_est = sqrt(Sigma_r_cons_save);

diff_x = estim(:,1);
upper_x = diff_x + num_sig_reg*unc_est(:,1);
lower_x = diff_x - num_sig_reg*unc_est(:,1);

diff_y = estim(:,2);
upper_y = diff_y + num_sig_reg*unc_est(:,2);
lower_y = diff_y - num_sig_reg*unc_est(:,2);

diff_z = estim(:,3);
upper_z = diff_z + num_sig_reg*unc_est(:,3);
lower_z = diff_z - num_sig_reg*unc_est(:,3);

% Plot
figure()
plot(t,diff_x,'LineWidth',.7)
grid on
hold on
plot(t,diff_y,'LineWidth',.7)
plot(t,diff_z,'LineWidth',.7)
plot([0,t(end)],[estim_true(1),estim_true(1)], 'b--')
plot([0,t(end)],[estim_true(2),estim_true(2)], 'r--')
plot([0,t(end)],[estim_true(3),estim_true(3)], 'k--')
if plot_uncertainty
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_x(1:crop_size:length(estim))',fliplr(lower_x(1:crop_size:length(estim))')],'b','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_y(1:crop_size:length(estim))',fliplr(lower_y(1:crop_size:length(estim))')],'r','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_z(1:crop_size:length(estim))',fliplr(lower_z(1:crop_size:length(estim))')],'y','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
end
legend('r1_x','r1_y','r1_z','r1_x true','r1_y true','r1_z true','unc_x','unc_y','unc_z')   
xlabel('Time (s)')
ylabel('Estimation r1 (m)')
title('Consensus estimate of r_1')

% Plot errors
if plot_errors
low_border = 1.5e-1;
unc = Sigma_r_cons_save;
if length(estim_true) > 1
    diff = VecNorm((estim-estim_true')',2);
    unc_est = sum(sqrt(unc),2);
else
    diff = abs((estim-estim_true')');
    unc_est = sqrt(unc);
end
upper = diff + num_sig_reg*unc_est';
lower = low_border*ones(1,length(estim));
figure()
semilogy(t,diff)
grid on
hold on
% fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper(1:crop_size:length(estim)),fliplr(lower(1:crop_size:length(estim)))],'k','FaceAlpha',.15,'EdgeAlpha',.4,'LineStyle','--','EdgeColor','k')
% fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper(1:crop_size:length(estim)),fliplr(diff(1:crop_size:length(estim)))],'k','FaceAlpha',.15,'EdgeAlpha',.4,'LineStyle','--','EdgeColor','k')
xlabel('Time (s)')
ylabel('YLABEL')
legend('diff','unc')
title('Error of consensus estimate of r_1')
end


%% r2 consensus estimate
estim = r_cons_save(:,4:6);
estim_true = O_r_2;
unc_est = sqrt(Sigma_r_cons_save);

diff_x = estim(:,1);
upper_x = diff_x + num_sig_reg*unc_est(:,1);
lower_x = diff_x - num_sig_reg*unc_est(:,1);

diff_y = estim(:,2);
upper_y = diff_y + num_sig_reg*unc_est(:,2);
lower_y = diff_y - num_sig_reg*unc_est(:,2);

diff_z = estim(:,3);
upper_z = diff_z + num_sig_reg*unc_est(:,3);
lower_z = diff_z - num_sig_reg*unc_est(:,3);

% Plot
figure()
plot(t,diff_x,'LineWidth',.7)
grid on
hold on
plot(t,diff_y,'LineWidth',.7)
plot(t,diff_z,'LineWidth',.7)
plot([0,t(end)],[estim_true(1),estim_true(1)], 'b--')
plot([0,t(end)],[estim_true(2),estim_true(2)], 'r--')
plot([0,t(end)],[estim_true(3),estim_true(3)], 'k--')
if plot_uncertainty
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_x(1:crop_size:length(estim))',fliplr(lower_x(1:crop_size:length(estim))')],'b','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_y(1:crop_size:length(estim))',fliplr(lower_y(1:crop_size:length(estim))')],'r','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_z(1:crop_size:length(estim))',fliplr(lower_z(1:crop_size:length(estim))')],'y','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
end
legend('r2_x','r2_y','r2_z','r2_x true','r2_y true','r2_z true','unc_x','unc_y','unc_z')   
xlabel('Time (s)')
ylabel('Estimation r2 (m)')
title('Consensus estimate of r_2')

% Plot errors
if plot_errors
low_border = 1.5e-1;
unc = Sigma_r_cons_save;
if length(estim_true) > 1
    diff = VecNorm((estim-estim_true')',2);
    unc_est = sum(sqrt(unc),2);
else
    diff = abs((estim-estim_true')');
    unc_est = sqrt(unc);
end
% upper = diff + num_sig_reg*unc_est';
% lower = low_border*ones(1,length(estim));
figure()
semilogy(t,diff)
grid on
hold on
xlabel('Time (s)')
ylabel('YLABEL')
legend('diff','unc')
title('Error of consensus estimate of r_2')
end
end

if plot_m
%% mo consensus estimate
estim = mo_cons_save;
estim_true = m_o;

unc_est = sqrt(Sigma_mo_cons_save);

upper_x = estim + num_sig_reg*unc_est;
lower_x = estim - num_sig_reg*unc_est;

% Plot
figure()
plot(t,estim,'LineWidth',.7)
grid on
hold on
plot([0,t(end)],[estim_true,estim_true], 'b--')
fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_x(1:crop_size:length(estim))',fliplr(lower_x(1:crop_size:length(estim))')],'b','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
legend('m_o','m_o true','unc_m')
xlabel('Time (s)')
ylabel('Estimation m_o (kg)')
title('Consensus estimate for m_o')

% Plot errors
if plot_errors
low_border = 1e-1;
unc = Sigma_mo_cons_save;
if length(estim_true) > 1
    diff = VecNorm((estim-estim_true')');
    unc_est = sum(sqrt(unc),2);
else
    diff = abs((estim-estim_true')');
    unc_est = sqrt(unc);
end
% upper = diff + num_sig_reg*unc_est';
% lower = low_border*ones(1,length(estim));
figure()
semilogy(t,diff)
grid on
hold on
xlabel('time (s)')
ylabel('YLABEL')
legend('diff','unc')
title('Consensus estimation for mo')
end
end

if plot_j
%% jo estimated by manipulator 1
estim = jo_dyn_man1_save;
estim_true = j_guess;

% Plot
figure()
plot(t,estim(:,1),'b','LineWidth',.7)
grid on
hold on
plot(t,estim(:,2),'r','LineWidth',.7)
plot(t,estim(:,3),'k','LineWidth',.7)
plot(t,estim(:,4),'m','LineWidth',.7)
plot(t,estim(:,5),'c','LineWidth',.7)
plot(t,estim(:,6),'g','LineWidth',.7)
plot([0,t(end)],[estim_true(1), estim_true(1)], 'b--')
plot([0,t(end)],[estim_true(2), estim_true(2)], 'r--')
plot([0,t(end)],[estim_true(3), estim_true(3)], 'k--')
plot([0,t(end)],[estim_true(4), estim_true(4)], 'm--')
plot([0,t(end)],[estim_true(5), estim_true(5)], 'c--')
plot([0,t(end)],[estim_true(6), estim_true(6)], 'g--')

legend('jo xx','jo xy','jo xz','jo yy','jo yz','jo zz','jo xx est.','jo xy est.','jo xz est.','jo yy est.','jo yz est.','jo zz est.')
xlabel('time (s)')
ylabel('Estimation jo (kg m^2)')
title('Estimation for jo from Model 1')


%% jo estimated by manipulator 2
estim = jo_dyn_man2_save;
estim_true = j_guess;

% Plot
figure()
plot(t,estim(:,1),'b','LineWidth',.7)
grid on
hold on
plot(t,estim(:,2),'r','LineWidth',.7)
plot(t,estim(:,3),'k','LineWidth',.7)
plot(t,estim(:,4),'m','LineWidth',.7)
plot(t,estim(:,5),'c','LineWidth',.7)
plot(t,estim(:,6),'g','LineWidth',.7)
plot([0,t(end)],[estim_true(1), estim_true(1)], 'b--')
plot([0,t(end)],[estim_true(2), estim_true(2)], 'r--')
plot([0,t(end)],[estim_true(3), estim_true(3)], 'k--')
plot([0,t(end)],[estim_true(4), estim_true(4)], 'm--')
plot([0,t(end)],[estim_true(5), estim_true(5)], 'c--')
plot([0,t(end)],[estim_true(6), estim_true(6)], 'g--')

legend('jo xx','jo xy','jo xz','jo yy','jo yz','jo zz','jo xx est.','jo xy est.','jo xz est.','jo yy est.','jo yz est.','jo zz est.')
xlabel('time (s)')
ylabel('Estimation jo (kg m^2)')
title('Estimation for jo from Model 2')


%% jo consensus estimate
estim = jo_cons_save;
estim_true = j_guess;
unc_est = sqrt(Sigma_jo_cons_save);

diff_x = estim(:,1);
upper_x = diff_x + num_sig_reg*unc_est(:,1);
lower_x = diff_x - num_sig_reg*unc_est(:,1);

diff_y = estim(:,4);
upper_y = diff_y + num_sig_reg*unc_est(:,4);
lower_y = diff_y - num_sig_reg*unc_est(:,4);

diff_z = estim(:,6);
upper_z = diff_z + num_sig_reg*unc_est(:,6);
lower_z = diff_z - num_sig_reg*unc_est(:,6);

% Plot
figure()
plot(t,estim(:,1),'b','LineWidth',.7)
grid on
hold on
plot(t,estim(:,2),'r','LineWidth',.7)
plot(t,estim(:,3),'k','LineWidth',.7)
plot(t,estim(:,4),'m','LineWidth',.7)
plot(t,estim(:,5),'c','LineWidth',.7)
plot(t,estim(:,6),'g','LineWidth',.7)
plot([0,t(end)],[estim_true(1), estim_true(1)], 'b--')
plot([0,t(end)],[estim_true(2), estim_true(2)], 'r--')
plot([0,t(end)],[estim_true(3), estim_true(3)], 'k--')
plot([0,t(end)],[estim_true(4), estim_true(4)], 'm--')
plot([0,t(end)],[estim_true(5), estim_true(5)], 'c--')
plot([0,t(end)],[estim_true(6), estim_true(6)], 'g--')
if plot_uncertainty
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_x(1:crop_size:length(estim))',fliplr(lower_x(1:crop_size:length(estim))')],'b','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_y(1:crop_size:length(estim))',fliplr(lower_y(1:crop_size:length(estim))')],'r','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
    fill([t(1:crop_size:end),fliplr(t(1:crop_size:end))],[upper_z(1:crop_size:length(estim))',fliplr(lower_z(1:crop_size:length(estim))')],'y','FaceAlpha',.1,'EdgeAlpha',.25,'LineStyle','--','EdgeColor','k')
end
legend('jo xx','jo xy','jo xz','jo yy','jo yz','jo zz','jo xx est.','jo xy est.','jo xz est.','jo yy est.','jo yz est.','jo zz est.')
xlabel('time (s)')
ylabel('Estimation jo (kg m^2)')
title('Consensus estimation for jo')
end

%% Model accuracies
if plot_model_acc
% Translational model
figure();
plot(t, VecNorm(t_man1_save'), 'r'); hold on
plot(t, VecNorm(PhiTheta_man1_save'), 'r--'); hold on
plot(t, VecNorm(t_man2_save'), 'b'); hold on
plot(t, VecNorm(PhiTheta_man2_save'), 'b--'); hold on
legend({'$t_1$', '$\Phi_1*\theta_1$', '$t_2$', '$\Phi_2*\theta_2$'}, 'Interpreter', 'latex', 'FontSize', 12);
title('Accuracy of translational model');

% Detailed examination for manipulator 1
figure();
plot(t, t_man1_save(:,1), 'r'); hold on
plot(t, t_man1_save(:,2), 'g'); hold on
plot(t, t_man1_save(:,3), 'b'); hold on
plot(t, PhiTheta_man1_save(:,1), 'r--');
plot(t, PhiTheta_man1_save(:,2), 'g--');
plot(t, PhiTheta_man1_save(:,3), 'b--');
legend({'$t_{1,1}$', '$t_{1,2}$', '$t_{1,3}$', '$\phi_{1,1}^T\theta_1$', '$\phi_{1,2}^T\theta_1$', '$\phi_{1,3}^T\theta_1$'}, 'Interpreter', 'latex', 'FontSize', 12);
title('Components of trans. model for manipulator 1');

% Detailed examination for manipulator 2
figure();
plot(t, t_man2_save(:,1), 'r'); hold on
plot(t, t_man2_save(:,2), 'g'); hold on
plot(t, t_man2_save(:,3), 'b'); hold on
plot(t, PhiTheta_man2_save(:,1), 'r--');
plot(t, PhiTheta_man2_save(:,2), 'g--');
plot(t, PhiTheta_man2_save(:,3), 'b--');
legend({'$t_{2,1}$', '$t_{2,2}$', '$t_{2,3}$', '$\phi_{2,1}^T\theta_2$', '$\phi_{2,2}^T\theta_2$', '$\phi_{2,3}^T\theta_2$'}, 'Interpreter', 'latex', 'FontSize', 12);
title('Components of trans. model for manipulator 2');

% Rotational model
figure();
plot(t, VecNorm(to_man1_save'), 'r'); hold on
plot(t, VecNorm(PhiThetao_man1_save'), 'r--'); hold on
plot(t, VecNorm(to_man2_save'), 'b'); hold on
plot(t, VecNorm(PhiThetao_man2_save'), 'b--'); hold on
legend({'$t_1$', '$\Phi_1\theta_1$', '$t_2$', '$\Phi_2\theta_2$'}, 'Interpreter', 'latex', 'FontSize', 12);
title('Accuracy of rotational model');

% Detailed examination for manipulator 1
figure();
plot(t, to_man1_save(:,1), 'r'); hold on
plot(t, to_man1_save(:,2), 'g'); hold on
plot(t, to_man1_save(:,3), 'b'); hold on
plot(t, PhiThetao_man1_save(:,1), 'r--');
plot(t, PhiThetao_man1_save(:,2), 'g--');
plot(t, PhiThetao_man1_save(:,3), 'b--');
legend({'$t_{1,1}^r$', '$t_{1,2}^r$', '$t_{1,3}^r$', '$\phi_{1,1}^{rT}\theta_1^r$', '$\phi_{1,2}^{rT}\theta_1^r$', '$\phi_{1,3}^{rT}\theta_1^r$'}, 'Interpreter', 'latex', 'FontSize', 12);
title('Components of rot. model for manipulator 1');

% Detailed examination for manipulator 2
figure();
plot(t, to_man2_save(:,1), 'r'); hold on
plot(t, to_man2_save(:,2), 'g'); hold on
plot(t, to_man2_save(:,3), 'b'); hold on
plot(t, PhiThetao_man2_save(:,1), 'r--');
plot(t, PhiThetao_man2_save(:,2), 'g--');
plot(t, PhiThetao_man2_save(:,3), 'b--');
legend({'$t_{2,1}$'}, {'$t_{2,2}$', '$t_{2,3}$', '$\phi_{2,1}^T\theta_2$', '$\phi_{2,2}^T\theta_2$', '$\phi_{2,3}^T\theta_2$'}, 'Interpreter', 'latex', 'FontSize', 12);
title('Components of rot. model for manipulator 2');
end

% REMOVE %%%%%
figure();
plot(t, to2_man1_save(:,1), 'r'); hold on
plot(t, to_man1_save(:,1), 'r--'); hold on
plot(t, to2_man1_save(:,2), 'g'); hold on
plot(t, to2_man1_save(:,3), 'b'); hold on
plot(t, to_man1_save(:,2), 'g--'); hold on
plot(t, to_man1_save(:,3), 'b--'); hold on
legend('forces','Tracking errors');
title('Target calculation via tracking errors and forces for man 1');

figure();
plot(t, to2_man2_save(:,1), 'r'); hold on
plot(t, to_man2_save(:,1), 'r--'); hold on
plot(t, to2_man2_save(:,2), 'g'); hold on
plot(t, to2_man2_save(:,3), 'b'); hold on
plot(t, to_man2_save(:,2), 'g--'); hold on
plot(t, to_man2_save(:,3), 'b--'); hold on
legend('forces','Tracking errors');
title('Target calculation via tracking errors and forces for man 2');

figure();
plot(t, imp1_l_save(:,1), 'r'); hold on
plot(t, imp1_r_save(:,1), 'r--'); hold on
plot(t, imp1_l_save(:,2), 'g'); hold on
plot(t, imp1_r_save(:,2), 'g--'); hold on
plot(t, imp1_l_save(:,3), 'b'); hold on
plot(t, imp1_r_save(:,3), 'b--'); hold on
legend('Errors','Forces');
title('Impedance equation (translation)');

figure();
plot(t, imp1_l_save(:,4), 'r'); hold on
plot(t, imp1_r_save(:,4), 'r--'); hold on
plot(t, imp1_l_save(:,5), 'g'); hold on
plot(t, imp1_r_save(:,5), 'g--'); hold on
plot(t, imp1_l_save(:,6), 'b'); hold on
plot(t, imp1_r_save(:,6), 'b--'); hold on
legend('Errors','Forces');
title('Impedance equation (rotation)');
%%%%%%%%%%%%%
