% F-16 EKF simulation — AoA estimation.

clc; clear; close all;

fprintf('F-16 EKF Simulation\n\n');

dt    = 0.01;   % 100 Hz
t_end = 30;
t_vec = 0:dt:t_end;
N     = length(t_vec);

% Trim initial conditions (from F16_Trim)
V0     = 100;
alpha0 = deg2rad(8.2816);
theta0 = alpha0;
u0     = V0 * cos(alpha0);
w0     = V0 * sin(alpha0);
x0     = [u0; 0; w0; 0; 0; 0; 0; theta0; 0; 0; 0; 6000];
u_ctrl = [-0.3957; 0; 0; 6279.17];

% Truth simulation
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~, X_truth] = ode45(@(t,x) F16_EOM(t,x,u_ctrl,[]), t_vec, x0, opts);
fprintf('Truth model complete. %d steps.\n\n', N);

% EKF initialisation
ekf = F16_EKF_Init(x0, dt);
fprintf('\n');

alpha_truth_hist = zeros(N,1);
alpha_ekf_hist   = zeros(N,1);
alpha_sens_hist  = zeros(N,1);
innov_hist       = zeros(N,6);
P_diag_hist      = zeros(N,9);

fprintf('%-8s %-12s %-12s %-12s\n','Time','Truth AoA','EKF AoA','P trace');
fprintf('%s\n', repmat('-',1,48));

for i = 1:N
    x_true    = X_truth(i,:)';
    xdot_true = F16_EOM(t_vec(i), x_true, u_ctrl, []);
    sens      = F16_Sensors(x_true, xdot_true, []);

    % Measurement vector — order must match EKF_meas z_pred: [p, q, r, ax, V, az]
    z = [sens.p; sens.q; sens.r; sens.ax; sens.V; sens.az];

    ekf = F16_EKF(ekf, z, u_ctrl, dt);

    alpha_truth_hist(i) = sens.alpha_true;
    alpha_ekf_hist(i)   = ekf.alpha;
    alpha_sens_hist(i)  = sens.alpha;
    innov_hist(i,:)     = ekf.innov';
    P_diag_hist(i,:)    = ekf.P_diag';

    if mod(i-1, 500) == 0
        fprintf('  %-8.1f %-12.4f %-12.4f %-12.6f\n', ...
                t_vec(i), alpha_truth_hist(i), alpha_ekf_hist(i), trace(ekf.P));
    end
end

% Ignore first 2 s transient for RMS computation
idx      = t_vec > 2.0;
err_sens = alpha_sens_hist - alpha_truth_hist;
err_ekf  = alpha_ekf_hist  - alpha_truth_hist;

fprintf('\n========= RESULTS (after t=2s) =========\n');
fprintf('Sensor RMS error : %.4f deg\n', rms(err_sens(idx)));
fprintf('EKF    RMS error : %.4f deg\n', rms(err_ekf(idx)));
if rms(err_sens(idx)) > 0
    fprintf('Improvement      : %.1f%%\n', ...
            (1 - rms(err_ekf(idx))/rms(err_sens(idx)))*100);
end
fprintf('=========================================\n\n');

figure('Name','F-16 EKF Results','Position',[50 50 1200 800]);

subplot(2,2,1);
plot(t_vec, alpha_truth_hist, 'b-', 'LineWidth', 2); hold on;
plot(t_vec, alpha_sens_hist,  'r.', 'MarkerSize', 3);
plot(t_vec, alpha_ekf_hist,   'g-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('\alpha [deg]');
title('Angle of Attack Estimation');
legend('Truth','Sensor','EKF','Location','best');
ylim([6 11]); grid on;

subplot(2,2,2);
plot(t_vec, err_sens, 'r-', 'LineWidth', 1); hold on;
plot(t_vec, err_ekf,  'g-', 'LineWidth', 2);
yline(0,'k--');
xlabel('Time [s]'); ylabel('Error [deg]');
title('AoA Estimation Error');
legend('Sensor error','EKF error','Location','best');
ylim([-1.5 1.5]); grid on;

subplot(2,2,3);
plot(t_vec, P_diag_hist(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_vec, P_diag_hist(:,3), 'r-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Variance');
title('State Uncertainty (P diagonal)');
legend('P_{uu}','P_{ww}','Location','best'); grid on;

subplot(2,2,4);
plot(t_vec, innov_hist(:,2), 'b-', 'LineWidth', 1);
yline(0,'k--');
xlabel('Time [s]'); ylabel('Innovation [rad/s]');
title('Pitch Gyro Innovation (healthy = zero-mean noise)');
grid on;

sgtitle('F-16 EKF — AoA Estimation Results','FontWeight','bold');