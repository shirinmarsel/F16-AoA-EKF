% Test F16 sensor models

clc; clear; close all;

% Run truth simulation first
% Use trim conditions
V0     = 100;
alpha0 = deg2rad(8.2816);
theta0 = alpha0;
u0     = V0 * cos(alpha0);
v0     = 0;
w0     = V0 * sin(alpha0);
x0     = [u0;v0;w0; 0;0;0; 0;theta0;0; 0;0;6000];
u_ctrl = [-0.3957; 0; 0; 6279.17];

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t, X] = ode45(@(t,x) F16_EOM(t,x,u_ctrl,[]), [0,30], x0, opts);

% Generate sensor measurements at each time step
N = length(t);

% Preallocate
alpha_true_hist  = zeros(N,1);
alpha_sens_hist  = zeros(N,1);
q_true_hist      = zeros(N,1);
q_sens_hist      = zeros(N,1);
h_true_hist      = zeros(N,1);
h_sens_hist      = zeros(N,1);
V_true_hist      = zeros(N,1);
V_sens_hist      = zeros(N,1);

for i = 1:N
    % Get truth derivatives at this point
    xdot_i = F16_EOM(t(i), X(i,:)', u_ctrl, []);
    
    % Get sensor measurements
    sens_i = F16_Sensors(X(i,:)', xdot_i, []);
    
    % Store
    alpha_true_hist(i) = sens_i.alpha_true;
    alpha_sens_hist(i) = sens_i.alpha;
    q_true_hist(i)     = sens_i.q_true;
    q_sens_hist(i)     = sens_i.q;
    h_true_hist(i)     = sens_i.h_true;
    h_sens_hist(i)     = sens_i.h;
    V_true_hist(i)     = sens_i.V_true;
    V_sens_hist(i)     = sens_i.V;
end

% Plot truth vs measurements
figure('Name','F-16 Sensor Models','Position',[100 100 1000 700]);

subplot(2,2,1);
plot(t, alpha_true_hist, 'b-', 'LineWidth', 2); hold on;
plot(t, alpha_sens_hist, 'r.', 'MarkerSize', 4);
xlabel('Time [s]'); ylabel('\alpha [deg]');
title('AoA: Truth vs Sensor');
legend('Truth','Sensor','Location','best');
grid on;

subplot(2,2,2);
plot(t, rad2deg(q_true_hist), 'b-', 'LineWidth', 2); hold on;
plot(t, rad2deg(q_sens_hist), 'r.', 'MarkerSize', 4);
xlabel('Time [s]'); ylabel('q [deg/s]');
title('Pitch Rate: Truth vs Sensor');
legend('Truth','Sensor','Location','best');
grid on;

subplot(2,2,3);
plot(t, h_true_hist, 'b-', 'LineWidth', 2); hold on;
plot(t, h_sens_hist, 'r.', 'MarkerSize', 4);
xlabel('Time [s]'); ylabel('h [m]');
title('Altitude: Truth vs Sensor');
legend('Truth','Sensor','Location','best');
grid on;

subplot(2,2,4);
plot(t, V_true_hist, 'b-', 'LineWidth', 2); hold on;
plot(t, V_sens_hist, 'r.', 'MarkerSize', 4);
xlabel('Time [s]'); ylabel('V [m/s]');
title('Airspeed: Truth vs Sensor');
legend('Truth','Sensor','Location','best');
grid on;

sgtitle('F-16 Brick 4 — Sensor Models: Truth vs Measurements', ...
        'FontWeight','bold');

fprintf('Sensor test complete.\n');
fprintf('Expected: red dots scattered around blue line\n');
fprintf('AoA bias: %.1f deg  noise: %.1f deg\n', 0.1, 0.2);
fprintf('Alt bias: %.1f m    noise: %.1f m\n',   2.0, 1.5);