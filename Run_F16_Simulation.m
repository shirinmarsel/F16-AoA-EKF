% F-16 6-DoF truth model simulation.
% Straight and level flight: theta = alpha (no climb), phi=0, beta=0.

clc; clear; close all;

V0     = 100;                   % [m/s]
alpha0 = deg2rad(8.2816);       % [rad]  — from trim solution
theta0 = alpha0;                % level flight assumption

u0  = V0 * cos(alpha0);
v0  = 0;
w0  = V0 * sin(alpha0);
h0  = 6000;   % [m]

x0 = [u0; v0; w0; 0; 0; 0; 0; theta0; 0; 0; 0; h0];

% Trim control inputs — from F16_Trim solution
u_ctrl = [-0.3957; 0; 0; 6279.17];   % [de_deg; da_deg; dr_deg; Thrust_N]

opts    = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, X]  = ode45(@(t,x) F16_EOM(t, x, u_ctrl, []), [0, 30], x0, opts);

u_hist     = X(:,1);
w_hist     = X(:,3);
theta_hist = X(:,8);
h_hist     = X(:,12);
V_hist     = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2);

alpha_truth = rad2deg(atan2(w_hist, u_hist));

figure('Name','F-16 6-DoF Truth Simulation','Position',[100 100 1000 700]);

subplot(2,2,1);
plot(t, alpha_truth, 'b-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('\alpha [deg]'); title('Angle of Attack (Truth)'); grid on;

subplot(2,2,2);
plot(t, rad2deg(theta_hist), 'r-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('\theta [deg]'); title('Pitch Angle'); grid on;

subplot(2,2,3);
plot(t, h_hist, 'g-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Altitude [m]'); title('Altitude'); grid on;

subplot(2,2,4);
plot(t, V_hist, 'm-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('V [m/s]'); title('Airspeed'); grid on;

sgtitle('F-16 6-DoF Truth Model', 'FontWeight', 'bold');

fprintf('Simulation complete.\n');
fprintf('Final altitude  : %.1f m\n',   h_hist(end));
fprintf('Final AoA       : %.2f deg\n', alpha_truth(end));
fprintf('Final airspeed  : %.1f m/s\n', V_hist(end));
