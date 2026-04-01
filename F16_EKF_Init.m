function ekf = F16_EKF_Init(x0_truth, dt)
% Initialise the 9-state EKF for F-16 AoA estimation.
%
% STATE VECTOR (9 states):
%   [u, v, w, p, q, r, phi, theta, psi]
%
% Design choice — 9 states, not 15:
%   Gyro and accel biases are unobservable in straight and level trimmed flight.
%   Including unobservable states causes EKF divergence.
%   Sensor biases are absorbed into a conservatively tuned R instead.
%
% AoA is derived post-hoc: alpha = atan2(w, u) — not a filter state.

n = 9;

% Initial state — 1% perturbation on u, w, theta to simulate realistic init uncertainty
ekf.x_hat    = zeros(n, 1);
ekf.x_hat(1) = x0_truth(1) * 1.01;   % u
ekf.x_hat(2) = x0_truth(2);           % v
ekf.x_hat(3) = x0_truth(3) * 1.01;   % w
ekf.x_hat(4) = x0_truth(4);           % p
ekf.x_hat(5) = x0_truth(5);           % q
ekf.x_hat(6) = x0_truth(6);           % r
ekf.x_hat(7) = x0_truth(7);           % phi
ekf.x_hat(8) = x0_truth(8) * 1.01;   % theta
ekf.x_hat(9) = x0_truth(9);           % psi

% Initial covariance P0 — diagonal, set to expected initial error squared
P0_diag = [
    (2.0)^2;    % u     [m/s]
    (1.0)^2;    % v     [m/s]
    (1.0)^2;    % w     [m/s]
    (0.05)^2;   % p     [rad/s]
    (0.05)^2;   % q     [rad/s]
    (0.05)^2;   % r     [rad/s]
    (0.02)^2;   % phi   [rad]
    (0.02)^2;   % theta [rad]
    (0.05)^2;   % psi   [rad]
];
ekf.P = diag(P0_diag);

% Process noise Q — represents unmodelled turbulence and aerodynamic disturbances.
% Scaled by dt (continuous-to-discrete: Q_d = sigma^2 * dt).
% Small values reflect high trust in the physics model.
Q_diag = [
    (0.5)^2  * dt;    % u
    (0.5)^2  * dt;    % v
    (0.5)^2  * dt;    % w
    (0.05)^2 * dt;    % p
    (0.05)^2 * dt;    % q
    (0.05)^2 * dt;    % r
    (0.005)^2 * dt;   % phi
    (0.005)^2 * dt;   % theta
    (0.005)^2 * dt;   % psi
];
ekf.Q = diag(Q_diag);

% Measurement noise R — z = [p, q, r, ax, V, az]
% Design choice: R set larger than true sensor noise to make filter conservative
% (trust physics model over sensors). Sensor biases are implicitly absorbed here.
R_diag = [
    (0.05)^2;   % p    [rad/s]
    (0.05)^2;   % q    [rad/s]
    (0.05)^2;   % r    [rad/s]
    (0.1)^2;    % ax   [m/s^2]
    (1.0)^2;    % V    [m/s]   — conservatively large; pitot noise ~0.3 m/s
    (0.1)^2;    % az   [m/s^2]
];
ekf.R = diag(R_diag);

ekf.dt     = dt;
ekf.n      = n;
ekf.alpha  = rad2deg(atan2(ekf.x_hat(3), ekf.x_hat(1)));
ekf.K      = zeros(n, 6);
ekf.innov  = zeros(6, 1);
ekf.P_diag = diag(ekf.P);

fprintf('EKF initialised (9-state, 6-measurement).\n');
fprintf('  Initial alpha : %.4f deg\n', ekf.alpha);
fprintf('  P trace       : %.4f\n',     trace(ekf.P));
end
