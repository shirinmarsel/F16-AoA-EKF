function sens = F16_Sensors(x, xdot, params)
% F-16 sensor models with additive Gaussian noise and constant bias.
% Each sensor: measurement = truth + bias + sigma*randn

u     = x(1);   v = x(2);   w  = x(3);
p     = x(4);   q = x(5);   r  = x(6);
phi   = x(7);   theta = x(8);
h     = x(12);

udot = xdot(1);  vdot = xdot(2);  wdot = xdot(3);

% Sensor noise parameters
sig_accel  = 0.02;   bias_accel = 0.01;   % accelerometer [m/s^2]
sig_gyro   = 0.003;  bias_gyro  = 0.005;  % gyroscope     [rad/s]
sig_V      = 0.3;    bias_V     = 0.5;    % airspeed      [m/s]
sig_alpha  = 0.2;    bias_alpha = 0.1;    % AoA vane      [deg]
sig_h      = 1.5;    bias_h     = 2.0;    % barometer     [m]

g0 = 9.80665;
gx = -g0 * sin(theta);
gy =  g0 * cos(theta) * sin(phi);
gz =  g0 * cos(theta) * cos(phi);

% Specific force = non-gravitational force / mass = Fx/m
% Derived from EOM: udot = (r*v - q*w) + Fx/m + gx  =>  Fx/m = udot - (r*v - q*w) - gx
ax_true = udot - (r*v - q*w) - gx;
ay_true = vdot - (p*w - r*u) - gy;
az_true = wdot - (q*u - p*v) - gz;

sens.ax = ax_true + bias_accel + sig_accel * randn;
sens.ay = ay_true + bias_accel + sig_accel * randn;
sens.az = az_true + bias_accel + sig_accel * randn;

sens.p = p + bias_gyro + sig_gyro * randn;
sens.q = q + bias_gyro + sig_gyro * randn;
sens.r = r + bias_gyro + sig_gyro * randn;

V_true = sqrt(u^2 + v^2 + w^2);
sens.V = V_true + bias_V + sig_V * randn;

alpha_true_deg = rad2deg(atan2(w, u));
sens.alpha     = alpha_true_deg + bias_alpha + sig_alpha * randn;

sens.h = h + bias_h + sig_h * randn;

% Truth values stored for post-processing comparison
sens.ax_true    = ax_true;
sens.ay_true    = ay_true;
sens.az_true    = az_true;
sens.p_true     = p;
sens.q_true     = q;
sens.r_true     = r;
sens.V_true     = V_true;
sens.alpha_true = alpha_true_deg;
sens.h_true     = h;
end
