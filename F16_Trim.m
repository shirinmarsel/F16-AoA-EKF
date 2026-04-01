% Find F-16 trim condition for straight and level flight.
% Solves for [alpha, de, Thrust] given fixed V and h.
% Assumptions: beta=0, p=q=r=0, v=0, theta=alpha (no climb angle).

clc; clear;

V_trim = 100;    % [m/s]
h_trim = 6000;   % [m]

x_guess = [0.0; -2.0; 4000];   % [alpha_deg, de_deg, Thrust_N]

opts   = optimset('Display','iter','TolFun',1e-10,'TolX',1e-10);
x_trim = fsolve(@(x) trim_residual(x, V_trim, h_trim), x_guess, opts);

theta_trim_rad = deg2rad(x_trim(1));  % theta = alpha for level flight [rad]

fprintf('\n======= TRIM SOLUTION =======\n');
fprintf('Airspeed  : %.2f m/s\n',  V_trim);
fprintf('Altitude  : %.0f m\n',    h_trim);
fprintf('Alpha     : %.4f deg\n',  x_trim(1));
fprintf('Elevator  : %.4f deg\n',  x_trim(2));
fprintf('Thrust    : %.2f N\n',    x_trim(3));
fprintf('Theta     : %.4f deg\n',  rad2deg(theta_trim_rad));

res = trim_residual(x_trim, V_trim, h_trim);
fprintf('\nResiduals (should all be ~0):\n');
fprintf('  u_dot : %e\n', res(1));
fprintf('  w_dot : %e\n', res(2));
fprintf('  q_dot : %e\n', res(3));
fprintf('=============================\n');


function residual = trim_residual(x, V, h)
% Equilibrium residuals for straight and level flight.
% Assumes: beta=0, p=q=r=0, v=0, theta=alpha.

alpha_deg = x(1);
de_deg    = x(2);
Thrust    = x(3);

m    = 9295.44;  g0 = 9.80665;
S    = 27.87;    cbar = 3.450;  Iyy = 75673.6;

alpha_rad = deg2rad(alpha_deg);
theta     = alpha_rad;

u = V * cos(alpha_rad);
w = V * sin(alpha_rad);

atm  = F16_Atmosphere(h);
qbar = 0.5 * atm.rho * V^2;
aero = F16_AeroTables(alpha_deg, 0, de_deg, 0, 0);

Fx = qbar * S * aero.Cx + Thrust;
Fz = qbar * S * aero.Cz;
M  = qbar * S * cbar * aero.Cm;

gx = -g0 * sin(theta);
gz =  g0 * cos(theta);   % phi=0 assumed (wings level)

residual(1) = Fx/m + gx;   % u_dot = 0
residual(2) = Fz/m + gz;   % w_dot = 0
residual(3) = M / Iyy;     % q_dot = 0
end
