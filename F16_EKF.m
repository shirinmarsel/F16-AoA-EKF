function ekf = F16_EKF(ekf, z, u_ctrl, dt)
% Extended Kalman Filter for F-16 AoA estimation.
% Called once per time step: PREDICT → UPDATE → extract AoA.
%
% Inputs:
%   z = [p; q; r; ax; V; az]  measurement vector [6x1]
%   u_ctrl = [de; da; dr; dT]
%
% Output: ekf.alpha = estimated AoA [deg]

n      = ekf.n;
n_meas = 6;

% --- PREDICT ---
% Euler integration: x^- = x + f(x,u)*dt
f0          = EKF_dynamics(ekf.x_hat, u_ctrl);
x_hat_minus = ekf.x_hat + f0 * dt;

% Jacobian via numerical finite difference; discretised as F_d = I + F_c*dt
eps   = 1e-6;
F_jac = zeros(n, n);
for j = 1:n
    x_p        = ekf.x_hat;
    x_p(j)     = x_p(j) + eps;
    F_jac(:,j) = (EKF_dynamics(x_p, u_ctrl) - f0) / eps;
end
F_d = eye(n) + F_jac * dt;

P_minus = F_d * ekf.P * F_d' + ekf.Q;

% --- UPDATE ---
[z_pred, H] = EKF_meas(x_hat_minus, u_ctrl, n_meas);
innovation  = z - z_pred;

S = H * P_minus * H' + ekf.R;
K = (P_minus * H') / S;     % using / instead of inv() for numerical stability
K = max(min(K, 2.0), -2.0); % clip gain to prevent wild corrections on filter startup

ekf.x_hat = x_hat_minus + K * innovation;

% Joseph form — more numerically stable than simple (I-KH)*P
IKH   = eye(n) - K * H;
ekf.P = IKH * P_minus * IKH' + K * ekf.R * K';
ekf.P = 0.5 * (ekf.P + ekf.P');  % enforce symmetry

% Enforce minimum diagonal variance to prevent filter death
min_var = [1e-4; 1e-4; 1e-4; 1e-7; 1e-7; 1e-7; 1e-7; 1e-7; 1e-7];
for k = 1:n
    if ekf.P(k,k) < min_var(k)
        ekf.P(k,k) = min_var(k);
    end
end

% --- AoA ESTIMATE ---
% AoA is derived from estimated velocity states, not estimated directly
ekf.alpha  = rad2deg(atan2(ekf.x_hat(3), ekf.x_hat(1)));

ekf.K      = K;
ekf.innov  = innovation;
ekf.P_diag = diag(ekf.P);
end


function xdot = EKF_dynamics(x_hat, u_ctrl)
% 9-state dynamics model used inside the EKF predict step.
% Assumption: atmosphere evaluated at fixed h_ref = 6000 m (nominal cruise altitude).
% Justified for short (30 s) near-trim simulation where altitude variation is small.

u   = x_hat(1);  v = x_hat(2);  w = x_hat(3);
p   = x_hat(4);  q = x_hat(5);  r = x_hat(6);
phi = x_hat(7);  theta = x_hat(8);

de = u_ctrl(1);  da = u_ctrl(2);  dr = u_ctrl(3);  dT = u_ctrl(4);

m    = 9295.44;  Ixx = 12874.8;  Iyy = 75673.6;
Izz  = 85552.1;  Ixz = 1331.4;   g0  = 9.80665;
S    = 27.87;    b   = 9.144;     cbar = 3.450;

Gamma = Ixx*Izz - Ixz^2;
c1 = ((Iyy-Izz)*Izz - Ixz^2) / Gamma;
c2 = ((Ixx-Iyy+Izz)*Ixz)     / Gamma;
c3 = Izz / Gamma;   c4 = Ixz / Gamma;
c5 = (Izz-Ixx) / Iyy;  c6 = Ixz / Iyy;  c7 = 1/Iyy;
c8 = ((Ixx-Iyy)*Ixx + Ixz^2) / Gamma;   c9 = Ixx / Gamma;

h_ref = 6000;
atm   = F16_Atmosphere(h_ref);
V     = max(sqrt(u^2 + v^2 + w^2), 1.0);

alpha_deg = rad2deg(atan2(w, u));
beta_deg  = rad2deg(asin(max(min(v/V, 1), -1)));
qbar      = 0.5 * atm.rho * V^2;
aero      = F16_AeroTables(alpha_deg, beta_deg, de, da, dr);

Fx = qbar*S*aero.Cx + dT;
Fy = qbar*S*aero.Cy;
Fz = qbar*S*aero.Cz;
La = qbar*S*b*aero.Cl;
Ma = qbar*S*cbar*aero.Cm;
Na = qbar*S*b*aero.Cn;

gx = -g0*sin(theta);
gy =  g0*cos(theta)*sin(phi);
gz =  g0*cos(theta)*cos(phi);

udot     = r*v - q*w + Fx/m + gx;
vdot     = p*w - r*u + Fy/m + gy;
wdot     = q*u - p*v + Fz/m + gz;
pdot     = (c1*r + c2*p)*q + c3*La + c4*Na;
qdot     =  c5*p*r - c6*(p^2-r^2) + c7*Ma;
rdot     = (c8*p - c2*r)*q + c4*La + c9*Na;
phidot   = p + (q*sin(phi) + r*cos(phi))*tan(theta);
thetadot = q*cos(phi) - r*sin(phi);
psidot   = (q*sin(phi) + r*cos(phi)) / cos(theta);

xdot = [udot; vdot; wdot; pdot; qdot; rdot; phidot; thetadot; psidot];
end


function [z_pred, H] = EKF_meas(x_hat, u_ctrl, n_meas)
% Predicted measurement function h(x) and numerical Jacobian H [6x9].
% Measurement vector: z = [p, q, r, ax, V, az]
%   ax = Fx/m  (specific force x — non-gravitational only)
%   az = Fz/m  (specific force z)
% Assumption: atmosphere at fixed h_ref = 6000 m (same as EKF_dynamics).

u   = x_hat(1);  v = x_hat(2);  w = x_hat(3);
p   = x_hat(4);  q = x_hat(5);  r = x_hat(6);

de = u_ctrl(1);  da = u_ctrl(2);  dr = u_ctrl(3);  dT = u_ctrl(4);

m    = 9295.44;  S = 27.87;  b = 9.144;  cbar = 3.450;

h_ref = 6000;
atm   = F16_Atmosphere(h_ref);
V     = max(sqrt(u^2 + v^2 + w^2), 1.0);

alpha_deg = rad2deg(atan2(w, u));
beta_deg  = rad2deg(asin(max(min(v/V,1),-1)));
qbar      = 0.5 * atm.rho * V^2;
aero      = F16_AeroTables(alpha_deg, beta_deg, de, da, dr);

Fx = qbar*S*aero.Cx + dT;
Fz = qbar*S*aero.Cz;

z_pred = [p; q; r; Fx/m; V; Fz/m];

% Numerical H Jacobian [6x9]
n_states = length(x_hat);
H        = zeros(n_meas, n_states);
eps      = 1e-6;
for j = 1:n_states
    x_p    = x_hat;
    x_p(j) = x_p(j) + eps;

    u_p = x_p(1);  v_p = x_p(2);  w_p = x_p(3);
    p_p = x_p(4);  q_p = x_p(5);  r_p = x_p(6);

    V_p     = max(sqrt(u_p^2+v_p^2+w_p^2), 1.0);
    alpha_p = rad2deg(atan2(w_p, u_p));
    beta_p  = rad2deg(asin(max(min(v_p/V_p,1),-1)));
    qbar_p  = 0.5 * atm.rho * V_p^2;
    aero_p  = F16_AeroTables(alpha_p, beta_p, de, da, dr);

    Fx_p = qbar_p*S*aero_p.Cx + dT;
    Fz_p = qbar_p*S*aero_p.Cz;

    H(:,j) = ([p_p; q_p; r_p; Fx_p/m; V_p; Fz_p/m] - z_pred) / eps;
end
end
