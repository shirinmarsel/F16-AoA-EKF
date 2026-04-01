function xdot = F16_EOM(t, x, u_ctrl, params)

x = x(:);

% State: [u v w p q r phi theta psi x_E y_E h]
u = x(1);  v = x(2);  w = x(3);
p = x(4);  q = x(5);  r = x(6);
phi = x(7);  theta = x(8);  psi = x(9);
x_E = x(10);  y_E = x(11);  h = x(12);
h = max(h, 0);  % floor altitude to keep atmosphere model valid

de = u_ctrl(1);  da = u_ctrl(2);  dr = u_ctrl(3);  dT = u_ctrl(4);

% Mass properties — NASA F-16 Model, Table 2, p.18
m    = 9295.44;   % [kg]
Ixx  = 12874.8;   % [kg.m^2]
Iyy  = 75673.6;   % [kg.m^2]
Izz  = 85552.1;   % [kg.m^2]
Ixz  = 1331.4;    % [kg.m^2]  — Ixz coupling retained in moment equations
g0   = 9.80665;   % [m/s^2]

% Reference geometry
S    = 27.87;   % wing area  [m^2]
b    = 9.144;   % wingspan   [m]
cbar = 3.450;   % MAC        [m]

% Inertia coupling constants (Stevens & Lewis formulation)
Gamma = Ixx*Izz - Ixz^2;
c1 = ((Iyy - Izz)*Izz - Ixz^2) / Gamma;
c2 = ((Ixx - Iyy + Izz)*Ixz)   / Gamma;
c3 = Izz / Gamma;
c4 = Ixz / Gamma;
c5 = (Izz - Ixx) / Iyy;
c6 = Ixz / Iyy;
c7 = 1 / Iyy;
c8 = ((Ixx - Iyy)*Ixx + Ixz^2) / Gamma;
c9 = Ixx / Gamma;

atm   = F16_Atmosphere(h);
V     = sqrt(u^2 + v^2 + w^2);
alpha = atan2(w, u);        % [rad]
beta  = asin(v / V);        % [rad]
qbar  = 0.5 * atm.rho * V^2;

aero = F16_AeroTables(rad2deg(alpha), rad2deg(beta), de, da, dr);

% Aerodynamic forces and moments — rolling/yawing use b, pitching uses cbar
Fa_x = qbar * S * aero.Cx;
Fa_y = qbar * S * aero.Cy;
Fa_z = qbar * S * aero.Cz;
La   = qbar * S * b    * aero.Cl;
Ma   = qbar * S * cbar * aero.Cm;
Na   = qbar * S * b    * aero.Cn;

% Thrust acts along body x-axis only
Fx = Fa_x + dT;
Fy = Fa_y;
Fz = Fa_z;

% Gravity projected into body frame
gx = -g0 * sin(theta);
gy =  g0 * cos(theta) * sin(phi);
gz =  g0 * cos(theta) * cos(phi);

% Force equations
udot = r*v - q*w + Fx/m + gx;
vdot = p*w - r*u + Fy/m + gy;
wdot = q*u - p*v + Fz/m + gz;

% Moment equations
pdot = (c1*r + c2*p)*q + c3*La + c4*Na;
qdot =  c5*p*r - c6*(p^2 - r^2) + c7*Ma;
rdot = (c8*p - c2*r)*q + c4*La + c9*Na;

% Kinematic equations — 3-2-1 Euler; singularity at theta = ±90 deg (not reached in level flight)
phidot   = p + (q*sin(phi) + r*cos(phi))*tan(theta);
thetadot = q*cos(phi) - r*sin(phi);
psidot   = (q*sin(phi) + r*cos(phi)) / cos(theta);

% Navigation equations — full DCM body-to-Earth rotation
xEdot =  u*cos(theta)*cos(psi) ...
       + v*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) ...
       + w*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));

yEdot =  u*cos(theta)*sin(psi) ...
       + v*(sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)) ...
       + w*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));

hdot  =  u*sin(theta) ...
       - v*sin(phi)*cos(theta) ...
       - w*cos(phi)*cos(theta);

xdot = [udot; vdot; wdot; ...
        pdot; qdot; rdot; ...
        phidot; thetadot; psidot; ...
        xEdot; yEdot; hdot];
end
