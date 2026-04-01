function atm = F16_Atmosphere(h)
% ISA standard atmosphere model. Valid range: -2000 m to 86000 m.

R_air  = 287.058;   % [J/(kg.K)]
gamma  = 1.4;
g0     = 9.80665;   % [m/s^2]

% Sutherland's law constants for dry air
mu_ref = 1.716e-5;  % [Pa.s]
T_ref  = 273.15;    % [K]
S_suth = 110.4;     % [K]

% ISA layers: [base altitude (m), base temp (K), base pressure (Pa), lapse rate (K/m)]
atm_layers = [
        0,    288.15,    101325.000,   -0.0065;
    11000,    216.65,     22632.100,    0.0000;
    20000,    216.65,      5474.890,    0.0010;
    32000,    228.65,       868.019,    0.0028;
    47000,    270.65,       110.906,    0.0000;
    51000,    270.65,        66.939,   -0.0028;
    71000,    214.65,         3.956,   -0.0020;
];

layer_bases = atm_layers(:, 1);
layer_idx   = find(layer_bases <= h, 1, 'last');
if isempty(layer_idx),                   layer_idx = 1; end
if layer_idx > size(atm_layers, 1),      layer_idx = size(atm_layers, 1); end

h_b = atm_layers(layer_idx, 1);
T_b = atm_layers(layer_idx, 2);
P_b = atm_layers(layer_idx, 3);
L_b = atm_layers(layer_idx, 4);

T = T_b + L_b * (h - h_b);

if abs(L_b) > 1e-10
    % Gradient layer: P = P_b * (T_b/T)^(g0 / R*L)
    P = P_b * (T_b / T)^(g0 / (R_air * L_b));
else
    % Isothermal layer: P = P_b * exp(-g0*(h-h_b) / R*T_b)
    P = P_b * exp(-g0 * (h - h_b) / (R_air * T_b));
end

rho = P / (R_air * T);
a   = sqrt(gamma * R_air * T);
mu  = mu_ref * (T/T_ref)^1.5 * (T_ref + S_suth) / (T + S_suth);
nu  = mu / rho;

atm.h    = h;
atm.h_ft = h * 3.28084;
atm.T    = T;
atm.T_C  = T - 273.15;
atm.p    = P;
atm.rho  = rho;
atm.a    = a;
atm.mu   = mu;
atm.nu   = nu;
end
