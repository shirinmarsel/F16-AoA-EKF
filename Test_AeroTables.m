clc; clear;

% Test at F-16 trim condition (approximately)
% Straight and level flight at moderate speed:
% alpha ~ 2 deg, beta = 0, elevator ~ -2 deg

aero = F16_AeroTables(2.0, 0.0, -2.0, 0.0, 0.0);

fprintf('Aero Coefficients at alpha=2, beta=0, de=-2:\n');
fprintf('Cx = %.4f\n', aero.Cx);
fprintf('Cy = %.4f  (should be 0 - symmetric flight)\n', aero.Cy);
fprintf('Cz = %.4f\n', aero.Cz);
fprintf('Cl = %.4f  (should be 0 - symmetric flight)\n', aero.Cl);
fprintf('Cm = %.4f\n', aero.Cm);
fprintf('Cn = %.4f  (should be 0 - symmetric flight)\n', aero.Cn);