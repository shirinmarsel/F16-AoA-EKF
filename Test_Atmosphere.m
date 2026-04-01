clc;    % clear command window
clear;  % clear all variables

fprintf('=Testing F16_Atmosphere= \n\n');

% Test 1 - Sea level
atm = F16_Atmosphere(0);
fprintf('Sea Level (h = 0 m): \n');
fprintf('Temperature : %.2f K   \n', atm.T);
fprintf('Pressure    : %.1f Pa  \n', atm.p);
fprintf('Density     : %.4f kg/m3 \n', atm.rho);
fprintf('Speed sound : %.2f m/s \n\n', atm.a);

% Test 2 - Tropopause base
atm = F16_Atmosphere(11000);
fprintf('Tropopause (h = 11000 m): \n');
fprintf('Temperature : %.2f K   \n', atm.T);
fprintf('Pressure    : %.1f Pa  \n', atm.p);
fprintf('Density     : %.4f kg/m3 \n', atm.rho);
fprintf('Speed sound : %.2f m/s \n\n', atm.a);

% Test 3 - F-16 typical cruise altitude (~40,000 ft = 12,192 m)
atm = F16_Atmosphere(12192);
fprintf('F-16 Cruise (h = 12192 m / 40000 ft): \n');
fprintf('Temperature : %.2f K\n', atm.T);
fprintf('Pressure    : %.1f Pa\n', atm.p);
fprintf('Density     : %.4f kg/m3\n', atm.rho);
fprintf(['Mach 1 speed: %.2f m/s \n\n'], atm.a);

fprintf('=Done= \n');