% Input parameters
L = 1; % Length of the panel (m)
R = 5; % Radius of curvature (m)
t_skin = 0.001; % Thickness of the skins (m)
t_core = 0.01; % Thickness of the core (m)
E_skin = 70e9; % Young's modulus of the skins (Pa)
E_core = 2e9; % Young's modulus of the core (Pa)
v_skin = 0.33; % Poisson's ratio of the skins
v_core = 0.1; % Poisson's ratio of the core
alpha = 23e-6; % Coefficient of thermal expansion (1/K)
T = 100; % Temperature change (K)

% Calculation of the thermal stress distribution
[~, ~, ~, ~, sigma_xx, ~, ~] = thermal_stress(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, alpha, T);

% Plotting of the thermal stress distribution
x = linspace(0, L, 100); % Position along the length of the panel (m)
y = linspace(-R, R, 100); % Position along the width of the panel (m)
[X, Y] = meshgrid(x, y);
figure;
contourf(X, Y, sigma_xx/1e6, 50);
colorbar;
title('Thermal Stress Distribution (MPa)');
xlabel('Length (m)');
ylabel('Width (m)');

% Calculation of the natural frequencies and mode shapes
N = 10; % Number of terms in the approximation
[u_n, w_n, omega] = curved_panel_vibration(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, N);

% Plotting of the first mode shape
x = linspace(0, L, 100); % Position along the length of the panel (m)
u = u_n(1)*sin(pi*x/L);
w = w_n(1)*sin(pi*x/L);
figure;
subplot(2,1,1);
plot(x, u);
title('Displacement of the First Mode Shape');
xlabel('Length (m)');
ylabel('Displacement (m)');
subplot(2,1,2);
plot(x, w);
title('Deflection of the First Mode Shape');
xlabel('Length (m)');
ylabel('Deflection (m)');


function [u, w, phi, kappa, sigma_xx, sigma_yy, tau_xy] = thermal_stress(L, R,t_skin, t_core, E_skin, E_core, v_skin, v_core, alpha, T)
% Calculation of the thermal stress distribution in a curved sandwich panel beam
% Inputs:
% L - Length of the panel (m)
% R - Radius of curvature (m)
% t_skin - Thickness of the skins (m)
% t_core - Thickness of the core (m)
% E_skin - Young's modulus of the skins (Pa)
% E_core - Young's modulus of the core (Pa)
% v_skin - Poisson's ratio of the skins
% v_core - Poisson's ratio of the core
% alpha - Coefficient of thermal expansion (1/K)
% T - Temperature change (K)
% Outputs:
% u - Displacement in the x-direction (m)
% w - Displacement in the y-direction (m)
% phi - Rotation (rad)
% kappa - Curvature (1/m)
% sigma_xx - Stress in the x-direction (Pa)
% sigma_yy - Stress in the y-direction (Pa)
% tau_xy - Shear stress (Pa)

% Calculation of the moment of inertia and area of the panel
I = (pi/4)*((R+t_skin)^4 - R^4 - (R-t_skin)^4) + ...
    t_core*(pi/4)*((R+t_skin+t_core)^4 - (R+t_skin)^4 - (R-t_skin)^4 + (R-t_skin-t_core)^4);
A = (pi/4)*((R+t_skin)^2 - (R-t_skin)^2) + t_core*(pi/4)*((R+t_skin+t_core)^2 - (R+t_skin)^2 - (R-t_skin)^2 + (R-t_skin-t_core)^2);

% Calculation of the thermal strain
epsilon_T = alpha*T;

% Calculation of the thermal stress in the skins
sigma_T_skin = E_skin*epsilon_T;

% Calculation of the thermal stress in the core
sigma_T_core = E_core*epsilon_T;

% Calculation of the bending moments and curvature
M_x = -sigma_T_skin*t_skin*(2*R-t_skin)/2;
M_y = -sigma_T_skin*t_skin^2/2;
M_xc = -sigma_T_core*t_core*(2*R-t_skin-t_core)/2;
kappa = (M_x+M_xc)/(E_skin*I);

% Calculation of the displacement, rotation, and shear stress
u = -M_y/(E_skin*I)*R;
w = -M_x/(E_skin*I)*R;
phi = -kappa*R;
tau_xy = (E_skin/(2*(1+v_skin)))*(-kappa*t_skin);

% Calculation of the stress in the skins
sigma_xx = -E_skin*kappa*y + sigma_T_skin;
sigma_yy = -E_skin*phi*x;

% Calculation of the stress inthe core
sigma_xx = sigma_xx + E_core*kappa*y;
sigma_yy = sigma_yy + E_core*phi*x;

end

function [u_n, w_n, omega] = curved_panel_vibration(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, N)
% Calculation of the natural frequencies and mode shapes of a curved sandwich panel beam
% Inputs:
% L - Length of the panel (m)
% R - Radius of curvature (m)
% t_skin - Thickness of the skins (m)
% t_core - Thickness of the core (m)
% E_skin - Young's modulus of the skins (Pa)
% E_core - Young's modulus of the core (Pa)
% v_skin - Poisson's ratio of the skins
% v_core - Poisson's ratio of the core
% N - Number of terms in the approximation
% Outputs:
% u_n - Coefficients of the displacement function (m)
% w_n - Coefficients of the deflection function (m)
% omega - Natural frequencies (rad/s)

% Calculation of the moment of inertia and area of the panel
I = (pi/4)*((R+t_skin)^4 - R^4 - (R-t_skin)^4) + t_core*(pi/4)*((R+t_skin+t_core)^4 - (R+t_skin)^4 - (R-t_skin)^4 + (R-t_skin-t_core)^4);
A = (pi/4)*((R+t_skin)^2 - (R-t_skin)^2) + t_core*(pi/4)*((R+t_skin+t_core)^2 - (R+t_skin)^2 - (R-t_skin)^2 + (R-t_skin-t_core)^2);

% Calculation of the stiffness matrix
K = zeros(2*N, 2*N);
for n = 1:N
    K(2*n-1, 2*n-1) = E_skin*A*(n*pi/L)^3 + E_core*t_core*A*(n*pi/L)^3;
    K(2*n-1, 2*n) = -E_skin*I*(n*pi/L)*((n*pi/L)^2*(1-v_skin)+n^2*pi^2*v_skin/L^2) - E_core*t_core*I*(n*pi/L)*((n*pi/L)^2*(1-v_core)+n^2*pi^2*v_core/L^2);
    K(2*n, 2*n-1) = -E_skin*I*(n*pi/L)*((n*pi/L)^2*(1-v_skin)+n^2*pi^2*v_skin/L^2) - E_core*t_core*I*(n*pi/L)*((n*pi/L)^2*(1-v_core)+n^2*pi^2*v_core/L^2);
    K(2*n, 2*n) = E_skin*I*(n*pi/L)^2*((n*pi/L)^2*(1-v_skin)+2*n^2*pi^2*v_skin/L^2) + E_core*t_core*I*(n*pi/L)^2*((n*pi/L)^2*(1-v_core)+2*n^2*pi^2*v_core/L^2);
end

% Calculation of the mass matrix
M = rho*A*L/6*[2*eye(N), zeros(N); zeros(N), 2*eye(N)];

% Calculation of the eigenvalues and eigenvectors
[V, D] = eig(K, M);
omega = sqrt(diag(D));
u_n = V(1:2:end, :);
w_n = V(2:2:end, :);

% Sorting of the eigenvalues and eigenvectors in ascending order
[omega, idx] = sort(omega);
u_n = u_n(:, idx);
w_n = w_n(:, idx);

% Normalization of the eigenvectors
for n = 1:N
    u_n(:, n) = u_n(:, n)/norm(u_n(:, n));
    w_n(:, n) = w_n(:, n)/norm(w_n(:, n));

end
% Calculation of the first mode shape
x = linspace(0, L, 100); % Position along the length of the panel (m)
u = u_n(1)*sin(pi*x/L); % Displacement of the first mode shape (m)
w = w_n(1)*sin(pi*x/L); % Deflection of the first mode shape (m)

% Plotting of the first mode shape
figure;
subplot(2,1,1);
plot(x, u);
title('Displacement of the First Mode Shape');
xlabel('Length (m)');
ylabel('Displacement (m)');
subplot(2,1,2);
plot(x, w);
title('Deflection of the First Mode Shape');
xlabel('Length (m)');
ylabel('Deflection (m)');

% Calculation of the first mode shape
x = linspace(0, L, 100); % Position along the length of the panel (m)
u = u_n(1)*sin(pi*x/L); % Displacement of the first mode shape (m)
w = w_n(1)*sin(pi*x/L); % Deflection of the first mode shape (m)

% Plotting of the first mode shape
figure;
subplot(2,1,1);
plot(x, u);
title('Displacement of the First Mode Shape');
xlabel('Length (m)');
ylabel('Displacement (m)');
subplot(2,1,2);
plot(x, w);
title('Deflection of the First Mode Shape');
xlabel('Length (m)');
ylabel('Deflection (m)');


% Input parameters
L = 1; % Length of the panel (m)
R = 5; % Radius of curvature (m)
t_skin = 0.001; % Thickness of the skins (m)
t_core = 0.01; % Thickness of the core (m)
E_skin = 70e9; % Young's modulus of the skins (Pa)
E_core = 2e9; % Young's modulus of the core (Pa)
v_skin = 0.33; % Poisson's ratio of the skins
v_core = 0.1; % Poisson's ratio of the core
alpha = 23e-6; % Coefficient of thermal expansion (1/K)
T = 100; % Temperature change (K)

% Calculation of the natural frequencies andmode shapes without the core
N = 10; % Number of terms in the approximation
[u_n1, w_n1, omega1] = curved_panel_vibration(L, R, t_skin, 0, E_skin, 0, v_skin, 0, N);

% Calculation of the natural frequencies and mode shapes with the core
[u_n2, w_n2, omega2] = curved_panel_vibration(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, N);

% Comparison of the natural frequencies
fprintf('Natural frequencies without core: %.2f, %.2f, %.2f rad/s\n', omega1(1), omega1(2), omega1(3));
fprintf('Natural frequencies with core: %.2f, %.2f, %.2f rad/s\n', omega2(1), omega2(2), omega2(3));

% Calculation of the thermal stress distribution without the core
[~, ~, ~, ~, sigma_xx1, ~, ~] = thermal_stress(L, R, t_skin, 0, E_skin, 0, v_skin, 0, alpha, T);

% Calculation of the thermal stress distribution with the core
[~, ~, ~, ~, sigma_xx2, ~, ~] = thermal_stress(L, R, t_skin, t_core, E_skin, E_core, v_skin,v_core, alpha, T);

% Comparison of the thermal stress distributions
x = linspace(0, L, 100); % Position along the length of the panel (m)
y = linspace(-R, R, 100); % Position along the width of the panel (m)
[X, Y] = meshgrid(x, y);
figure;
subplot(2,1,1);
contourf(X, Y, sigma_xx1/1e6, 50);
colorbar;
title('Thermal Stress Distribution without Core (MPa)');
xlabel('Length (m)');
ylabel('Width (m)');
subplot(2,1,2);
contourf(X, Y, sigma_xx2/1e6, 50);
colorbar;
title('Thermal Stress Distribution with Core (MPa)');
xlabel('Length (m)');
ylabel('Width (m)');

% Comparison of the first mode shapes
x = linspace(0, L, 100); % Position along the length of the panel (m)
u1 = u_n1(1)*sin(pi*x/L); % Displacement of the first mode shape without the core (m)
w1 = w_n1(1)*sin(pi*x/L); % Deflection of the first mode shape without the core (m)
u2 = u_n2(1)*sin(pi*x/L); % Displacement of the first mode shape with the core (m)
w2 = w_n2(1)*sin(pi*x/L); % Deflection of the first mode shape with the core (m)

figure;
subplot(2,1,1);
plot(x, u1, x, u2);
legend('Without Core', 'With Core');
title('Displacement of the First Mode Shape');
xlabel('Length (m)');
ylabel('Displacement (m)');
subplot(2,1,2);
plot(x, w1, x, w2);
legend('Without Core', 'With Core');
title('Deflection of the First Mode Shape');
xlabel('Length (m)');
ylabel('Deflection (m)');


% Input parameters
L = 1; % Length of the panel (m)
R = 5; % Radius of curvature (m)
t_skin = 0.001; % Thickness of the skins (m)
t_core = 0.01; % Thickness of the core (m)
E_skin = 70e9; % Young's modulus of the skins (Pa)
E_core = 2e9; % Young's modulus of the core (Pa)
v_skin = 0.33; % Poisson's ratio of the skins
v_core = 0.1; % Poisson's ratio of the core
rho_skin = 2700; % Density of the skins (kg/m^3)
rho_core = 1000; % Density of the core (kg/m^3)
N = 10; % Number of terms in the approximation
alpha = 23e-6; % Coefficient of thermal expansion (1/K)
T = 100; % Temperature change (K)

% Calculation of the natural frequencies and mode shapes without the core
[u_n1, w_n1, omega1] = curved_panel_vibration(L, R, t_skin, 0, E_skin, 0, v_skin, 0,N);

% Calculation of the natural frequencies and mode shapes with the core
[u_n2, w_n2, omega2] = curved_panel_vibration(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, N);

% Calculation of the thermal stress distribution without the core
[~, ~, ~, ~, sigma_xx1, ~, ~] = thermal_stress(L, R, t_skin, 0, E_skin, 0, v_skin, 0, alpha, T);

% Calculation of the thermal stress distribution with the core
[~, ~, ~, ~, sigma_xx2, ~, ~] = thermal_stress(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, alpha, T);

% Calculation of the equivalent bending stiffness without the core
k_bend1 = bending_stiffness(t_skin, E_skin, v_skin, L, R);

% Calculation of the equivalent bending stiffness with the core
k_bend2 = bending_stiffness_sandwich(t_skin, t_core, E_skin, E_core, v_skin, v_core, L, R);

% Calculation of the effective stiffness and damping without the core
[keff1, c_eff1] = effective_stiffness_damping(omega1, rho_skin, k_bend1, L, R);

% Calculation of the effective stiffness and damping withthe core
[keff2, c_eff2] = effective_stiffness_damping_sandwich(omega2, rho_skin, rho_core, k_bend2, L, R);

% Comparison of the natural frequencies
fprintf('Natural frequencies without core: %.2f, %.2f, %.2f rad/s\n', omega1(1), omega1(2), omega1(3));
fprintf('Natural frequencies with core: %.2f, %.2f, %.2f rad/s\n', omega2(1), omega2(2), omega2(3));

% Comparison of the effective stiffness and damping
fprintf('Effective stiffness without core: %.2f N/m\n', keff1);
fprintf('Effective stiffness with core: %.2f N/m\n', keff2);
fprintf('Effective damping without core: %.2f Ns/m\n', c_eff1);
fprintf('Effective damping with core: %.2f Ns/m\n', c_eff2);

% Comparison of the thermal stress distributions
x = linspace(0, L, 100); % Position along the length of the panel (m)
y = linspace(-R, R, 100); % Position along the width of the panel (m)
[X, Y] = meshgrid(x, y);
figure;
subplot(2,1,1);
contourf(X, Y, sigma_xx1/1e6, 50);
colorbar;
title('Thermal Stress Distribution without Core (MPa)');
xlabel('Length (m)');
ylabel('Width (m)');
subplot(2,1,2);
contourf(X, Y, sigma_xx2/1e6, 50);
colorbar;
title('Thermal Stress Distribution with Core (MPa)');
xlabel('Length (m)');
ylabel('Width (m)');

% Comparison of the first mode shapes
x = linspace(0, L, 100); % Position along the length of the panel (m)
u1 = u_n1(1)*sin(pi*x/L); % Displacement of the first mode shape without the core (m)
w1 = w_n1(1)*sin(pi*x/L); % Deflection of the first mode shape without the core (m)
u2 = u_n2(1)*sin(pi*x/L); % Displacement of the first mode shape with the core (m)
w2 = w_n2(1)*sin(pi*x/L); % Deflection of the first mode shape with the core (m)

figure;
subplot(2,1,1);
plot(x, u1, x, u2);
legend('Without Core', 'With Core');
title('Displacement of the First Mode Shape');
xlabel('Length (m)');
ylabel('Displacement (m)');
subplot(2,1,2);
plot(x, w1, x, w2);
legend('Without Core', 'With Core');
title('Deflection of the First Mode Shape');
xlabel('Length (m)');
ylabel('Deflection (m)');

% Input parameters
L = 1; % Length of the panel (m)
R = 5; % Radius of curvature (m)
t_skin = 0.001; % Thickness of the skins(m)
t_core = 0.01; % Thickness of the core (m)
E_skin = 70e9; % Young's modulus of the skins (Pa)
E_core = 2e9; % Young's modulus of the core (Pa)
v_skin = 0.33; % Poisson's ratio of the skins
v_core = 0.1; % Poisson's ratio of the core
rho_skin = 2700; % Density of the skins (kg/m^3)
rho_core = 1000; % Density of the core (kg/m^3)
N = 10; % Number of terms in the approximation
alpha = 23e-6; % Coefficient of thermal expansion (1/K)
T = 100; % Temperature change (K)

% Calculation of the natural frequencies and mode shapes without the core
[u_n1, w_n1, omega1] = curved_panel_vibration(L, R, t_skin, 0, E_skin, 0, v_skin, 0,N);

% Calculation of the natural frequencies and mode shapes with the core
[u_n2, w_n2, omega2] = curved_panel_vibration(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, N);

% Calculation of the thermal stress distribution without the core

% Comparison of the natural frequencies
fprintf('Natural frequencies without core: %.2f, %.2f, %.2f rad/s\n', omega1(1), omega1(2), omega1(3));
fprintf('Natural frequencies with core: %.2f, %.2f, %.2f rad/s\n', omega2(1), omega2(2), omega2(3));

% Calculation of the equivalent bending stiffness without the core
k_bend1 = bending_stiffness(t_skin, E_skin, v_skin, L, R);

% Calculation of the equivalent bending stiffness with the core
k_bend2 = bending_stiffness_sandwich(t_skin, t_core, E_skin, E_core, v_skin, v_core, L, R);

% Calculation of the effective stiffness and damping without the core
[keff1, c_eff1] = effective_stiffness_damping(omega1, rho_skin, k_bend1, L, R);

% Calculation of the effective stiffness and damping with the core
[keff2, c_eff2] = effective_stiffness_damping_sandwich(omega2, rho_skin, rho_core, k_bend2, L, R);

% Comparison of the effective stiffness and damping
fprintf('Effective stiffness without core: %.2f N/m\n', keff1);


% Calculation of the thermal stress distribution with the core
[~, ~, ~, ~, sigma_xx2, ~, ~] = thermal_stress(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, alpha, T);

% Calculation of the displacement and deflection of the panel without the core
x = linspace(0, L, 100); % Position along the length of the panel (m)
u1 = u_n1(1)*sin(pi*x/L); % Displacement of the first mode shape without the core (m)
w1 = w_n1(1)*sin(pi*x/L); % Deflection of the first mode shape without the core (m)
u2 = u_n2(1)*sin(pi*x/L); % Displacement of the first mode shape with the core (m)
w2 = w_n2(1)*sin(pi*x/L); % Deflection of the first mode shape with the core (m)

% Plotting the natural frequencies and mode shapes
figure;
subplot(2,2,1);
bar(omega1(1:3));
title('Natural Frequencies without Core (rad/s)');
xlabel('Mode Number');
ylabel('Frequency');
subplot(2,2,2);
bar(omega2(1:3));
title('Natural Frequencies with Core (rad/s)');
xlabel('Mode Number');
ylabel('Frequency');
subplot(2,2,3);
x = linspace(0, L, 100); % Position along the length of the panel (m)
plot(x, u1, x, w1);
legend('Displacement', 'Deflection');
title('First Mode Shape without Core');
xlabel('Length (m)');
ylabel('Displacement or Deflection (m)');
subplot(2,2,4);
x = linspace(0, L, 100); % Position along the length of the panel (m)
plot(x, u2, x, w2);
legend('Displacement', 'Deflection');
title('First Mode Shape with Core');
xlabel('Length (m)');
ylabel('Displacement or Deflection (m)');

% Plotting the thermal stress distribution
figure;
subplot(2,1,1);
contourf(X, Y, sigma_xx1/1e6, 50);
colorbar;
title('Thermal Stress Distribution without Core (MPa)');
xlabel('Length (m)');
ylabel('Width (m)');
subplot(2,1,2);
contourf(X, Y, sigma_xx2/1e6, 50);
colorbar;
title('Thermal Stress Distribution with Core (MPa)');
xlabel('Length (m)');
ylabel('Width (m)');

% Calculation of the effective stiffness and damping without the core


% Comparison of the effective stiffness and damping
fprintf('Effective stiffness without core: %.2f N/m\n', keff1);
fprintf('Effective stiffness with core: %.2f N/m\n', keff2);
fprintf('Effective damping without core: %.2f Ns/m\n', c_eff1);
fprintf('Effective damping with core: %.2f Ns/m\n', c_eff2);

% Calculation of the displacement and deflection of the panel without the core
x = linspace(0, L, 100); % Position along the length of the panel (m)
u1 = u_n1(1)*sin(pi*x/L); % Displacement of the first mode shape without the core (m)
w1 = w_n1(1)*sin(pi*x/L); % Deflection of the first mode shape without the core (m)
u2 = u_n2(1)*sin(pi*x/L); % Displacement of the first mode shape with the core (m)
w2 = w_n2(1)*sin(pi*x/L); % Deflection of the first mode shape with the core (m)

% Calculation of the thermal stress distribution without the core
[~, ~, ~, ~, sigma_xx1, ~, ~] = thermal_stress(L, R, t_skin, 0, E_skin, 0, v_skin, 0, alpha, T);

% Calculation of the thermalstress distribution with the core
[~, ~, ~, ~, sigma_xx2, ~, ~] = thermal_stress(L, R, t_skin, t_core, E_skin, E_core, v_skin, v_core, alpha, T);

% Calculation of the displacement and deflection of the panel without the core
x = linspace(0, L, 100); % Position along the length of the panel (m)
u1 = u_n1(1)*sin(pi*x/L); % Displacement of the first mode shape without the core (m)
w1 = w_n1(1)*sin(pi*x/L); % Deflection of the first mode shape without the core (m)
u2 = u_n2(1)*sin(pi*x/L); % Displacement of the first mode shape with the core (m)
w2 = w_n2(1)*sin(pi*x/L); % Deflection of the first mode shape with the core (m)

% Calculation of the thermal stress distribution without the core
[~, ~, ~, ~, sigma_xx1, ~, ~] = thermal_stress(L, R, t_skin, 0, E_skin, 0, v_skin, 0, alpha, T);

% Calculation of the thermal stress distribution with the core
[~, ~, ~, ~, sigma_xx2, ~, ~] = thermal_stress(L, R, t_skin, t_core, skin, E_core, v_skin, v_core, alpha, T);

% Calculation of the displacement and deflection of the panel without the core
x = linspace(0, L, 100); % Position along the length of the panel (m)
u1 = u_n1(1)*sin(pi*x/L); % Displacement of the first mode shape without the core (m)
w1 = w_n1(1)*sin(pi*x/L); % Deflection of the first mode shape without the core (m)
u2 = u_n2(1)*sin(pi*x/L); % Displacement of the first mode shape with the core (m)