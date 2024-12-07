%%
close all; clc; clear;
%% Input data

% Physical
R = 6.5; % Length of the balde (radius) [m]
y0 = 0.5; % Radial position of the root of the blade [m]
Nb = 5; % Number of blades
c = 0.6; % Chord size [m]
twist_type = "hyperbolic";
theta_tw = deg2rad(-20); % Linear twist coefficient [rad] (theta(r) = theta0 + r*theta_tw) (angle increase between tip and root of the blade)
vz = 46/60; % Axial velocity (for later calculation of lambda_c (inflow ratio)) [m/s]

Cl_alpha = 2*pi; % Lift slope [1/rad]
Cd0 = 0.009; % Profile's parasitic drag coefficient (zero lift drag coefficient)
mass = 8000; % Required mass to hover [kg]
rho = 1.225; % Air density [kg/m^3] (1.18955 @ 1000ft)
rpm = 270; % Rotations per minute of the rotor
g = 9.81; % Gravity [m/s^2] (9.803572201306 @ 1000ft)

Cl_max = 1.55; % Maximum lift coefficient of NACA 23012
alpha_max_deg = 15; % Maximum angle of attack of NACA 23012 [deg]
P_available = 2*1279*1e3;

% Numerical
n = 500; % number of elements
max_F_error = 1e-12; % Maximum error until convergence of the Prandl's correction factor
max_CT_error = 1e-12; % Maximum error until convergence of the thrust coefficient

%% Plots
set(0, 'DefaultTextInterpreter', 'latex');

colorPalette = [
    239 71 111;    % RGB for ef476f
    255 209 102;   % RGB for ffd166
    6 214 160;     % RGB for 06d6a0
    17 138 178;    % RGB for 118ab2
    7 59 76  ;     % RGB for 073b4c
    83 46 99;
    239 71 111;    % RGB for ef476f
    255 209 102;   % RGB for ffd166
    6 214 160;     % RGB for 06d6a0
    17 138 178;    % RGB for 118ab2
    7 59 76        % RGB for 073b4c
    83 46 99;
] / 255; % Normalize RGB values to [0, 1] range

legend_entries = {};
P_req_list = {};

for vz=0:5:20
    [r,theta0,theta_vec,alpha_vec,Cl_vec,dCT_vec,dCP_vec,CP,P_req] = bemt_solver_vectorized(R, y0, Nb, c, theta_tw, vz, Cl_alpha, Cd0, mass, rho, rpm, g, n, max_F_error, max_CT_error, twist_type);
    alpha_vec_deg = rad2deg(alpha_vec); % Angle of attack at each blade element [deg]
    P_req_list{end+1} = P_req;
    
    figure(1);
    colororder(colorPalette);
    plot(r,Cl_vec,'LineWidth',1.3);
    hold on;

    figure(2);
    colororder(colorPalette);
    plot(r,dCP_vec,'LineWidth',1.3)
    hold on;

    figure(3);
    colororder(colorPalette);
    plot(r,dCT_vec,'LineWidth',1.3)
    hold on;
    
    %legend_entry = '$v_z = 1000$ ft/min';
    legend_entry = sprintf('$v_z = %.0f$ m/s', vz);
    legend_entries{end+1} = legend_entry;
end

%CP_available = P_available/(rho*pi*R^2*(rpm*2*pi/60*R)^3);
P_excess = P_available - cell2mat(P_req_list)
figure(1);
yline(Cl_max,'LineWidth',1,'LineStyle','--','Color','black');
legend_entries{end+1} = '$C_{l_{max}}$';
title('Distribution of lift coefficient along the blade')
xlabel('$r$ $[]$');
ylabel('$C_l$ $[]$');
grid minor, grid on;
legend(legend_entries,'Interpreter','latex', 'Location','best')

figure(2);
title('Distribution of the required power coefficient')
xlabel('$r$ $[]$');
ylabel('$C_P$ $[]$');
grid minor, grid on;
legend(legend_entries{1:end-1},'Interpreter','latex', 'Location','best')

figure(3);
title('Distribution of the required thrust coefficient')
xlabel('$r$ $[]$');
ylabel('$C_T$ $[]$');
grid minor, grid on;
legend(legend_entries{1:end-1},'Interpreter','latex', 'Location','best')

%%
P_excess = P_available - cell2mat(P_req_list)