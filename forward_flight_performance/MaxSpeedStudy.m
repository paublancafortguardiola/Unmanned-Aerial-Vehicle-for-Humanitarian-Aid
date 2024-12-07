%% FORWARD FLIGHT PERFORMANCE ANALYSIS

%% Problem Data
%General
m = 8000; % kg Total mass
f = 1.85; % m^2 
l_t = 8; % m Main-tail rotor distance

% Main Rotor (NACA 23012) 
Rm = 6.5; % m Radius
Nbm = 5; % number of blades
rpm = 270;
cm = 0.6; % m Chord
eta_m = 0.97; % Efficiency
k_m = 1.15; % Non-idealities factor
Cdo_m = 0.00991; % Parasitic drag coefficient

rho = 0.5; % kg/m^3 Density
n = 0.6;
P_eng = 2558e3*(rho/1.225)^n; % W Total Power

% Tail Rotor
Rt = 1.3; % m Radius
Nbt = 10; % number of blades
rpm_t = 2000;
ct = 0.2; % m Chord
eta_t = 0.97; % Efficiency
k_t = 1.15; % Non-idealities factor
Cdo_t = 0.00618  ; % Parasitic drag coefficient

%Tail/Main Rotor Power ratio
ratio = 0.15;
P_eng_mr_av = P_eng*(1-ratio);
P_sh_mr_av = P_eng_mr_av*eta_m;

v_inf=5:10:350; % Speeds (km/h)
Pc = zeros(length(v_inf),1);
Pi = zeros(length(v_inf),1);
Po = zeros(length(v_inf),1);
Pf = zeros(length(v_inf),1);
P_sh = zeros(length(v_inf),1);
P_sh_av = zeros(length(v_inf),1);
v_inf_ms = zeros(length(v_inf),1);


for i=1:length(v_inf)
    v_inf_ms(i)=v_inf(i)/3.6;
    % Call the function
    [T_sol(i), alpha_sol, lambda_i_sol, Pc(i), Pi(i), Po(i), Pf(i), P_sh(i), P_sh_av(i)] = CalculationsForward(l_t,Rt, Nbt, rpm_t, ct, eta_t, k_t, Cdo_t, Cdo_m, rho, v_inf_ms(i), m, f, P_eng, Rm, Nbm, rpm, cm, eta_m, k_m);
end


%% Power Analysis Graph
figure;
hold on;

% Potencia del eje
plot(v_inf, P_sh, 'k-', 'LineWidth', 2, 'DisplayName', 'Main Rotor Shaft Power');

% Potencia del eje disponible como línea discontinua
plot(v_inf, P_sh_av, 'r--', 'LineWidth', 2, 'DisplayName', 'Shaft Power Available');

% Personalización de la gráfica
xlabel('Forward Flight Speed, v_{inf} (km/h)', 'FontSize', 12);
ylabel('Power (W)', 'FontSize', 12);
title('Power vs Forward Flight Speed at sea-level', 'FontSize', 14);

% Límites de los ejes
xlim([0, 350]); % Limitar el eje x entre 0 y 310

% Leyenda más grande
lgd = legend('show', 'Location', 'best'); % Crear la leyenda
lgd.FontSize = 14; % Aumentar el tamaño de la fuente de la leyenda

grid on;
hold off;

%% Max speeds

% Data
heights = [0, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000];
v_maxs = [310, 308, 307, 305, 303.5, 296, 287, 275, 245];
P_avs = [1, 0.953, 0.908, 0.864, 0.822, 0.742, 0.669, 0.601, 0.559] * P_eng; % Power in W

% Create the first plot for speed
figure;
yyaxis left; % Left y-axis
plot(heights, v_maxs, '-', 'LineWidth', 2, 'Color', 'b');
ylabel('Maximum Speed, v_{max} (km/h)', 'FontSize', 14);
ylim([240, 310]); % Set limits for speed axis

% Create the second plot for power
yyaxis right; % Right y-axis
plot(heights, P_avs / 1e6, '-', 'LineWidth', 2, 'Color', 'r'); % Convert power to MW for clarity
ylabel('Available Power, P_{av} (MW)', 'FontSize', 14);
ylim([0, P_eng / 1e6]); % Set limits for power axis

% Common x-axis
xlabel('Altitude (m)', 'FontSize', 14);

% Title and grid
title('Maximum Speed and Available Power vs Altitude', 'FontSize', 20);
grid on;

% Style enhancements
set(gca, 'FontSize', 12); % Adjust tick label font