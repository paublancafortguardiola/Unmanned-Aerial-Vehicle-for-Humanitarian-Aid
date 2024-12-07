%% FORWARD FLIGHT PERFORMANCE ANALYSIS

%% Problem Data
%General
m = 6100; % kg Total mass
f = 1.85; % m^2 
l_t = 8; % m Main-tail rotor distance

%Engine
P_eng = 2558e3; % W Total Power

% Main Rotor
Rm = 6.5; % m Radius
Nbm = 5; % number of blades
rpm = 270;
cm = 0.6; % m Chord
eta_m = 0.97; % Efficiency
k_m = 1.15; % Non-idealities factor
Cdo_m = 0.009; % Parasitic drag coefficient

% Tail Rotor
Rt = 1.3; % m Radius
Nbt = 10; % number of blades
rpm_t = 2000;
ct = 0.2; % m Chord
eta_t = 0.97; % Efficiency
k_t = 1.15; % Non-idealities factor
Cdo_t = 0.00618; % Parasitic drag coefficient

rho = 1.225; % kg/m^3 Density

%Tail/Main Rotor Power ratio
ratio = 0.15;
P_eng_mr_av = P_eng*(1-ratio);
P_sh_mr_av = P_eng_mr_av*eta_m;

v_inf=2:2:304; % Speeds (km/h)
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

% Curvas de potencia
plot(v_inf, Pc, 'b-', 'LineWidth', 2, 'DisplayName', 'Climbing Power'); % Power for climb
plot(v_inf, Pi, 'g-', 'LineWidth', 2, 'DisplayName', 'Induced Power'); % Induced power
plot(v_inf, Po, 'm-', 'LineWidth', 2, 'DisplayName', 'Profile Power'); % Profile power
plot(v_inf, Pf, 'c-', 'LineWidth', 2, 'DisplayName', 'Forward Power'); % Parasite power

% Potencia del eje
plot(v_inf, P_sh, 'k-', 'LineWidth', 2, 'DisplayName', 'Main Rotor Shaft Power');

% Potencia del eje disponible como línea discontinua
plot(v_inf, P_sh_av, 'r--', 'LineWidth', 2, 'DisplayName', 'Shaft Power Available');

% Personalización de la gráfica
xlabel('Forward Flight Speed, v_{inf} (km/h)', 'FontSize', 12);
ylabel('Power (W)', 'FontSize', 12);
title('Power vs Forward Flight Speed at sea-level', 'FontSize', 14);

% Límites de los ejes
xlim([0, 310]); % Limitar el eje x entre 0 y 310

% Leyenda más grande
lgd = legend('show', 'Location', 'best'); % Crear la leyenda
lgd.FontSize = 14; % Aumentar el tamaño de la fuente de la leyenda

grid on;
hold off;

%% Excess Power Graph
% Calculate excess power
Excess_Power = P_sh_av - P_sh;

% Plot the excess power
figure;
hold on;

% Excess power curve
plot(v_inf, Excess_Power, 'b-', 'LineWidth', 2, 'DisplayName', 'Excess Power');

% Personalización de la gráfica
xlabel('Forward Flight Speed, v_{inf} (km/h)', 'FontSize', 12);
ylabel('Excess Power (W)', 'FontSize', 12);
title('Excess Power vs Forward Flight Speed', 'FontSize', 14);

% Límites de los ejes
xlim([0, 310]); % Limitar el eje x entre 0 y 310

grid on;
hold off;

%% Range and Endurance Plot

g = 9.81; % m/s^2
SFC = 0.546 / 1000; % kg/Wh
eta_p = 0.9;
Efficiency = 5;
W_0 = 8000; % kg (MTOW)
W_end = 6000; % kg (MZFW)
v_inf = (5:5:300) / 3.6; % m/s 

% Element-wise division for v_inf and scalar log term
Endurance = eta_p ./ (g * SFC * v_inf) * Efficiency * log(W_0 / W_end);

% Plot Endurance vs Velocity
figure;
plot(v_inf * 3.6, Endurance, 'LineWidth', 1.5);

% Add labels and title
xlabel('Velocity (km/h)', 'FontSize', 12);
ylabel('Endurance (hours)', 'FontSize', 12);
title('Endurance vs Velocity', 'FontSize', 14);

% Add grid
grid on;

% Add axis limits and ticks
xlim([0, 300]);
ylim([0, max(Endurance) + 0.5]);

% Add annotations
set(gca, 'FontSize', 12, 'GridLineStyle', '--'); % Customize axis and grid



