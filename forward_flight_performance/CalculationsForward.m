function [T_sol, alpha_sol, lambda_i_sol, Pc, Pi, Po, Pf,P_sh_mr, P_sh_mr_av] = CalculationsForward(l_t,Rt, Nbt, rpm_t, ct, eta_t, k_t, Cdo_t, Cdo_m, rho, v_inf, m, f, P_eng, Rm, Nbm, rpm, cm, eta_m, k_m)
    omega = 2 * pi / 60 * rpm;
    g = 9.81;
    A_m=pi*Rm^2;
    D = 1/2 * rho * v_inf^2 * f; % Drag
    V_tip=omega*Rm;
    solidity_m=Nbm*cm*Rm/A_m;
    
    % Define symbolic variables
    syms T alpha_sym lambda_i_sym

    % Forces equilibrium
    eq1 = T * cos(alpha_sym) == m * g;
    eq2 = T * sin(alpha_sym) == D;

    % Solve the equations
    solutions = solve([eq1, eq2], [T, alpha_sym]);

    % Extract numeric solutions
    T_values = double(solutions.T);
    alpha_values = double(solutions.alpha_sym);

    % Filter positive solutions
    positive_indices = T_values > 0 & alpha_values > 0;
    T_sol = T_values(positive_indices);
    alpha_sol = alpha_values(positive_indices);
    
    if v_inf==265/3.6
        alpha_sol*180/pi
    end

    % Ensure only one solution is returned
    if numel(T_sol) ~= 1 || numel(alpha_sol) ~= 1
        error('Multiple or no positive solutions found.');
    end

    Vih = sqrt(T_sol/(2*rho*A_m)); %Induced velocity
    lambda_ih=Vih/V_tip;
    Poh = rho*A_m*V_tip^3*1/8*solidity_m*Cdo_m; % Profile Power at hover
    
    %Forward flight calculations
    Vz = v_inf*sin(alpha_sol);
    mu = v_inf*cos(alpha_sol)/V_tip;

    eq3 = lambda_i_sym == lambda_ih^2 / sqrt(mu^2 + (mu * tan(alpha_sol) + lambda_i_sym)^2);

    % Solve numerically with a reasonable initial guess
    solution = vpasolve(eq3, lambda_i_sym, [0, inf]); % Assume lambda_i_sym > 0
    lambda_i_sol = double(solution);

    % Ensure a valid solution is found
    if isempty(lambda_i_sol)
        error('No valid solution found for lambda_i.');
    end

    %Power calculations
    Pc = T_sol*Vz; % Climbing Power
    Pi = k_m*T_sol*lambda_i_sol*V_tip; % Induced Power
    Po = Poh*(1+4.65*mu^2); % Profile Power
    Pf = 1/2*rho*V_tip^3*mu^3*f; % Forward Power
    
    P_sh_mr = Pc + Pi + Po + Pf; 
    P_mr = P_sh_mr/eta_m;

    %Tail Rotor
    M_mr = P_sh_mr/omega;
    T_t = M_mr/l_t; % Tail rotor thrust
    A_t = pi*Rt^2; % Tail rotor area
    Vih_t = sqrt(T_t/(2*rho*A_t)); % Induced speed
    Pih_t = k_t*T_t*Vih_t; % Induced power at hover
    solidity_t = Nbt*ct*Rt/A_t;
    omega_t = 2*pi/60*rpm_t; % Tail rotor rotational speed
    V_tip_t = omega_t*Rt;
    Poh_t = rho*A_t*V_tip_t^3*1/8*solidity_t*Cdo_t; % Profile Power at hover
    P_sh_t = Pih_t + Poh_t; % Shaft power tail rotor
    P_t = P_sh_t/eta_t;

    ratio = P_t/P_mr;
    P_eng_mr = P_eng*(1-ratio);
    P_sh_mr_av = P_eng_mr*eta_m;

end
