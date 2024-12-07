function [r,theta0,theta_vec,alpha_vec,Cl_vec,dCT_vec,dCP_vec,CP, P_req] = bemt_solver_vectorized(R, y0, Nb, c, theta_tw, vz, Cl_alpha, Cd0, mass, rho, rpm, g, n, max_F_error, max_CT_error, twist_type)
    % BEMT Calculator
    %% Previous calculations
    
    % Nodes calculation
    r0 = y0/R;
    [r, dr] = getNodes(r0,1,n);
    
    % Parameters calculation
    A = pi*R^2;
    sigma = Nb*c*R/A; % solidity [-];
    T = mass*g; % Thrust [N]
    omega = rpm*2*pi/60; % Angular velocity [rad/s]
    v_tip = omega*R; % Velocity of the tip of the rotor blade [m/s]
    lambda_c = vz/v_tip; % Axial climb inflow ratio
    CT_req = T/(rho*A*v_tip^2); % Required thrust coefficient [-]
    lambda_ih_req = sqrt(CT_req/2); % Induced hover inflow ratio using the required thrust coefficient
    lambda_c_hat_req = lambda_c/lambda_ih_req; % Axial climb inflow ratio wrt required induced inflow ratio on hover
    if ( (lambda_c_hat_req>-2) && (lambda_c_hat_req<0) )
        warning('lambda_c_hat = '+string(lambda_c_hat_req)+' and should not be between -2 and 0, try with another climb velocity or different RPMs')
    end
    lambda_i_hat_func = @(x) -x/2 + (1*(x>=0) - 1*(x<=-2))*sqrt( (x^2)/4 + 1*(x>=0) -1*(x<=-2) ); % Function of lambda_c_hat
    lambda_i_hat_req = lambda_i_hat_func(lambda_c_hat_req); % Iduced inflow ratio hat, using lamba_c_hat_req
    
    
    %% Solution algorithm
    % Variables initialization:
    lambda_vec = zeros(1,n); % Vector of inflow ratios at each blade element
    F_vec = zeros(1,n); % Vector of Prandtl's correction factor at each blade element
   
    if twist_type == "linear"
        % Starting theta0 estimation
        theta0_est = 6*CT_req/(sigma*Cl_alpha) - (3/4)*theta_tw + (3/2)*( lambda_c + sqrt(CT_req/2)*lambda_i_hat_req ); % Not well defined for when vz < 0, why?? (theta0_est should be positive for convergence when vz<0, otherwise lambda < 0 and f < 0 and F doesn't exist)
        % Function to calculate the blade angle at each blade element
        % (given a linear twist)
        get_theta_vec = @(r_, theta0_est_) theta0_est_ + r_.*theta_tw; 
    elseif twist_type == "hyperbolic"
        % Starting theta0 estimation
        theta0_est = 6*CT_req/(sigma*Cl_alpha) - (3/2)*theta_tw + (3/2)*( lambda_c + sqrt(CT_req/2)*lambda_i_hat_req ); % Not well defined for when vz < 0, why?? (theta0_est should be positive for convergence when vz<0, otherwise lambda < 0 and f < 0 and F doesn't exist)
        % Function to calculate the blade angle at each blade element
        % (given an hyperbolic twist)
        get_theta_vec = @(r_, theta0_est_) theta0_est_ + theta_tw./r_;
    end

    %%
    % Iterative algorithm until thrust coefficient convergence to the required
    CT_converged = false;
    iterations=0;
   
    while ~CT_converged
    iterations = iterations + 1;
    theta_vec = get_theta_vec(r,theta0_est);

    % For each blade element (vectorized)
    % Iterative algorithm to find the inflow ratio (lambda) and Prandtl's correction factor (F)
    F_est = ones(1,n); % Starting Prandtl's correction factor estimation
    F_error = max_F_error+1;
    % Convergence criteria on Prandtl's correction factor
    while any(F_error > max_F_error)
        lambda_vec = -( sigma*Cl_alpha./(16*F_est) - lambda_c/2 ) + sqrt( (sigma*Cl_alpha./(16*F_est) - lambda_c/2).^2 + sigma*Cl_alpha./(8*F_est).*theta_vec.*r );
        f = (Nb/2)*((1-r)./lambda_vec);
        F_vec = (2/pi)*acos(exp(-f));
        F_error = abs(F_vec - F_est);
        F_est = F_vec;
    end

    dCT_vec = 4*F_vec.*(lambda_vec - lambda_c).*lambda_vec.*r;
    CT = sum(dCT_vec.*dr,2);
    CT_error = abs(CT - CT_req);

    % Convergence criteria on thrust coefficient
    if CT_error > max_CT_error
        theta0 = theta0_est + 6*(CT_req -  CT)/(sigma*Cl_alpha) + (3*sqrt(2)/4)*(sqrt(CT_req) - sqrt(CT));
        theta0_est = theta0;
    else
        CT_converged = true;
    end
    end
    
    %% Final calculations
    Cl_vec = Cl_alpha*(theta_vec - lambda_vec./r); % Local lift coefficient at each blade element
    Cl = sum(Cl_vec.*dr); % Total blade lift coefficient
    lambda_i_vec = lambda_vec - lambda_c; % Induced inflow ratio at each blade element
    dCPi_vec = lambda_i_vec.*dCT_vec; % Induced power coefficient at each balde element
    dCPc_vec = lambda_c*dCT_vec; % Axial climb power coefficient at each blade element
    dCP0_vec = (sigma/2)*Cd0*r.^3; % Profile power coefficient at each blade element
    dCP_vec = dCPi_vec + dCPc_vec + dCP0_vec; % Total power coefficient at each blade element
    CP = sum(dCP_vec.*dr); % Power coefficient
    P_req = CP*(rho*A*v_tip^3); % Required power for the given flight condition

    phi_vec = lambda_vec./r; % Inflow angle at each blade element [rad]
    alpha_vec = theta_vec - phi_vec; % Angle of attack at each blade element [rad]
end