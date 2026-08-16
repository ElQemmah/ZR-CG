classdef RobustZRCommandGovernor < handle
    %% ZONOTOPIC ROBUST COMMAND GOVERNOR (ONLINE QP)
    % Version for Simulation Example 1

    properties
        M_c           % Center of the interval hull B_cl
        Delta         % The FULL radius matrix

        c_d           % Center of disturbance zonotope D
        Q_d           % Generators of disturbance zonotope D

        v_max         % The velocity constraint limit
        k0            % Prediction horizon
        Psi           % Reference weight matrix
        solver_name   % Name of the QP solver

        dim_x         % State dimension
        dim_g         % Command dimension

        T_j           % Constraint row vector [0 1 0]
        Tj_Mc         % T_j * M_c

        % Tj_Delta_abs = |T_j| * Delta
        Tj_Delta_abs  % 1x(n+p) row vector

        % Precomputed constant radius from disturbance
        const_radius_d

        % --- Useful for Terminal constraint ---
        delta_inf     % Infinite-horizon disturbance tail
        delta_k0      % Residual transient beyond k0
        Tj_Hc_Gam     % T_j*H^c*(I-Phi_c)^-1        (1 x n)
        cg_slope      % T_j*H^c*(I-Phi_c)^-1*G_c + T_j*L^c : equilibrium slope
        Gam_Gc        % (I-Phi_c)^-1*G_c            (n x p)
        absGam        % |(I-Phi_c)^-1|              (n x n)
        Wbar          % (I - |Gam_c|*Delta_Phi)^-1  (n x n, entrywise >= 0)
        Delta_Phi     % radius block on Phi         (n x n)
        Delta_G       % radius block on G           (n x p)
        absTjHc       % |T_j*H^c|                   (1 x n)
        use_terminal  % Enforce (29) when margins are provided
    end

    methods
        function obj = RobustZRCommandGovernor(system_interval_hull, disturbance_zonotope, vel_limit, horizon, weight_matrix, solver, delta_inf, delta_k0)

            if ~isa(system_interval_hull, 'intervalMatrix')
                error('Input must be an intervalMatrix (the B_cl set).');
            end

            % Get bounds and compute M_c, Delta
            intermediate_int_obj = system_interval_hull.int;
            L = intermediate_int_obj.inf;
            U = intermediate_int_obj.sup;
            obj.M_c = (L + U) / 2;
            obj.Delta = (U - L) / 2; % The full, coupled Delta

            obj.dim_x = size(obj.M_c, 1);
            obj.dim_g = size(obj.M_c, 2) - obj.dim_x;

            % Store disturbance zonotope properties
            obj.c_d = disturbance_zonotope.center;
            obj.Q_d = disturbance_zonotope.Z(:, 2:end); % Generators only

            obj.v_max = vel_limit;
            obj.k0 = horizon;
            obj.Psi = weight_matrix;
            obj.solver_name = solver;

            % --- Precompute FULLY-COUPLED terms ---
            obj.T_j = zeros(1, obj.dim_x);
            obj.T_j(2) = 1; % T_j = [0 1 0]
            obj.Tj_Mc = obj.T_j * obj.M_c; % 1x(n+p) vector

            % This is |T_j| * Delta from the proof
            obj.Tj_Delta_abs = abs(obj.T_j) * obj.Delta;

            % Precompute constant radius from disturbance
            obj.const_radius_d = sum(abs(obj.T_j * obj.Q_d));

            % Terminal constraint
            % Enforced only if the offline margins are supplied.
            obj.use_terminal = (nargin >= 8) && ~isempty(delta_inf) && ~isempty(delta_k0);
            if obj.use_terminal
                obj.delta_inf = delta_inf;
                obj.delta_k0  = delta_k0;

                n_  = obj.dim_x;
                Phi_c = obj.M_c(:, 1:n_);      G_c = obj.M_c(:, n_+1:end);
                obj.Delta_Phi = obj.Delta(:, 1:n_);
                obj.Delta_G   = obj.Delta(:, n_+1:end);

                Gam = (eye(n_) - Phi_c) \ eye(n_);       % Gamma_c
                obj.Gam_Gc = Gam * G_c;                  % nominal equilibrium map
                obj.absGam = abs(Gam);

                % H^c = I and L^c = 0 in this example
                obj.absTjHc   = abs(obj.T_j);
                obj.Tj_Hc_Gam = obj.T_j * Gam;
                obj.cg_slope  = obj.T_j * obj.Gam_Gc;    % T_j*c_g^(M_c) per unit g

                % (I - |Gamma_c|*Delta_Phi)^-1 : well posed iff rho(.) < 1
                Sm = obj.absGam * obj.Delta_Phi;
                if max(abs(eig(Sm))) >= 1
                    error(['rho(|Gamma_c|*Delta_Phi) >= 1: the bound (46) on ', ...
                           'delta_sigma is not well defined.']);
                end
                obj.Wbar = (eye(n_) - Sm) \ eye(n_);
            end
        end

        function [g, status] = compute_cmd(obj, x_t, r_t)

            yalmip('clear');
            g = sdpvar(obj.dim_g, 1);
            constraints = [];

            % k=0
            mu_k = x_t;

            % L_k_state: a cell array of affine expressions for
            % the generators of the full state
            L_k_state = {}; 

            for k = 1:obj.k0

                % 1. Form input center c_z 
                c_z = [mu_k; g];

                % 2. Propagate FULL STATE Center
                mu_k_next = obj.M_c * c_z + obj.c_d;

                % Get the affine center of the velocity
                nu_k = mu_k_next(2);

                % 3. Propagate Radius Terms 

                if isempty(L_k_state)
                    rad_Rz_state = zeros(obj.dim_x, 1);
                else
                    rad_Rz_state = sdpvar(obj.dim_x, 1);
                    L_km1_matrix = [L_k_state{:}];
                    abs_L_km1 = sdpvar(size(L_km1_matrix,1), size(L_km1_matrix,2));
                    constraints = [constraints, ...
                        abs_L_km1 >= L_km1_matrix, ...
                        abs_L_km1 >= -L_km1_matrix, ...
                        rad_Rz_state == sum(abs_L_km1, 2)];
                end

                % Build the new list of radius terms L_k_vel
                new_L_k_vel = {};

                % B1 terms: T_j * M_c * Q_z
                for i = 1:length(L_k_state)
                    Q_z_i = [L_k_state{i}; zeros(obj.dim_g, 1)];
                    new_L_k_vel{end+1} = obj.Tj_Mc * Q_z_i;
                end

                % B2 terms: T_j * Diag(Delta * |c_z|)
                abs_c_z = sdpvar(obj.dim_x + obj.dim_g, 1);
                constraints = [constraints, abs_c_z >= c_z, abs_c_z >= -c_z];
                for i = 1:(obj.dim_x + obj.dim_g)
                    new_L_k_vel{end+1} = obj.Tj_Delta_abs(i) * abs_c_z(i);
                end

                % B3 terms: T_j * Diag(Delta * rad(R_z))
                rad_Rz_input = [rad_Rz_state; zeros(obj.dim_g, 1)];
                for i = 1:(obj.dim_x + obj.dim_g)
                    new_L_k_vel{end+1} = obj.Tj_Delta_abs(i) * rad_Rz_input(i);
                end

                % 4. Build the Constraint for this step k
                if isempty(new_L_k_vel)
                    rad_k_dynamic = 0;
                else
                    rad_k_dynamic = sdpvar(1, 1);
                    abs_L_k = sdpvar(1, length(new_L_k_vel));
                    L_k_matrix = [new_L_k_vel{:}];

                    constraints = [constraints, ...
                        abs_L_k >= L_k_matrix, ...
                        abs_L_k >= -L_k_matrix, ...
                        rad_k_dynamic == sum(abs_L_k)];
                end

                total_radius_k = rad_k_dynamic + obj.const_radius_d;

                constraints = [constraints, ...
                    nu_k + total_radius_k <= obj.v_max, ...
                    nu_k - total_radius_k >= -obj.v_max];

                % 5. Update state and generators for next loop
                mu_k = mu_k_next;

                % Create a new list for the next generators of the state
                L_k_next_state = {};

                for i = 1:length(L_k_state)
                    Q_z_i = [L_k_state{i}; zeros(obj.dim_g, 1)];
                    L_k_next_state{end+1} = obj.M_c * Q_z_i;
                end
                Q_k_B2 = diag(obj.Delta * abs_c_z);
                for i = 1:size(Q_k_B2, 2)
                    L_k_next_state{end+1} = Q_k_B2(:, i);
                end
                Q_k_B3 = diag(obj.Delta * rad_Rz_input);
                for i = 1:size(Q_k_B3, 2)
                    L_k_next_state{end+1} = Q_k_B3(:, i);
                end
                % (Disturbance)
                for i = 1:size(obj.Q_d, 2)
                    L_k_next_state{end+1} = obj.Q_d(:, i);
                end

                % Update the main list for the next iteration
                L_k_state = L_k_next_state;
            end

            % Terminal constraint 
            % T_j*c_g^(M_c) + delta_sigma,j(g) <= h_j - delta_inf - delta_k0
            % delta_sigma,j(g) is rendered linear through epigraph reformulation used for the robust constraints 
            if obj.use_terminal
                abs_g = sdpvar(obj.dim_g, 1);
                constraints = [constraints, abs_g >= g, abs_g >= -g];

                xc      = obj.Gam_Gc * g;               % nominal equilibrium state
                abs_xc  = sdpvar(obj.dim_x, 1);
                constraints = [constraints, abs_xc >= xc, abs_xc >= -xc];

                % xbar_g >= |xhat_g| for every realisation in B_cl
                xbar_g = obj.Wbar * (abs_xc + obj.absGam * obj.Delta_G * abs_g);

                delta_sigma = obj.absTjHc * obj.absGam * ...
                              (obj.Delta_Phi * xbar_g + obj.Delta_G * abs_g);

                h_term = obj.v_max - obj.delta_inf - obj.delta_k0;
                cg_val = obj.cg_slope * g;              % T_j*c_g^(M_c)

                constraints = [constraints, ...
                    cg_val + delta_sigma <=  h_term, ...
                    cg_val - delta_sigma >= -h_term];
            end

            objective = (r_t - g)' * obj.Psi * (r_t - g);
            options = sdpsettings('verbose', 0, 'solver', obj.solver_name);
            status = optimize(constraints, objective, options);

            if status.problem ~= 0
                warning('YALMIP: %s', status.info);
            end

            g = double(g);
        end
    end
end


