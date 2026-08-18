function [Lambda_k_j, J_bcl_k_j, J_mcl_k_j, Z_bcl_all, Z_mcl_all] = ...
    analyze_conservatism(M_cl, box_B_cl, W, x_t, g, k0)
% ANALYZE_CONSERVATISM - Computes conservatism and returns reachable sets
%
% This function implements the simulation described in
% Section VII of the ZRCG paper to quantify conservatism.
%
% Please notice that the M_cl method (N_mcl) grows exponentially.
% Depending on the horizon, running this function may take an extremely long
% time or crash MATLAB due to memory limitations.

T_j = [0 1 0]; 

int_obj = box_B_cl.int;
L = int_obj.inf;
U = int_obj.sup;
M_c = (L + U) / 2;
Delta = (U - L) / 2;

C_cl = M_cl.center;
H_cl_list = M_cl.generator;
m_gen = length(H_cl_list);

c_d = W.center;
Q_d = W.Z(:, 2:end);
v_gen = size(Q_d, 2);

% Now J_{S^m} and J_{\hat{S}}
J_bcl_k_j = zeros(k0, 1);
J_mcl_k_j = zeros(k0, 1);

% Used to store zonotope objects
Z_bcl_all = cell(k0, 1);
Z_mcl_all = cell(k0, 1);

mu_k_bcl = x_t;
L_k_state_bcl = {};
mu_k_mcl = x_t;
L_k_state_mcl = {}; 

fprintf('Running conservatism analysis for k=1:%d...\n', k0);
%fprintf('WARNING: Classic (MCL) method has O((p_n+1)^k) generator growth.\n');

dim_g = length(g);

for k = 1:k0
    % (BCL - Lemma 2)
    c_z_bcl = [mu_k_bcl; g];
    mu_k_bcl_next = M_c * c_z_bcl + c_d;
    
    if isempty(L_k_state_bcl)
        rad_Rz_state_bcl = zeros(size(x_t));
    else
        rad_Rz_state_bcl = sum(abs([L_k_state_bcl{:}]), 2);
    end
    rad_Rz_input_bcl = [rad_Rz_state_bcl; zeros(dim_g, 1)];

    L_k_next_state_bcl = {};
    
    % B1
    for i = 1:length(L_k_state_bcl)
        Q_z_i = [L_k_state_bcl{i}; zeros(dim_g, 1)];
        L_k_next_state_bcl{end+1} = M_c * Q_z_i;
    end
    % B2
    Q_k_B2 = diag(Delta * abs(c_z_bcl));
    for i = 1:size(Q_k_B2, 2)
        L_k_next_state_bcl{end+1} = Q_k_B2(:, i);
    end
    % B3
    Q_k_B3 = diag(Delta * rad_Rz_input_bcl);
    for i = 1:size(Q_k_B3, 2)
        L_k_next_state_bcl{end+1} = Q_k_B3(:, i);
    end
    % B4
    for i = 1:v_gen
        L_k_next_state_bcl{end+1} = Q_d(:, i);
    end
    
    % Store the full zonotope object
    Z_bcl_all{k} = zonotope(mu_k_bcl_next, [L_k_next_state_bcl{:}]);
    
    % Compute Support Function
    nu_k_bcl = mu_k_bcl_next(2);
    rad_k_bcl = sum(abs(T_j * [L_k_next_state_bcl{:}]));
    J_bcl_k_j(k) = nu_k_bcl + rad_k_bcl;
    
    mu_k_bcl = mu_k_bcl_next;
    L_k_state_bcl = L_k_next_state_bcl;

    
    % (MCL)
    c_z_mcl = [mu_k_mcl; g];
    mu_k_mcl_next = C_cl * c_z_mcl + c_d;
    L_k_next_state_mcl = {};
    
    % B1
    for i = 1:length(L_k_state_mcl)
        Q_z_i = [L_k_state_mcl{i}; zeros(dim_g, 1)];
        L_k_next_state_mcl{end+1} = C_cl * Q_z_i;
    end
    % B2
    for i = 1:m_gen
        L_k_next_state_mcl{end+1} = H_cl_list{i} * c_z_mcl;
    end
    % B3
    for i = 1:m_gen
        for j = 1:length(L_k_state_mcl)
            Q_z_j = [L_k_state_mcl{j}; zeros(dim_g, 1)];
            L_k_next_state_mcl{end+1} = H_cl_list{i} * Q_z_j;
        end
    end
    % B4
    for i = 1:v_gen
        L_k_next_state_mcl{end+1} = Q_d(:, i);
    end
    
    % Store the full zonotope object
    Z_mcl_all{k} = zonotope(mu_k_mcl_next, [L_k_next_state_mcl{:}]);
    
    % Compute Support Function
    nu_k_mcl = mu_k_mcl_next(2);
    rad_k_mcl = sum(abs(T_j * [L_k_next_state_mcl{:}]));
    J_mcl_k_j(k) = nu_k_mcl + rad_k_mcl;
    
    mu_k_mcl = mu_k_mcl_next;
    L_k_state_mcl = L_k_next_state_mcl;
    
    fprintf('... completed step k=%d. (N_bcl=%d, N_mcl=%d)\n', ...
            k, length(L_k_state_bcl), length(L_k_state_mcl));
end

Lambda_k_j = J_bcl_k_j - J_mcl_k_j;
fprintf('Analysis complete.\n');

end

