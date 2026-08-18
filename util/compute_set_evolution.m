function Z_bcl_all = compute_set_evolution(box_B_cl, W, x_t, g, k0)
% COMPUTE_SET_EVOLUTION - Computes k0 reachable set evolution (Lemma 2)

% INPUTS:
%   box_B_cl: The interval hull <M_c, Delta>
%   W:        The disturbance zonotope <c_d, Q_d>
%   x_t:      The initial state 
%   g:        The fixed command to apply 
%   k0:       The prediction horizon 
%
% OUTPUTS:
%   Z_bcl_all: computed zonotope objects


int_obj = box_B_cl.int;
L = int_obj.inf;
U = int_obj.sup;
M_c = (L + U) / 2;
Delta = (U - L) / 2; 

c_d = W.center;
Q_d = W.Z(:, 2:end);
v_gen = size(Q_d, 2);

% Used to store zonotope objects
Z_bcl_all = cell(k0, 1); 

mu_k_bcl = x_t;
L_k_state_bcl = {}; % generator vectors

fprintf('Running set evolution for k=1:%d...\n', k0);
dim_g = length(g);

for k = 1:k0
    % (BCL - Lemma 2)
    c_z_bcl = [mu_k_bcl; g];
    
    % Center
    mu_k_bcl_next = M_c * c_z_bcl + c_d;
    
    % Generators
    if isempty(L_k_state_bcl)
        rad_Rz_state_bcl = zeros(size(x_t));
    else
        rad_Rz_state_bcl = sum(abs([L_k_state_bcl{:}]), 2);
    end
    rad_Rz_input_bcl = [rad_Rz_state_bcl; zeros(dim_g, 1)];

    L_k_next_state_bcl = {};
    

    for i = 1:length(L_k_state_bcl)
        Q_z_i = [L_k_state_bcl{i}; zeros(dim_g, 1)];
        L_k_next_state_bcl{end+1} = M_c * Q_z_i;
    end

    Q_k_B2 = diag(Delta * abs(c_z_bcl));
    for i = 1:size(Q_k_B2, 2)
        L_k_next_state_bcl{end+1} = Q_k_B2(:, i);
    end

    Q_k_B3 = diag(Delta * rad_Rz_input_bcl);
    for i = 1:size(Q_k_B3, 2)
        L_k_next_state_bcl{end+1} = Q_k_B3(:, i);
    end

    for i = 1:v_gen
        L_k_next_state_bcl{end+1} = Q_d(:, i);
    end
    
    % Store the full zonotope 
    Z_bcl_all{k} = zonotope(mu_k_bcl_next, [L_k_next_state_bcl{:}]);
    
    % Update state
    mu_k_bcl = mu_k_bcl_next;
    L_k_state_bcl = L_k_next_state_bcl;
    
    fprintf('... completed step k=%d. (N_bcl=%d)\n', ...
            k, length(L_k_state_bcl));
end

fprintf('Set evolution computation complete.\n');

end

