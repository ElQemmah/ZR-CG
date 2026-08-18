function Z_next = exact_paper_enclosure(M_cl, Z_prev, W)
    % Implements the exact unnumbered block matrix equation below Eq (29)
    
    c_k_1 = center(Z_prev);
    Q_k_1 = generators(Z_prev);
    N_k_1 = size(Q_k_1, 2);
    
    C_cl = center(M_cl);
    
    % Generators of the matZonotope
    if isprop(M_cl, 'generator')
        H_cl = M_cl.generator;
    elseif isprop(M_cl, 'G')
        H_cl = M_cl.G;
    else
        error('Could not find the generator property for M_cl.');
    end
    
    % In some versions of CORA, H_cl might be a 3D tensor instead of a cell. 
    % This is useful to ensure it behaves like a cell array of 2D matrices.
    if ~iscell(H_cl)
        num_gens = size(H_cl, 3);
        temp_cell = cell(1, num_gens);
        for i = 1:num_gens
            temp_cell{i} = H_cl(:, :, i);
        end
        H_cl = temp_cell;
    end
    
    p_n = length(H_cl);
    
    % Center
    c_new = C_cl * c_k_1;
    
    % C_cl * Q_{k-1}
    G_1 = C_cl * Q_k_1;
    
    % H_{cl}^{(i)} * c_{k-1}  (for i = 1...p_n)
    G_2 = zeros(size(C_cl, 1), p_n);
    for i = 1:p_n
        G_2(:, i) = H_cl{i} * c_k_1;
    end
    
    % Cross block H_{cl}^{(i)} * q_{k-1}^{(j)}
    G_3 = zeros(size(C_cl, 1), p_n * N_k_1);
    col_idx = 1;
    for i = 1:p_n
        % Multiplication (using vectors) for all q in H_i 
        G_3(:, col_idx : col_idx + N_k_1 - 1) = H_cl{i} * Q_k_1;
        col_idx = col_idx + N_k_1;
    end
    
    % Add Disturbance W
    c_w = center(W);
    Q_w = generators(W);
    
    % Final zonotope
    c_final = c_new + c_w;
    G_final = [G_1, G_2, G_3, Q_w];
    
    Z_next = zonotope([c_final, G_final]);
end

