% First define true system and physical disturbance W
% Second, collect data via standard simulation and compute M_cl (Lemma 1).
% Then visualize the model space in 2D and 3D (M_cl vs B_cl).
% Finally propagate the state-space reachable sets over multiple steps.


clear; clc; close all;

%% =========================================
% (Lemma 1)
% ==========================================
disp('Data Collection & Building M_cl');

% True unknown system
Phi_true = [0.5, 0.1; 
           -0.2, 0.4]; 
G_true   = [0.5; 
            0.3];
dim_x = 2; 
dim_v = 1;

% Disturbance Zonotope (W)
G_base = [0.3, -0.1; 
          0.1,  0.3];
rng(42); % Seed for reproducible shape
G_fatten = 0.05 * randn(2, 20); % 20 gens to keep M_cl size manageable

scale_factor = 0.02; 
G_noise = [G_base, G_fatten] * scale_factor;
W_dist = zonotope([0;0], G_noise); % True physical process noise

% Data Collection 
N_s = 10;
X_minus = zeros(dim_x, N_s);
X_plus  = zeros(dim_x, N_s);
U_minus = zeros(dim_v, N_s);

x_curr = [1.0; -1.0]; % Initial state for data collection
rng(100); % Seed for random inputs and noise realizations

% Data collection 
for i = 1:N_s
    u_curr = 15*rand(dim_v, 1) - 1; 
    w_curr = randPoint(W_dist);    

    x_next = Phi_true * x_curr + G_true * u_curr + w_curr;

    X_minus(:, i) = x_curr;
    U_minus(:, i) = u_curr;
    X_plus(:, i)  = x_next;

    x_curr = x_next; 
end

Z_data = [X_minus; U_minus];

% Disturbance Matrix Zonotope M_w
num_w_gens = size(W_dist.Z, 2) - 1;
c_w_mat = zeros(dim_x, N_s);
g_w_mat = cell(1, num_w_gens * N_s);

gen_idx = 1;
for i = 1:num_w_gens
    gen_vec = W_dist.Z(:, i+1);
    for j = 1:N_s
        temp_gen = zeros(dim_x, N_s);
        temp_gen(:, j) = gen_vec;
        g_w_mat{gen_idx} = temp_gen;
        gen_idx = gen_idx + 1;
    end
end
M_w = matZonotope(c_w_mat, g_w_mat);

% Lemma 1
M_cl_matZ = (X_plus + -1*M_w) * pinv(Z_data);
C_cl = center(M_cl_matZ);

% Extract generators
if isprop(M_cl_matZ, 'generator')
    H_cl_raw = M_cl_matZ.generator;
elseif isprop(M_cl_matZ, 'G')
    H_cl_raw = M_cl_matZ.G;
end

if ~iscell(H_cl_raw)
    num_gens = size(H_cl_raw, 3);
    H_cl_list = cell(num_gens, 1);
    for i = 1:num_gens
        H_cl_list{i} = H_cl_raw(:,:,i);
    end
else
    H_cl_list = H_cl_raw;
    num_gens = length(H_cl_list);
end

% Interval Hull B_cl 
Delta = zeros(size(C_cl));
for i = 1:num_gens
    Delta = Delta + abs(H_cl_list{i});
end
disp(['Successfully built data-driven M_cl with ', num2str(num_gens), ' generators.']);
%disp('Notice that the center of M_cl slightly differs from the True system due to realized noise!');


%% =========================================
% Figure 1 (2D projection, Phi_{1,1} and Phi_{2,1})
% ==========================================

C_plot_2d = [C_cl(1,1); C_cl(2,1)];
G_plot_2d = zeros(2, num_gens);
for i = 1:num_gens
    G_plot_2d(:, i) = [H_cl_list{i}(1,1); H_cl_list{i}(2,1)];
end

Z_mcl_2d = zonotope(C_plot_2d, G_plot_2d);
Z_bcl_2d = interval(Z_mcl_2d);

figure('Name', 'Model Space (2D)', 'Color', 'w', 'Position', [100, 100, 500, 450]);
hold on; 

plot(Z_bcl_2d, [1 2], 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r', 'LineWidth', 1.5);
plot(Z_mcl_2d, [1 2], 'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'b', 'LineWidth', 1.5);
plot(C_plot_2d(1), C_plot_2d(2), 'kx', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('Matrix Entry $\Phi_{1,1}$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Matrix Entry $\Phi_{2,1}$', 'Interpreter', 'latex', 'FontSize', 12);
title('Model Space: 2D Projection ($\mathcal{M}_{cl}$ vs $\mathcal{B}_{cl}$)', 'Interpreter', 'latex');
legend('Interval Hull $\mathcal{B}_{cl}$', 'Matrix Zonotope $\mathcal{M}_{cl}$', ...
       'Identified Center $M_c$', 'Interpreter', 'latex', 'Location', 'best');
grid on; axis equal;
xlim([Z_bcl_2d.inf(1) - 0.005,  Z_bcl_2d.sup(1) + 0.005])
ylim([Z_bcl_2d.inf(2) - 0.005,  Z_bcl_2d.sup(2) + 0.005])


area_mcl = volume(Z_mcl_2d);
area_bcl = volume(zonotope(Z_bcl_2d));
r_B = area_bcl / area_mcl;

disp(['Projection-Area Ratio (r_{BM})']);
disp(['Area of M_cl projection: ', num2str(area_mcl)]);
disp(['Area of B_cl projection: ', num2str(area_bcl)]);
disp(['Ratio r_BM: ', num2str(r_B)]);



%% =========================================
% Figure 2 (Predicted sets in x_1 and x_2)
% ==========================================
disp('Propagate Reachable Sets');
x_t = [0.0; 0.0];  
g_cmd = [0.25];     
% Initialize state Zonotopes
Z_exact = zonotope(x_t);
Z_proposed = zonotope(x_t);
Z_reduced = zonotope(x_t);
max_order = 15;  % Maximum zonotope order
K_steps = 3;    

figure('Name', 'State Space Evolution', 'Color', 'w', 'Position', [150, 400, 650, 550]);
hold on;
centers_x = x_t(1);
centers_y = x_t(2);
plot(x_t(1), x_t(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Initial State $x_t$ ($\mathcal{R}_0$)');

for k = 1:K_steps
    fprintf('\n--- Computing step k = %d ---\n', k);
    
    % Augment state sets with the command g
    c_ex = center(Z_exact); G_ex = generators(Z_exact);
    Z_aug_exact = zonotope([[c_ex; g_cmd], [G_ex; zeros(1, size(G_ex,2))]]);
    
    c_pr = center(Z_proposed); G_pr = generators(Z_proposed);
    Z_aug_prop = zonotope([[c_pr; g_cmd], [G_pr; zeros(1, size(G_pr,2))]]);
    
    c_rd = center(Z_reduced); G_rd = generators(Z_reduced);
    Z_aug_red = zonotope([[c_rd; g_cmd], [G_rd; zeros(1, size(G_rd,2))]]);
    
    % 1 - Exact Classic Pipeline 
    Z_exact = exact_enclosure(C_cl, H_cl_list, Z_aug_exact, W_dist);
    
    % 2 - Proposed Pipeline (Lemma 2: B_cl)
    Z_proposed = lemma2_enclosure(C_cl, Delta, Z_aug_prop, W_dist);
    
    % 3 - Reduced Pipeline (For comparison)
    Z_reduced_temp = exact_enclosure(C_cl, H_cl_list, Z_aug_red, W_dist);
    Z_reduced = reduce(Z_reduced_temp, 'girard', max_order);
    
    % --- Volume Comparison ---
    vol_ex = exact_2D_zonotope_volume(Z_exact);
    vol_rd = exact_2D_zonotope_volume(Z_reduced);
    vol_pr = exact_2D_zonotope_volume(Z_proposed);
    
    inc_rd = ((vol_rd - vol_ex) / vol_ex) * 100;
    inc_pr = ((vol_pr - vol_ex) / vol_ex) * 100;
    
    fprintf('Exact Volume: %.4f\n', vol_ex);
    fprintf('Reduced Pipeline Overestimation: %.2f%%\n', inc_rd);
    fprintf('Proposed Pipeline Overestimation: %.2f%%\n', inc_pr);
    
    c_new = center(Z_proposed);
    centers_x(end+1) = c_new(1);
    centers_y(end+1) = c_new(2);
    
    % --- Plotting ---
    h_prop = plot(Z_proposed, [1 2], 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r', 'LineWidth', 1.5);
    h_red = plot(Z_reduced, [1 2], 'FaceColor', 'none', 'EdgeColor', [0 0.5 0], 'LineWidth', 1.5, 'LineStyle', '-.');
    h_ex = plot(Z_exact, [1 2], 'FaceColor', 'none', 'EdgeColor', 'b', 'LineWidth', 1.5);
    if k > 1
        set(h_prop, 'HandleVisibility', 'off');
        set(h_red, 'HandleVisibility', 'off');
        set(h_ex, 'HandleVisibility', 'off');
    end
end

plot(centers_x, centers_y, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Nominal Trajectory Path');
xlabel('State $x_1$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('State $x_2$', 'Interpreter', 'latex', 'FontSize', 12);
title('State Space: Evolution of Reachable Sets ($k=1 \dots 3$)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Initial State $x_t$ ($\mathcal{R}_0$)', ...
       'Proposed $\mathcal{B}_{cl}$ ($\mathcal{O}(k)$)', ...
       'Reduced $\mathcal{M}_{cl}$ (Needs tuning)', ...
       'Exact $\mathcal{M}_{cl}$ (Intractable)', ...
       'Nominal Trajectory Path', ...
       'Interpreter', 'latex', 'Location', 'northwest');
grid on; axis equal;
disp('Execution completed successfully.');





% Useful functions

function Z_next = exact_enclosure(C_cl, H_cl_list, Z_prev, W)
    c_k_1 = center(Z_prev); Q_k_1 = generators(Z_prev);
    N_k_1 = size(Q_k_1, 2); p_n = length(H_cl_list);
    
    c_new = C_cl * c_k_1;
    G_1 = C_cl * Q_k_1;
    
    G_2 = zeros(size(C_cl, 1), p_n);
    for i = 1:p_n
        G_2(:, i) = H_cl_list{i} * c_k_1;
    end
    
    G_3 = zeros(size(C_cl, 1), p_n * N_k_1);
    col_idx = 1;
    for i = 1:p_n
        G_3(:, col_idx : col_idx + N_k_1 - 1) = H_cl_list{i} * Q_k_1;
        col_idx = col_idx + N_k_1;
    end
    
    Z_next = zonotope([c_new + center(W), [G_1, G_2, G_3, generators(W)]]);
end

function Z_next = lemma2_enclosure(M_c, Delta, Z_prev, W)
    c_prev = center(Z_prev); Q_prev = generators(Z_prev);
    c_new = M_c * c_prev;
    Q_1 = M_c * Q_prev;
    Q_2 = diag(Delta * abs(c_prev));
    Q_3 = diag(Delta * sum(abs(Q_prev), 2));
    
    Z_next = zonotope([c_new + center(W), [Q_1, Q_2, Q_3, generators(W)]]);
end



function vol = exact_2D_zonotope_volume(Z)
    % Compute exact 2D zonotope volume (polygon perimeter)
    G = generators(Z);
    if isempty(G)
        vol = 0;
        return;
    end
    
    % Extract first two dimensions
    G = G(1:2, :);
    
    % Sort all generators and their negations by angle so a 2D zonotope
    % boundary can be traced
    G_all = [G, -G];
    angles = atan2(G_all(2, :), G_all(1, :));
    
    % Sort vectors radially
    [~, sort_idx] = sort(angles);
    G_sorted = G_all(:, sort_idx);
    
    % Cumulative sum to generate perimeter vertices (head to tail)
    vertices = cumsum(G_sorted, 2);
    
    % Compute volume using standard polygon area
    vol = polyarea(vertices(1, :), vertices(2, :));
end




