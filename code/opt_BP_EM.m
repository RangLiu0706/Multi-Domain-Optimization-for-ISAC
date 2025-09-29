function [SINR, W_jn, alpha, beta, Pj, pk] = opt_BP_EM(theta_AP, Hjnk, Gijn, Phiij, P, sigma2_r, sigma2_c, Nsub, XPD, Rc,alg)
% BP-EM optimization
% Fixed: AP selection (geometric uniform distribution)
% Dynamic: Interference-aware subcarrier allocation, polarization, beamforming

[K, ~, M2, Ns, J] = size(Hjnk); M = M2/2;

% Algorithm parameters
maxIter = alg.maxIter; res_thr = alg.res_thr;

% Initialization
pjm = [1; 1]/sqrt(2);
Pj = repmat(kron(eye(M), pjm), [1, 1, J]);
pk = repmat([1; 1]/sqrt(2), 1, K);

% ============= Stage 1: Resource Allocation =============
% Step 1: Fixed geometric uniform AP selection
alpha = ones(J, 1);
num_rx = floor(J/2);
step = floor(J / num_rx);
rx_indices = 1:step:J;
rx_indices = rx_indices(1:num_rx);
alpha(rx_indices) = 0;

% Step 2: Interference-aware subcarrier block allocation
beta_sub = allocate_subcarrier(Gijn, alpha, Nsub);

% Convert block allocation to full subcarrier allocation
Ng = Ns/Nsub; beta = zeros(J, Ns);
for j = 1:J
    beta(j,:) = kron(beta_sub(j,:), ones(1, Ng));
end
beta_tilde = 1 - sum(beta, 1);

% Evaluate initial configuration
He_nk = compute_He(Hjnk, pk, Pj, alpha, beta_tilde);
[Wc_jn, Pmin] = initialize_MRT_grouped(He_nk, sigma2_c, alpha, beta_tilde, Rc,Nsub);
if Pmin >= P
    fprintf('Warning: Cannot satisfy comm constraints with fixed AP selection, using fallback\n');
end

% Initialize radar beamforming
[Wr_jn, SINR] = initialize_Wr_grouped(theta_AP, Gijn, Phiij, P-Pmin, sigma2_r, XPD, Pj, alpha, beta,Nsub);

% ============= Stage 2: Joint Beamforming and Polarization Design =============
% Initialize combined beamforming matrix
W_jn = zeros(J*M, M+K, Ns);
for n = 1:Ns
    W_jn(:, 1:K, n) = Wc_jn(:, :, n);
    W_jn(:, K+1:end, n) = Wr_jn(:, :, n);
end

% Initialize auxiliary variables
t_nk = compute_t(He_nk, Wc_jn, sigma2_c, beta_tilde);
c_nk = compute_c(He_nk, Wc_jn, t_nk, sigma2_c);
[u_in, gamma_in] = update_u(Pj, theta_AP, Phiij, alpha, beta, Wr_jn, Gijn, sigma2_r, XPD);
eta_in = compute_eta(theta_AP, XPD, Wr_jn, Phiij, Pj, alpha, beta, u_in);

% Iterative optimization
obj = zeros(maxIter, 1); obj(1) = sum(gamma_in(:));

for it = 2:maxIter
    % Update polarizations
    pk = update_pk(Hjnk, Pj, W_jn, t_nk, c_nk);
    Pj = update_rPj(Wr_jn, alpha, beta, eta_in, u_in, Pj, Phiij, gamma_in, Gijn, theta_AP, XPD);
    Pj = update_tPj(Hjnk, pk, W_jn, t_nk, c_nk, alpha, beta, Pj, Rc, sigma2_c);

    % Update equivalent channel
    He_nk = compute_He(Hjnk, pk, Pj, alpha, beta_tilde);

    % Update beamforming
    W_jn = update_W(He_nk, t_nk, c_nk, P, eta_in, gamma_in, alpha, beta, Rc, ...
        sigma2_c, u_in, Pj, Phiij, theta_AP, XPD, Gijn, Nsub);

    % Extract communication and radar beamformers
    Wc_jn = zeros(J*M, K, Ns);
    Wr_jn = zeros(J*M, M, Ns);
    for n = 1:Ns
        Wc_jn(:, :, n) = W_jn(:, 1:K, n);
        Wr_jn(:, :, n) = W_jn(:, K+1:end, n);
    end

    % Update auxiliary variables
    eta_in = compute_eta(theta_AP, XPD, Wr_jn, Phiij, Pj, alpha, beta, u_in);
    [u_in, gamma_in] = update_u(Pj, theta_AP, Phiij, alpha, beta, Wr_jn, Gijn, sigma2_r, XPD);
    t_nk = compute_t(He_nk, W_jn, sigma2_c, beta_tilde);
    c_nk = compute_c(He_nk, Wc_jn, t_nk, sigma2_c);

    % Check convergence
    obj(it) = sum(gamma_in(:));
    if it >= 6
        current = mean(obj(it-2:it));
        previous = mean(obj(it-5:it-3));
        rel_change = abs(current - previous) / (abs(previous) + 1e-10);
        if rel_change < res_thr
            break;
        end
    end
end
% figure;plot(obj(1:it)); grid on;
SINR = sum(gamma_in(:));
end
