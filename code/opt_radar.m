function [SINR, W_jn, alpha, beta, Pj] = opt_radar(theta_AP, Gijn, Phiij, P, sigma2_r, XPD, Nsub,alg)
% Optimize radar performance without communication constraints

[M2, ~, Ns, J, ~] = size(Gijn); M = M2/2;

% Algorithm parameters
maxIter = alg.maxIter; res_thr = alg.res_thr;

% Initialization
pjm = [1; 1]/sqrt(2);
Pj = repmat(kron(eye(M), pjm), [1, 1, J]);

% ============= Stage 1: Resource Allocation =============
% Step 1: Bistatic-aware AP selection
alpha = allocate_ap_radar(Gijn, Phiij, theta_AP, P, sigma2_r, XPD, Nsub);
% Step 2: Subcarrier block allocation
beta_sub = allocate_subcarrier_radar(Gijn, alpha, Nsub);

% Convert block allocation to full subcarrier allocation
Ng = Ns/Nsub; beta = zeros(J, Ns);
for j = 1:J
    beta(j,:) = kron(beta_sub(j,:), ones(1, Ng));
end

% Initialize radar beamforming
[W_jn, SINR] = initialize_Wr_grouped(theta_AP, Gijn, Phiij, P, sigma2_r, XPD, Pj, alpha, beta,Nsub);

% ============= Stage 2: Joint Beamforming and Polarization Design =============
% Extract radar beamforming component
Wr_jn = W_jn;  % For radar-only, W_jn contains only radar beamforming

% Initialize auxiliary variables
[u_in, gamma_in] = update_u(Pj, theta_AP, Phiij, alpha, beta, Wr_jn, Gijn, sigma2_r, XPD);
eta_in = compute_eta(theta_AP, XPD, Wr_jn, Phiij, Pj, alpha, beta, u_in);

% Iterative optimization
obj = zeros(maxIter, 1); obj(1) = sum(gamma_in(:));
for it = 2:maxIter
    % Update polarizations
    Pj = update_rPj(Wr_jn, alpha, beta, eta_in, u_in, Pj, Phiij, gamma_in, Gijn, theta_AP, XPD);
    Pj = update_tPj_radar(Wr_jn, alpha, beta, eta_in, u_in, Pj, Phiij, gamma_in, Gijn, theta_AP, XPD);

    % Update beamforming
    W_jn = update_W_radar(P, eta_in, gamma_in, alpha, beta, u_in, Pj, Phiij, theta_AP, XPD, Gijn, Nsub);
    Wr_jn = W_jn;

    % Update auxiliary variables
    eta_in = compute_eta(theta_AP, XPD, Wr_jn, Phiij, Pj, alpha, beta, u_in);
    [u_in, gamma_in] = update_u(Pj, theta_AP, Phiij, alpha, beta, Wr_jn, Gijn, sigma2_r, XPD);

    % Check convergence
    obj(it) = sum(gamma_in(:));
    if it >= 6
        current = mean(obj(it-2:it)); previous = mean(obj(it-5:it-3));
        rel_change = abs(current - previous) / (abs(previous) + 1e-10);
        if rel_change < res_thr
            break;
        end
    end
end

SINR = sum(gamma_in(:));
end

% ==================== Helper Functions ====================
function Pj = update_tPj_radar(Wr_jn, alpha, beta, eta_in, u_in, Pj, Phiij, gamma_in, Gijn, theta_AP, XPD)
% Update polarization at Tx APs
[J, Ns] = size(beta); M = size(Pj, 2);
p = zeros(2*J*M, 1);
for j = 1:J
    for m = 1:M
        p(2*M*(j-1)+2*(m-1)+1:2*M*(j-1)+2*m, 1) = Pj(2*(m-1)+1:2*m, m, j);
    end
end

idxTx = find(alpha); idxRx = find(~alpha);

a = zeros(2*J*M, 1);
for j = idxTx.'
    aj = zeros(2*M, 1);
    Atj = kron(exp(1j*pi*sin(-theta_AP(j))*(0:M-1).'), XPD);
    for i = idxRx.'
        Ari = kron(exp(1j*pi*sin(-theta_AP(i))*(0:M-1).'), XPD);
        for n = 1:Ns
            idx = find(beta(:,n));
            if isempty(idx) || idx ~= j
                continue;
            end

            aijn = 2*u_in(:,i,n)' * Pj(:,:,i).' * Ari * Phiij(:,:,J*(i-1)+j) * Atj.';
            aj = aj + kron(Wr_jn((j-1)*M+1:j*M, :, n)*eta_in(:,i,n), [1;1]) .* aijn.';
        end
    end

    a((j-1)*2*M+1:j*2*M, 1) = aj;
end
a = real(a);

B = zeros(2*J*M, 2*J*M);
for i = idxRx.'
    for n = 1:Ns
        j = find(beta(:,n));
        if isempty(j)
            continue;
        end
        bijn = sqrt(gamma_in(i,n)) * u_in(:,i,n)' * Pj(:,:,i).' * Gijn(:,:,n,i,j);
        Bijn = sparse(1:2*M, repelem(1:M, 2), bijn(:), 2*M, M);
        Tj = zeros(2*M, 2*J*M);
        Tj(:, (j-1)*2*M+1:j*2*M) = eye(2*M);
        Bin = Tj.' * Bijn * Wr_jn((j-1)*M+1:j*M, :, n);
        B = B + Bin * Bin';
    end
end
B = real(B+B')/2;

s1 = max(abs(a));
if s1 > 0
    a = a/s1; B = B/s1;
end

mask1 = logical(kron(alpha, ones(2*M,1)));
p1 = p(mask1); p1t = p(~mask1);
a1 = a(mask1); B1 = B(mask1, mask1);
temp = 2*p1t.' * B(~mask1, mask1);
a1 = a1 - temp.';
Mt2 = size(a1, 1); Mt = Mt2/2;

manifold = obliquefactory(2, Mt);
problem.M = manifold;
problem.cost = @(P) (P(:)).' * B1 * P(:) - real(a1.' * P(:));
problem.egrad = @(P) 2*reshape(B1*P(:), 2, Mt) - reshape(a1, 2, Mt);
options.verbosity = 0;

[P_opt, ~, ~] = conjugategradient(problem, reshape(p1, 2, Mt), options);
p1 = P_opt(:); p(mask1) = p1;
for j = 1:J
    if alpha(j) == 1
        for m = 1:M
            Pj(2*(m-1)+1:2*m, m, j) = p(2*M*(j-1)+2*(m-1)+1:2*M*(j-1)+2*m, 1);
        end
    end
end
end

function Wjn = update_W_radar(Ptot, eta_in, gamma_in, alpha, beta, u_in, Pj, Phiij, theta_AP, XPD, Gijn, Nsub)
% Update beamforming matrices
[J, Ns] = size(beta); M = size(Pj, 2); Ng = Ns/Nsub; Mt = J*M;

idxRx = find(~alpha); jtx = zeros(Nsub, 1);
An = zeros(M, M, Nsub); Bn = zeros(M, M, Nsub);
for g = 1:Nsub
    idx = (g-1)*Ng+1:g*Ng;
    idxTx = find(beta(:,g*Ng));
    if isempty(idxTx)
        continue;
    end
    jtx(g) = idxTx;
    Atj = kron(exp(1j*pi*sin(-theta_AP(idxTx))*(0:M-1).'), XPD);
    Temp = zeros(M, M);
    for irx = idxRx.'
        for n = idx
            Ari = kron(exp(1j*pi*sin(-theta_AP(irx))*(0:M-1).'), XPD);
            An(:,:,g) = An(:,:,g) + 2*eta_in(:,irx,n) * u_in(:,irx,n)' * ...
                Pj(:,:,irx).' * Ari * Phiij(:,:,(irx-1)*J+idxTx) * ...
                Atj.' * Pj(:,:,idxTx);

            ain = sqrt(gamma_in(irx,n)) * u_in(:,irx,n)' * Pj(:,:,irx).' * ...
                Gijn(:,:,n,irx,idxTx) * Pj(:,:,idxTx);
            Temp = Temp + ain' * ain;
        end
    end

    [Ve, De] = eigs(Temp, numel(idxRx));
    Bn(:,1:numel(idxRx),g) = Ve * sqrt(De);
end

scale = max(abs(Bn(:)));
if scale > max(abs(An(:)))/1e3
    An = An/(scale^2); Bn = Bn/scale;
else
    scale = max(abs(An(:)));
    An = An/scale; Bn = zeros(M, M, Nsub);
end

cvx_begin quiet
variable Wg(Mt, M, Nsub) complex
expression pow
expression ff
pow = 0;
ff = 0;
for g = 1:Nsub
    if jtx(g) > 0
        pow = pow + square_pos(norm(Wg((jtx(g)-1)*M+1:jtx(g)*M, :, g), 'fro'));
        ff = ff + real(trace(An(:,:,g) * Wg((jtx(g)-1)*M+1:jtx(g)*M, :, g))) - ...
            square_pos(norm(Bn(:,:,g)' * Wg((jtx(g)-1)*M+1:jtx(g)*M, :, g), 'fro'));
    end
end
maximize ff
subject to
pow <= Ptot/Ng;
cvx_end

Wjn = zeros(Mt, M, Ns);
for g = 1:Nsub
    idx = (g-1)*Ng+1:g*Ng;
    Wjn(:, :, idx) = repmat(Wg(:,:,g), 1, 1, Ng);
end
end

function alpha = allocate_ap_radar(Gijn, Phiij, theta_AP, P, sigma2_r, XPD, Nsub)
% AP selection for radar-only system
[M2, ~, ~, J, ~] = size(Gijn); M = M2/2;

% Initialize: all APs as Tx
alpha = ones(J, 1);
current_SINR = 0;
best_alpha = alpha;
best_SINR = 0;

% Fixed circular polarization
pjm = [1; 1]/sqrt(2);
Pj = repmat(kron(eye(M), pjm), [1, 1, J]);

min_tx = 2;  % Keep at least 2 Tx
max_rx = J - min_tx;

% Convert Tx to Rx incrementally
for num_rx = 1:max_rx
    tx_indices = find(alpha == 1);
    if length(tx_indices) <= min_tx
        break;
    end

    % Evaluate conversions
    scores = zeros(length(tx_indices), 1);
    for i = 1:length(tx_indices)
        alpha_temp = alpha;
        alpha_temp(tx_indices(i)) = 0;

        % Allocate subcarriers for new configuration
        beta_sub = allocate_subcarrier_radar(Gijn, alpha_temp, Nsub);

        % Convert to full allocation
        Ng = size(Gijn, 3) / Nsub;
        beta = zeros(J, size(Gijn, 3));
        for j = 1:J
            beta(j,:) = kron(beta_sub(j,:), ones(1, Ng));
        end

        % Compute radar SINR
        [~, SINR] = initialize_Wr_grouped(theta_AP, Gijn, Phiij, P, sigma2_r, XPD, Pj, alpha_temp, beta,Nsub);
        scores(i) = SINR;
    end

    % Select best conversion
    [best_score, best_idx] = max(scores);

    % Update if improvement
    if best_score > best_SINR * 1.001  %
        alpha(tx_indices(best_idx)) = 0;
        best_SINR = best_score;
        best_alpha = alpha;
    else
        % No improvement, restore best configuration and stop
        alpha = best_alpha;
        break;
    end
end
end

function beta_sub = allocate_subcarrier_radar(Gijn, alpha, Nsub)
% Subcarrier allocation for radar-only system
% All blocks allocated to radar with round-robin distribution

[~, ~, Ns, J, ~] = size(Gijn);
Ng = Ns/Nsub;

idxTx = find(alpha == 1);
idxRx = find(alpha == 0);
n_tx = length(idxTx);
n_rx = length(idxRx);

% Initialize
beta_sub = zeros(J, Nsub);

% Handle edge cases
if isempty(idxRx) || isempty(idxTx)
    return;  % Cannot perform radar without both Tx and Rx
end

% Compute interference table
InterferenceTable = zeros(n_tx, Nsub);
for b = 1:Nsub
    sc_indices = (b-1)*Ng+1:b*Ng;
    for tx_idx = 1:n_tx
        j_tx = idxTx(tx_idx);
        total_interf = 0;
        for n = sc_indices
            for i_rx = idxRx'
                total_interf = total_interf + norm(Gijn(:,:,n,i_rx,j_tx), 'fro')^2;
            end
        end
        InterferenceTable(tx_idx, b) = total_interf / (Ng * n_rx);
    end
end

% Sort blocks by average interference (optional for radar-only)
block_avg_interf = mean(InterferenceTable, 1);
[~, sorted_blocks] = sort(block_avg_interf);

% All blocks for radar
radar_blocks = sorted_blocks;  % Use all Nsub blocks

% Round-robin allocation for perfect load balancing
for idx = 1:Nsub
    b = radar_blocks(idx);
    tx_idx = mod(idx-1, n_tx) + 1;  % Cycle through Tx APs
    beta_sub(idxTx(tx_idx), b) = 1;
end

end

