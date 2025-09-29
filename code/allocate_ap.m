function alpha = allocate_ap(Hjnk, Gijn, Phiij, theta_AP, P, sigma2_r, sigma2_c, Nsub, XPD, Rc, K, M)
% AP selection for ISAC system using greedy approach

J = size(Hjnk, 5);
Ns = size(Hjnk, 4);
min_tx = max(ceil(K/M),2);  % Minimum Tx needed for communication

% Initialize: all APs as Tx
alpha = ones(J, 1);

% Fixed circular polarization for evaluation
pjm = [1; 1]/sqrt(2);
Pj = repmat(kron(eye(M), pjm), [1, 1, J]);
pk = repmat([1; 1]/sqrt(2), 1, K);

% Compute initial SINR (all Tx, no Rx = no radar)
current_SINR = 0;  % No radar initially

% Main optimization loop
for iter = 1:(J - min_tx)
    tx_indices = find(alpha == 1);

    % Check constraint
    if length(tx_indices) <= min_tx
        break;
    end

    % Evaluate each possible Tx->Rx conversion
    scores = zeros(length(tx_indices), 1);

    for i = 1:length(tx_indices)
        % Try converting tx_indices(i) to Rx
        alpha_temp = alpha;
        alpha_temp(tx_indices(i)) = 0;

        % Check if we have at least one Rx
        if sum(alpha_temp == 0) == 0
            scores(i) = -inf;  % No Rx means no radar
            continue;
        end

        % Allocate subcarriers for the NEW configuration
        beta_sub = allocate_subcarrier(Gijn, alpha_temp, Nsub);

        % Convert block allocation to full subcarrier allocation
        Ng = Ns/Nsub;
        beta = zeros(J, Ns);
        for j = 1:J
            beta(j,:) = kron(beta_sub(j,:), ones(1, Ng));
        end
        beta_tilde = 1 - sum(beta, 1);

        % Check communication constraint
        He_nk = compute_He(Hjnk, pk, Pj, alpha_temp, beta_tilde);
        [~, Pmin] = initialize_MRT_grouped(He_nk, sigma2_c, alpha, beta_tilde, Rc,Nsub);

        if Pmin >= P
            scores(i) = -inf;  % Cannot meet communication requirement
            continue;
        end

        % Compute radar SINR with remaining power
        [~, SINR] = initialize_Wr_grouped(theta_AP, Gijn, Phiij, P-Pmin, ...
            sigma2_r, XPD, Pj, alpha_temp, beta,Nsub);
        scores(i) = SINR;
    end

    % Select best conversion
    [best_score, best_idx] = max(scores);

    % Check if we have valid scores
    if best_score == -inf
        break;  % No valid conversions possible
    end

    % Accept if improvement (for first conversion, any positive SINR is good)
    if iter == 1
        % First conversion: accept if SINR > 0
        if best_score > 0
            alpha(tx_indices(best_idx)) = 0;
            current_SINR = best_score;
        else
            break;
        end
    else
        % Subsequent conversions: require improvement
        if best_score > current_SINR * 1.001
            alpha(tx_indices(best_idx)) = 0;
            current_SINR = best_score;
        else
            break;  % No significant improvement
        end
    end
end

% Ensure minimum diversity (at least 2 Rx for better performance)
num_rx = sum(~alpha);
num_tx = sum(alpha);

if num_rx < 2 && num_tx > min_tx
    tx_indices = find(alpha == 1);
    % Convert the last Tx to Rx
    alpha(tx_indices(end)) = 0;
end

end
