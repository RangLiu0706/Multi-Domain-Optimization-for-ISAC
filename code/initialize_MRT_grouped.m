function [Wjn, Pmin] = initialize_MRT_grouped(He_nk, sigma2, alpha, beta_tilde, Rc, Nsub)
% MRT-based initial communication beamformer with subcarrier grouping
% Each group of Ng consecutive subcarriers uses the same beamformer
% ================ Parameter Extraction ================
[K, JM, Ns] = size(He_nk);
Ng = Ns / Nsub;
J = size(alpha, 1);
M = JM/J;
mask = logical(kron(alpha.', ones(1,M)));
He = He_nk(:, mask, :);
Mt = size(He, 2);
beta_tilde = logical(beta_tilde(:).');

% ================ Step 1: Compute Group-based MRT Beamformers ================
Wjn = zeros(JM, K, Ns);
Wg_unit = zeros(Mt, K, Nsub);  % One beamformer per group

for g = 1:Nsub
    idx = (g-1)*Ng+1:g*Ng;  % Subcarrier indices for this group

    % Check if this group is for communication
    if ~beta_tilde(g*Ng)  % Use last subcarrier of group as indicator
        continue;
    end

    % Aggregate channel for the group
    H_group = mean(He(:,:,idx), 3);

    % MRT beamformer for the group
    W_group = zeros(Mt, K);
    for k = 1:K
        w_k = H_group(k,:)';
        W_group(:,k) = w_k / norm(w_k);
    end

    Wg_unit(:,:,g) = W_group;

    % Apply same beamformer to all subcarriers in group
    for n = idx
        if beta_tilde(n)
            Wjn(mask, :, n) = W_group;
        end
    end
end

% ================ Step 2: Find Minimum Power via Binary Search ================
scale_min = 0;
scale_max = 100;
tolerance = 1e-6;
max_iter = 50;

for iter = 1:max_iter
    scale = (scale_min + scale_max) / 2;

    % Compute sum-rate with grouped beamformers
    sum_rate = compute_grouped_sum_rate(He, Wjn(mask,:,:), scale, beta_tilde, sigma2, Nsub);

    if sum_rate >= Rc
        scale_max = scale;
    else
        scale_min = scale;
    end

    if abs(scale_max - scale_min) < tolerance
        break;
    end
end

optimal_scale = scale_max;

% ================ Step 3: Apply Power ================
% Power calculation considering grouping
num_active_groups = sum(beta_tilde(1:Ng:Ns));  % Count active groups
Pmin = optimal_scale^2 * K * sum(beta_tilde);  % Total power
Wjn = optimal_scale * Wjn;
end

function sum_rate = compute_grouped_sum_rate(He, W, scale, beta_tilde, sigma2, Nsub)
[K, ~, Ns] = size(He);
Ng = Ns / Nsub;
sum_rate = 0;

for g = 1:Nsub
    idx = (g-1)*Ng+1:g*Ng;

    if ~beta_tilde(g*Ng)
        continue;
    end

    % Sum rate over all subcarriers in the group
    for n = idx
        if ~beta_tilde(n)
            continue;
        end

        H = He(:,:,n);
        Wn = scale * W(:,:,n);
        HW = H * Wn;

        for k = 1:K
            signal_power = abs(HW(k,k))^2;
            interference_power = sum(abs(HW(k,:)).^2) - signal_power;
            SINR = signal_power / (interference_power + sigma2);
            sum_rate = sum_rate + log2(1 + SINR);
        end
    end
end
end
