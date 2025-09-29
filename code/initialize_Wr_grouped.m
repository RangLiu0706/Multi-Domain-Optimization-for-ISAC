function [Wjn, SINR] = initialize_Wr_grouped(theta_AP, Gijn, Phiij, P, sigma2_r, XPD, Pj, alpha, beta, Nsub)
% Sensing beamformer initialization with subcarrier grouping

% ================ System Parameters ================
[M2, ~, Ns, J, ~] = size(Gijn);
M = M2/2;
Ng = Ns / Nsub;

idx_Rx = find(~alpha);
Nr = numel(idx_Rx) * M;

% Find sensing groups
idxSenseGroup = [];
for g = 1:Nsub
    idx = (g-1)*Ng+1:g*Ng;
    if any(any(beta(:,idx)))
        idxSenseGroup = [idxSenseGroup, g];
    end
end
Nsub_sense = numel(idxSenseGroup);

% ================ Step 1: Design Group Beamformers ================
g = zeros(Nsub_sense, 1);
wUnit = cell(Nsub_sense, 1);
HtStore = cell(Nsub_sense, 1);
HIStore = cell(Nsub_sense, 1);
maskTxStore = cell(Nsub_sense, 1);

for gg = 1:Nsub_sense
    group = idxSenseGroup(gg);
    idx = (group-1)*Ng+1:group*Ng;

    % Find transmitting BS for this group
    tx_count = sum(beta(:,idx), 2);
    [~, j_tx] = max(tx_count);

    maskTx = false(1, J*M);
    maskTx((j_tx-1)*M+1:j_tx*M) = true;
    Nt = M;

    % Aggregate channels for the group
    Ht_group = zeros(Nr, Nt);
    HI_group = zeros(Nr, Nt);

    active_count = 0;
    for n = idx
        if beta(j_tx, n)
            active_count = active_count + 1;

            Ht_n = zeros(Nr, Nt);
            HI_n = zeros(Nr, Nt);

            for ir = 1:numel(idx_Rx)
                r = idx_Rx(ir);
                rows = (ir-1)*M+1:ir*M;

                HI_n(rows, :) = Pj(:,:,r).' * Gijn(:,:,n,r,j_tx) * Pj(:,:,j_tx);

                Ar = kron(exp(1j*pi*sin(-theta_AP(r))*(0:M-1).'), XPD);
                At = kron(exp(1j*pi*sin(-theta_AP(j_tx))*(0:M-1).'), XPD);
                Phi = Phiij(:,:,(r-1)*J + j_tx);
                Ht_n(rows, :) = Pj(:,:,r).' * Ar * Phi * At.' * Pj(:,:,j_tx);
            end

            Ht_group = Ht_group + Ht_n;
            HI_group = HI_group + HI_n;
        end
    end

    if active_count > 0
        Ht_group = Ht_group / active_count;
        HI_group = HI_group / active_count;
    end

    HtStore{gg} = Ht_group;
    HIStore{gg} = HI_group;
    maskTxStore{gg} = maskTx;

    rankHI = rank(HI_group, 1e-8);
    if rankHI < Nt
        % Null space method with error checking
        V0 = null(HI_group, 'r');

        if isempty(V0)
            % No null space found, fall back to regularized solution
            P_eq = P / Nsub_sense;
            alpha_rzf = Nt * sigma2_r / P_eq;
            [~, ~, Vh] = svd(Ht_group, 'econ');
            h_eff = Vh(:,1);
            w_raw = (HI_group'*HI_group + alpha_rzf*eye(Nt)) \ h_eff;
        else
            % Project target channel to null space
            Ht_proj = Ht_group * V0;

            % Check if projection is non-zero
            if norm(Ht_proj, 'fro') < 1e-10
                % Target channel orthogonal to null space, use regularized solution
                P_eq = P / Nsub_sense;
                alpha_rzf = Nt * sigma2_r / P_eq;
                [~, ~, Vh] = svd(Ht_group, 'econ');
                if isempty(Vh)
                    % Fallback: use random direction
                    w_raw = randn(Nt, 1) + 1j*randn(Nt, 1);
                else
                    h_eff = Vh(:,1);
                    w_raw = (HI_group'*HI_group + alpha_rzf*eye(Nt)) \ h_eff;
                end
            else
                % Normal null space processing
                [~, ~, Vh] = svd(Ht_proj, 'econ');
                if isempty(Vh) || size(Vh,2) < 1
                    % SVD failed, use first null space vector
                    w_raw = V0(:,1);
                else
                    w_raw = V0 * Vh(:,1);
                end
            end
        end
    else
        % Full rank case - regularized solution
        P_eq = P / Nsub_sense;
        alpha_rzf = Nt * sigma2_r / P_eq;

        % Check if Ht_group is valid
        if norm(Ht_group, 'fro') < 1e-10
            % Degenerate case: no target signal
            w_raw = randn(Nt, 1) + 1j*randn(Nt, 1);
        else
            [~, ~, Vh] = svd(Ht_group, 'econ');
            if isempty(Vh)
                w_raw = randn(Nt, 1) + 1j*randn(Nt, 1);
            else
                h_eff = Vh(:,1);
                w_raw = (HI_group'*HI_group + alpha_rzf*eye(Nt)) \ h_eff;
            end
        end
    end

    % Normalize beamformer
    if norm(w_raw) < 1e-10
        w_raw = randn(Nt, 1) + 1j*randn(Nt, 1);  % Random initialization if all else fails
    end
    wUnit{gg} = w_raw / norm(w_raw);

    % Compute channel quality with numerical stability
    signal_power = norm(Ht_group * wUnit{gg})^2;
    interference_power = norm(HI_group * wUnit{gg})^2;
    g(gg) = signal_power / (interference_power + sigma2_r);
end

% ================ Step 2: Power Allocation ================
if sum(g) > 0
    pAlloc = P * g ./ sum(g);
else
    % Equal allocation if all channels are bad
    pAlloc = P * ones(Nsub_sense, 1) / Nsub_sense;
end

% ================ Step 3: Apply to All Subcarriers ================
Wjn = zeros(J*M, M, Ns);
SINR = 0;

for gg = 1:Nsub_sense
    group = idxSenseGroup(gg);
    idx = (group-1)*Ng+1:group*Ng;

    if pAlloc(gg) > 0
        w_group = sqrt(pAlloc(gg)/Ng) * wUnit{gg};

        % Apply same beamformer to all subcarriers in group
        for n = idx
            Wtmp = zeros(J*M, 1);
            Wtmp(maskTxStore{gg}) = w_group;
            Wjn(:,1,n) = Wtmp;
        end

        % Compute SINR for the group
        Ht = HtStore{gg};
        HI = HIStore{gg};

        for ir = 1:numel(idx_Rx)
            rows = M*(ir-1)+1:ir*M;
            HI_r = HI(rows,:);
            Ht_r = Ht(rows,:);

            % Check for valid channels
            if norm(HI_r, 'fro') > 0 || norm(Ht_r, 'fro') > 0
                R = HI_r * (w_group*w_group') * HI_r' + sigma2_r*eye(M);
                if cond(R) > 1e12
                    R = R + 1e-6 * trace(R) / M * eye(M);
                end

                u = R \ (Ht_r * w_group);
                if norm(u) > 0
                    SINR_r = real(w_group' * Ht_r' * u) * Ng;
                    SINR = SINR + SINR_r;
                end
            end
        end
    end
end
end
