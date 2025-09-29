function beta_sub = allocate_subcarrier(Gijn, alpha, Nsub)
% Simple 50-50 allocation with load balancing for ISAC system

[~, ~, Ns, J, ~] = size(Gijn); Ng = Ns/Nsub;

idxTx = find(alpha == 1); idxRx = find(alpha == 0);
n_tx = length(idxTx); n_rx = length(idxRx);

% Initialize
beta_sub = zeros(J, Nsub);

% Handle edge cases
if isempty(idxRx) || isempty(idxTx)
    return;  % All blocks for communication
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

% Sort blocks by average interference
block_avg_interf = mean(InterferenceTable, 1);
[~, sorted_blocks] = sort(block_avg_interf);

% 50-50 split: low interference blocks for radar
num_radar_blocks = floor(Nsub / 2);
radar_blocks = sorted_blocks(1:num_radar_blocks);

% Round-robin allocation (most uniform)
for idx = 1:num_radar_blocks
    b = radar_blocks(idx);
    tx_idx = mod(idx-1, n_tx) + 1;  % Cycle through Tx APs
    beta_sub(idxTx(tx_idx), b) = 1;
end

end