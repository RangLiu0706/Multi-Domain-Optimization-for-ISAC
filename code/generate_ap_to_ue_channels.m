function Hjnk = generate_ap_to_ue_channels(ap_pos, ue_pos, params)
% Generate AP-to-UE communication channels
J = params.J; K = params.K; M = params.M; Ns = params.Ns;

Hjnk = complex(zeros(K, 2, 2*M, Ns, J));
for j = 1:J
    for k = 1:K
        d_jk = norm(ap_pos(j,:) - ue_pos(k,:));
        PL_dB = params.FSPL_1m + 10*params.nPLE_NLOS*log10(d_jk) + ...
            params.sigma_SF_UE * randn();
        gain = 10^(-PL_dB/20);
        for n = 1:Ns
            W = (randn(2*M, 2) + 1j*randn(2*M, 2)) / sqrt(2);
            Hjnk(k,:,:,n,j) = (gain * params.Rsqrt * W * params.Rsqrt_rx.').';
        end
    end
end
end