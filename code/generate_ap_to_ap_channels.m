function Gijn = generate_ap_to_ap_channels(ap_pos, params)
% Generate inter-AP interference channels
J = params.J; M = params.M; Ns = params.Ns;

Gijn = complex(zeros(2*M, 2*M, Ns, J, J));
for i = 1:J
    for j = [1:i-1, i+1:J]   
        d_ij = norm(ap_pos(i,:) - ap_pos(j,:));
        PL_dB = params.FSPL_1m + 10*params.nPLE_NLOS*log10(d_ij) + ...
            params.sigma_SF_AP * randn();
        gain = 10^(-PL_dB/20);
        for n = 1:Ns
            W = (randn(2*M) + 1j*randn(2*M)) / sqrt(2);
            Gijn(:,:,n,i,j) = gain * params.Rsqrt * W * params.Rsqrt.';
        end
    end
end
end