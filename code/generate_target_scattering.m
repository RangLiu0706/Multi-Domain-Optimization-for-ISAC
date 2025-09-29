function Phiij = generate_target_scattering(ap_pos, target_pos, params)
% Generate polarization-dependent target scattering matrices
J = params.J; XPD_main = 10^(20/10);

Phiij = complex(zeros(2, 2, J, J));
for i = 1:J
    d_i = norm(ap_pos(i,:) - target_pos);
    g_i = 10^(-(params.FSPL_1m + 20*log10(d_i))/20);
    v_i = [1; 1/XPD_main * exp(1j*2*pi*rand)];
    for j = 1:J
        d_j = norm(ap_pos(j,:) - target_pos);
        g_j = 10^(-(params.FSPL_1m + 20*log10(d_j))/20);
        v_j = [1; 1/XPD_main * exp(1j*2*pi*rand)];

        vec_i = ap_pos(i,:) - target_pos;
        vec_j = ap_pos(j,:) - target_pos;
        beta = acosd(dot(vec_i, vec_j) / (d_i * d_j));

        alpha_ij = sqrt(params.sigma_b2) * cosd(beta/2)^params.nPat * ...
            raylrnd(1) * exp(1j*2*pi*rand);

        Phiij(:,:,i,j) = g_i * g_j * alpha_ij * (v_i * v_j.');
    end
end
end