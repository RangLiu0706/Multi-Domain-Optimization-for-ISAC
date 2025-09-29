function [u_in, gamma_in] = update_u(Pj, theta_AP, Phiij, alpha, beta, Wr_jn, Gijn, sigma2_r, XPD)
[J, Ns] = size(beta); M = size(Pj, 2); idxRx = find(~alpha);
u_in = zeros(M, J, Ns); gamma_in = zeros(J, Ns);
for irx = idxRx.'
    Ar = kron(exp(1j*pi*sin(-theta_AP(irx))*(0:M-1).'), XPD);
    for n = 1:Ns
        jtx = find(beta(:,n));
        if isempty(jtx) || alpha(irx) == 1
            continue;
        end

        Bin = zeros(M);  Cin = zeros(M);
        maskTxRows = (jtx-1)*M+1:jtx*M;
        Wr_n_j = Wr_jn(maskTxRows, :, n);

        At = kron(exp(1j*pi*sin(-theta_AP(jtx))*(0:M-1).'), XPD);
        Phi = Phiij(:,:,(irx-1)*J + jtx);

        Bin = Bin + Pj(:,:,irx).' * Ar * Phi * At.' * Pj(:,:,jtx) * Wr_n_j;
        Cin = Cin + Pj(:,:,irx).' * Gijn(:,:,n,irx,jtx) * Pj(:,:,jtx) * Wr_n_j;

        Bin = Bin / sqrt(sigma2_r);
        Cin = Cin / sqrt(sigma2_r);

        R_i = Cin * Cin' + eye(M);
        [V, D] = eig((R_i\Bin) * (Bin'/R_i));
        [~, idx] = max(real(diag(D)));
        u_opt = V(:,idx); u_opt = u_opt / norm(u_opt);

        sig = norm(u_opt' * Bin)^2;
        intf = norm(u_opt' * Cin)^2;
        gamma = sig / (intf + 1);

        u_in(:, irx, n) = u_opt;
        gamma_in(irx, n) = gamma;
    end
end
end
