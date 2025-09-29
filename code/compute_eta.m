function eta_in = compute_eta(theta_AP, XPD, Wr_jn, Phiij, Pj, alpha, beta, u_in)
[J, Ns] = size(beta); M = size(Pj, 2); idxRx = find(~alpha);
eta_in = zeros(M, J, Ns);
for irx = idxRx.'
    Ar = kron(exp(1j*pi*sin(-theta_AP(irx))*(0:M-1).'), XPD);
    for n = 1:Ns
        jtx = find(beta(:,n));
        if isempty(jtx) || alpha(irx) == 1
            continue;
        end
        Atj = kron(exp(1j*pi*sin(-theta_AP(jtx))*(0:M-1).'), XPD);
        Temp = Phiij(:,:,J*(irx-1)+jtx) * Atj.' * Pj(:,:,jtx) * alpha(jtx) * ...
            beta(jtx,n) * Wr_jn(M*(jtx-1)+1:jtx*M, :, n);
        eta_in(:,irx,n) = (u_in(:,irx,n)' * Pj(:,:,irx).' * Ar * Temp)';
    end
end
end
