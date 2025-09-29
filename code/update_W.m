function Wjn = update_W(He_nk, t_nk, c_nk, Ptot, eta_in, gamma_in, alpha, beta, Rc, sigma2_c, u_in, Pj, Phiij, theta_AP, XPD, Gijn, Nsub)
[K, Mt, Ns] = size(He_nk); J = size(alpha, 1); M = Mt/J; Ng = Ns/Nsub;
beta_tilde = 1 - sum(beta, 1);
const = Rc*log(2) + sigma2_c*sum(abs(c_nk(:)).^2) + sum(t_nk(:)) - sum(log(1+t_nk(:)));

H1 = zeros(K, Mt, Nsub); H2 = zeros(Mt, Mt, K, Nsub);
for g = 1:Nsub
    idx = (g-1)*Ng+1:g*Ng;
    for k = 1:K
        Temp = zeros(Mt, Mt);
        for n = idx
            H1(k,:,g) = H1(k,:,g) + 2*sqrt(1+t_nk(n,k)) * c_nk(n,k)' * He_nk(k,:,n);
            Temp = Temp + abs(c_nk(n,k))^2 * He_nk(k,:,n)' * He_nk(k,:,n);
        end
        Temp = (Temp + Temp')/2;
        [U, S] = svd(Temp, 'econ');
        H2(:,:,k,g) = sqrt(S) * U';
    end
end

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
    [U, S, ~] = svd(Temp);
    % Only keep significant singular values
    tol = 1e-10 * S(1,1);
    num_sig = sum(diag(S) > tol);
    num_sig = min(num_sig, numel(idxRx));  % At most numel(idxRx) meaningful dimensions

    if num_sig > 0
        Bn(:,1:num_sig,g) = U(:,1:num_sig) * sqrt(S(1:num_sig,1:num_sig));
    end
end

scale = max(abs(Bn(:)));
if scale > max(abs(An(:)))/1e3
    An = An/(scale^2); Bn = Bn/scale;
else
    scale = max(abs(An(:)));
    An = An/scale; Bn = zeros(M, M, Nsub);
end

cvx_begin quiet
variable Wg(Mt, M+K, Nsub) complex
expression cqos
expression pow
expression ff
cqos = 0; pow = 0; ff = 0;
for g = 1:Nsub
    if beta_tilde(g*Ng) == 1
        for k = 1:K
            cqos = cqos + real(H1(k,:,g)*Wg(:,k,g)) - ...
                square_pos(norm(H2(:,:,k,g)*Wg(:,1:K,g), 'fro'));
        end
        pow = pow + square_pos(norm(Wg(:,1:K,g), 'fro'));
    else
        pow = pow + square_pos(norm(Wg((jtx(g)-1)*M+1:jtx(g)*M, K+1:end, g), 'fro'));
        ff = ff + real(trace(An(:,:,g)*Wg((jtx(g)-1)*M+1:jtx(g)*M, K+1:end, g))) - ...
            square_pos(norm(Bn(:,:,g)'*Wg((jtx(g)-1)*M+1:jtx(g)*M, K+1:end, g), 'fro'));
    end
end
maximize ff
subject to
cqos >= const;
pow <= Ptot/Ng;
cvx_end

Wjn = zeros(Mt, K+M, Ns);
for g = 1:Nsub
    idx = (g-1)*Ng+1:g*Ng;
    Wjn(:,:,idx) = repmat(Wg(:,:,g), 1, 1, Ng);
end
end
