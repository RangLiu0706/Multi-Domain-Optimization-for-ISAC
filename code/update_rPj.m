function Pj = update_rPj(Wr_jn, alpha, beta, eta_in, u_in, Pj, Phiij, gamma_in, Gijn, theta_AP, XPD)
[J, Ns] = size(beta); M = size(Wr_jn, 2);
a = zeros(2*J*M, 1); B = zeros(2*J*M, 2*J*M);
idxRx = find(~alpha); nrAP = numel(idxRx);
for i = idxRx.'
    ai = zeros(2*M, 1);  Bi = zeros(2*M, 2*M);
    Ari = kron(exp(1j*pi*sin(-theta_AP(i))*(0:M-1).'), XPD);
    for n = 1:Ns
        j = find(beta(:,n));
        if isempty(j)
            continue;
        end

        Atj = kron(exp(1j*pi*sin(-theta_AP(j))*(0:M-1).'), XPD);
        ain = Ari * Phiij(:,:,J*(i-1)+j) * Atj.' * Pj(:,:,j) * ...
            alpha(j) * beta(j,n) * Wr_jn((j-1)*M+1:j*M, :, n) * eta_in(:,i,n);
        Bin = Gijn(:,:,n,i,j) * Pj(:,:,j) * alpha(j) * beta(j,n) * ...
            Wr_jn((j-1)*M+1:j*M, :, n);

        ai = ai + 2*ain .* kron(conj(u_in(:,i,n)), [1;1]);
        Bi = Bi + gamma_in(i,n) * kron(diag(u_in(:,i,n)'), eye(2)) * ...
            Bin * Bin' * kron(diag(u_in(:,i,n)), eye(2));
    end

    a((i-1)*2*M+1:i*2*M, 1) = ai;
    B((i-1)*2*M+1:i*2*M, (i-1)*2*M+1:i*2*M) = Bi;
end

scale = max(abs(a));
if scale > 0
    a = real(a/scale); B = real(B/scale);
else
    a = real(a); B = real(B);
end

mask = logical(kron(alpha(:)<1, true(2*M,1)));
a = a(mask); B = B(mask, mask);

manifold = obliquefactory(2, nrAP*M);
problem.M = manifold;
problem.cost = @(P) reshape(P, [], 1)' * B * reshape(P, [], 1) - real(reshape(P, [], 1)' * a);
problem.egrad = @(P) 2 * reshape(B * reshape(P, [], 1), 2, nrAP*M) - reshape(a, 2, nrAP*M);
options.verbosity = 0;

Pini = zeros(2, nrAP*M);
for idxj = 1:nrAP
    j = idxRx(idxj);
    for m = 1:M
        Pini(:, (idxj-1)*M+m) = Pj((m-1)*2+1:m*2, m, j);
    end
end

[P_opt, ~, ~] = conjugategradient(problem, Pini, options);
p = P_opt(:);
for idxj = 1:nrAP
    j = idxRx(idxj);
    for m = 1:M
        Pj(2*(m-1)+1:2*m, m, j) = p((idxj-1)*2*M+2*(m-1)+1:(idxj-1)*2*M+2*m, 1);
    end
end
end
