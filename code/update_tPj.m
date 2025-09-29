function Pj = update_tPj(Hjnk, pk, W_jn, t_nk, c_nk, alpha, beta, Pj, Rc, sigma2_c)
[K, ~, M2, Ns, J] = size(Hjnk); M = M2/2;
p = zeros(2*J*M, 1);
for j = 1:J
    for m = 1:M
        p(2*M*(j-1)+2*(m-1)+1:2*M*(j-1)+2*m, 1) = Pj(2*(m-1)+1:2*m, m, j);
    end
end

const = Rc*log(2) + sigma2_c*sum(abs(c_nk(:)).^2) + sum(t_nk(:)) - sum(log(1+t_nk(:)));
idxTx = find(alpha);
f = zeros(2*J*M, 1);
for j = idxTx.'
    aj = zeros(2*M, 1);
    for n = 1:Ns
        if sum(beta(:,n)) == 1
            continue;
        end
        for k = 1:K
            t1 = 2*sqrt(1+t_nk(n,k)) * c_nk(n,k)' * (squeeze(Hjnk(k,:,:,n,j))).' * pk(:,k);
            aj = aj + t1 .* kron(W_jn(M*(j-1)+1:j*M, k, n), [1;1]);
        end
    end
    f(2*M*(j-1)+1:2*j*M, 1) = aj;
end
f = real(f);
F = zeros(2*J*M, 2*J*M);
for n = 1:Ns
    if sum(beta(:,n)) == 1
        continue;
    end
    for k = 1:K
        for kt = 1:K
            b_nkkt = zeros(2*J*M, 1);
            for j = idxTx.'
                t2 = abs(c_nk(n,k)) * (squeeze(Hjnk(k,:,:,n,j))).' * pk(:,k);
                b_nkkt(2*M*(j-1)+1:2*j*M, 1) = t2 .* kron(W_jn(M*(j-1)+1:j*M, kt, n), [1;1]);
            end
            F = F + b_nkkt * b_nkkt';
        end
    end
end
F = real(F+F')/2;
f = f/abs(const); F = F/abs(const);

mask2 = logical(kron(alpha, ones(2*M,1)));
f2 = f(mask2); F2 = F(mask2, mask2); p2 = p(mask2);
temp = 2*p(~mask2).' * F(~mask2, mask2);
f2 = f2 - temp.';
Mt2 = size(f2, 1); Mt = Mt2/2;

manifold = obliquefactory(2, Mt);
problem.M = manifold;
problem.cost = @(P) (P(:)).' * F2 * P(:) - real(f2.' * P(:));
problem.egrad = @(P) 2*reshape(F2*P(:), 2, Mt) - reshape(f2, 2, Mt);
options.verbosity = 0;

[P_opt, ~, ~] = conjugategradient(problem, reshape(p2, 2, Mt), options);
p2 = P_opt(:); p(mask2) = p2;
for j = 1:J
    if alpha(j) == 1
        for m = 1:M
            Pj(2*(m-1)+1:2*m, m, j) = p(2*M*(j-1)+2*(m-1)+1:2*M*(j-1)+2*m, 1);
        end
    end
end
end
