function c_nk = compute_c(He, Wjn, t_nk, sigma2)
[K, ~, Ns] = size(He);
c_nk = zeros(Ns, K);
for n = 1:Ns
    for k = 1:K
        c_nk(n,k) = sqrt(1 + t_nk(n,k)) * He(k,:,n) * Wjn(:,k,n) / ...
            (norm(He(k,:,n) * Wjn(:,:,n))^2 + sigma2);
    end
end
end
