function t_nk = compute_t(He, Wjn, sigma2, beta_tilde)
[K, ~, Ns] = size(He);
t_nk = zeros(Ns, K);
for n = 1:Ns
    if beta_tilde(n) == 1
        for k = 1:K
            t_nk(n,k) = abs(He(k,:,n) * Wjn(:,k,n))^2 / ...
                (norm(He(k,:,n) * Wjn(:,:,n))^2 - abs(He(k,:,n) * Wjn(:,k,n))^2 + sigma2);
        end
    end
end
end
