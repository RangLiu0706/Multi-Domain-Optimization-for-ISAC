function He_nk = compute_He(Hjnk, pk, Pj, alpha, beta_tilde)
[K, ~, M2, Ns, J] = size(Hjnk); M = M2/2;
He_nk = zeros(K, J*M, Ns);
for j = 1:J
    for n = 1:Ns
        for k = 1:K
            He_nk(k, (j-1)*M+1:j*M, n) = beta_tilde(n) * alpha(j) * pk(:,k).' * ...
                squeeze(Hjnk(k, :, :, n, j)) * Pj(:,:,j);
        end
    end
end
end
