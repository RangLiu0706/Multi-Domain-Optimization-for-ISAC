function pk = update_pk(Hjnk, Pj, W_jn, t_nk, c_nk)
[K, ~, M2, Ns, J] = size(Hjnk); M = M2/2;
pk = zeros(2, K);
for k = 1:K
    ak = zeros(2, 1); Bk = zeros(2, 2);
    for n = 1:Ns
        akj = zeros(2, 1);
        for j = 1:J
            akj = akj + squeeze(Hjnk(k, :, :, n, j)) * Pj(:,:,j) * W_jn(M*(j-1)+1:j*M, k, n);
        end

        ak = ak + 2*sqrt(1 + t_nk(n,k)) * c_nk(n,k)' * akj;

        for kt = 1:K
            bnkkt = zeros(2, 1);
            for j = 1:J
                bnkkt = bnkkt + abs(c_nk(n,k)) * squeeze(Hjnk(k, :, :, n, j)) * ...
                    Pj(:,:,j) * W_jn(M*(j-1)+1:j*M, kt, n);
            end
            Bk = Bk + bnkkt * bnkkt';
        end
    end

    a = real(ak); B = real(Bk);
    a1 = a(1); a2 = a(2);
    b1 = B(1,1); b2 = B(1,2); b4 = B(2,2);

    f = @(theta) -(a1*cos(theta) + a2*sin(theta) - ...
        b1*(cos(theta))^2 - 2*b2*cos(theta)*sin(theta) - b4*(sin(theta))^2);

    [theta_opt, ~] = fminbnd(f, -pi/2, pi/2);
    pk(:,k) = [cos(theta_opt); sin(theta_opt)];
end
end
