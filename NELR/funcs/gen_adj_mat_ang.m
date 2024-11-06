function [A] = gen_adj_mat_ang(Zn, alpha)
    nv = max(size(Zn));
    n = size(Zn{1},1);
    Zs = zeros(n, n);
    for idx = 1 : nv 
        Zs = Zs + (Zn{idx}); 
    end
    [U, s, V] = svd(Zs);
    s = diag(s);
    r = sum(s>1e-6);
    U = U(:, 1 : r);
    s = diag(s(1 : r));
    V = V(:, 1 : r);
    M = U * s.^(1/2);
    mm = normr(M);
    rs = mm * mm';
    A = rs.^(2 * alpha);
end

