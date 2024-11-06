function [W] = getSamMatByDis(X, k)
    [dim, n] = size(X);
    X = X';
    distmat = squareform(pdist(X));
    distmat(distmat == 0) = 9999999999999;
    for i = 1:1:n
        distmat(i,i) = 0;
    end
    samMat = 1./(distmat+1);
    
    W = zeros(n,n);
    if k ~= 0
        [a, idx] = sort(samMat, 2, 'descend'); % sort each row
        for i = 1:n
            id = idx(i,1:k+1);
            W(i,id) = samMat(i, id);
        end
    else
        W = samMat;
    end
end


