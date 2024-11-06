function [NewSamOut, recordW] = K3_Neig(oldSamMat, threK, max_ly)
    n = size(oldSamMat,1);
    [~,idx] = sort(oldSamMat,2,'descend');
    k_ids = zeros(n, n);
    neigNums = sum(oldSamMat > 0)';
    for i = 1:n
        % k_ids是一个n*k的矩阵，对于每个样本，存放着最近的K个节点的编号
        if neigNums(i) > 1
            k_ids(i, 1:neigNums(i) - 1) = idx(i,2:neigNums(i));
        end
        neigNums(i) = neigNums(i) - 1;  % 自身不算邻居，放这里减去
    end

    k_ids_next_t = k_ids;
    neigNums_next_t = neigNums;
    NewSamOut = oldSamMat;
    recordW{1} = NewSamOut;

    for ll = 1:max_ly
        neigNums_t = neigNums_next_t;
        k_ids_t = k_ids_next_t;
        kSamMat_t = NewSamOut;
        recordW{ll} = NewSamOut;
        cov_end = true;
        for i = 1:n
            for j=1:n
                if (kSamMat_t(i,j) == 0)  % 对于本来不是邻居的节点，检查他们有没有相交的点
                    [H,m,d]=intersect(k_ids_t(i,1:neigNums_t(i)),k_ids_t(j,1:neigNums_t(j)));
                    Hn = length(H); % 邻域重复点个数 
                    if (Hn >= threK)
                        NewSamOut(i,j) = (mean(kSamMat_t(i,k_ids_t(i,m)) + kSamMat_t(j,k_ids_t(j,d))) ./ 2) * 1;
%                          fprintf("扩展节点%d->%d %f\n", i, j, NewSamOut(i,j));
%                          disp(H);
                        k_ids_next_t(i,neigNums_next_t(i)+1) = j;
                        neigNums_next_t(i) = neigNums_next_t(i)+1;
                        cov_end = false;
                    end
                end
            end
        end
        if cov_end
            disp('break');
            break;
        end
    end





