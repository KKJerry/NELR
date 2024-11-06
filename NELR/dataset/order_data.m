% def order_sam_for_diag(x, y):
%     """
%     rearrange samples
%     :param x: feature sets
%     :param y: ground truth
%     :return:
%     """
%     x_new = np.zeros(x.shape)
%     y_new = np.zeros(y.shape)
%     start = 0
%     for i in np.unique(y):
%         idx = np.nonzero(y == i)
%         stop = start + idx[0].shape[0]
%         x_new[start:stop] = x[idx]
%         y_new[start:stop] = y[idx]
%         start = stop
%     return x_new, y_new


function [newX, newGt] = order_data(X, gt)
    nv = size(X, 2);
    newX = cell(1,nv);
    newGt = zeros(size(gt));
    for v = 1 : nv
        start = 1;
        newX{v} = zeros(size(X{v}));
        for i = unique(gt)
            idx = gt == i;
            stop = start + sum(idx) - 1;
            newX{v}(:,start:stop) = X{v}(:,idx);
            newGt(start:stop) = gt(idx);
            start = stop + 1;
        end
    end
end

