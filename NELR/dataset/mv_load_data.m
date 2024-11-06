function [X,gt] = mv_load_data(dataname)
    load(dataname);
    X = full(data);
    if iscell(truelabel)
        gt = double(truelabel{1});
    else
        gt = double(truelabel);
    end
    if size(gt,1) == 1
        gt = gt';
    end
end

