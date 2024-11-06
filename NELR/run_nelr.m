clear;
clc;

originalPath = path;
cleanupObj = onCleanup(@() path(originalPath));

addpath('./cal_funcs')
addpath('./funcs')
addpath('./dataset')

resultdir = 'result/';
if(~exist('result','file'))
    mkdir('result');
end

%   3sources  Wiki  BBC  BBCSport  movie  100leaves
%% data
dtname = '3sources';

[X, gt] = mv_load_data(dtname);
[X, gt] = order_data(X, gt');
gt = gt';
n = size(X{1}, 2);
nv = size(X, 2);
K = length(unique(gt));
for nv_idx = 1 : nv 
    if ~strcmp(dtname, 'caltech') && ~strcmp(dtname, 'xx')
        disp('norm data')
        for idx = 1 : n
            sample = X{nv_idx}(:, idx);
            X{nv_idx}(:, idx) = sample ./ max(1e-12, norm(sample));
        end
    end

    if strcmp(dtname, '3sources')
        new_dim = 160;  % 3sources
        [eigen_vector, ~] = f_pca(X{nv_idx}, new_dim);
        X{nv_idx} = eigen_vector' *  X{nv_idx};
    end
    
    if strcmp(dtname, 'BBC')
        new_dim = 600;  % BBC
        [eigen_vector, ~] = f_pca(X{nv_idx}, new_dim);
        X{nv_idx} = eigen_vector' *  X{nv_idx};
    end



    if strcmp(dtname, 'BBCSport')
        new_dim = 500;  % BBCSport
        [eigen_vector, ~] = f_pca(X{nv_idx}, new_dim);
        X{nv_idx} = eigen_vector' *  X{nv_idx};
    end


    if strcmp(dtname, 'movie')
        new_dim = 600;  % movie
        [eigen_vector, ~] = f_pca(X{nv_idx}, new_dim);
        X{nv_idx} = eigen_vector' *  X{nv_idx};
    end

end



repNum = 10;
filePrefix = 'nelr_';
resultname = [resultdir, filePrefix, char(dtname),'_result_',datestr(clock,'yyyy_mm_dd_HH_MM_SS'),'_.mat'];

%% params


if strcmp(dtname, 'Wiki')
    nei_num = 6; lambda = 1; eta = 0.1; alpha = 5; thK = nei_num / 2; lyMAX = 5; tol = 1e-4;%Wiki inf
elseif strcmp(dtname, '3sources')
    nei_num = 6; lambda = 0.5; eta = 0.1; alpha = 5; thK = nei_num / 2; lyMAX = 5; tol = 1e-4; %3sources inf
elseif strcmp(dtname, 'BBCSport')
    nei_num = 4; lambda = 1; eta = 0.1; alpha = 4; thK = nei_num / 2; lyMAX = 5;  tol = 1e-4;%BBCSport inf
elseif strcmp(dtname, 'BBC')
    nei_num = 4; lambda = 0.05; eta = 0.1; alpha = 1; thK = nei_num / 2; lyMAX = 5;  tol = 1e-4; %BBC inf
elseif strcmp(dtname, 'movie')
    nei_num = 1; lambda = 0.1; eta = 0.1; alpha = 1; thK = nei_num / 2; lyMAX = 5;  tol = 1e-4; %movie inf
elseif strcmp(dtname, '100leaves')
    nei_num = 7; lambda = 0.01; eta = 5000; alpha = 1; thK = nei_num / 2; lyMAX = 5; tol = 1; %100leaves inf
end



An = cell(1, nv);
Aout = cell(1, nv);
normA = cell(1, nv);
eyeA = cell(1, nv);
parfor nv_idx = 1 : nv     
    An{nv_idx} = getSamMatByDis(X{nv_idx}, nei_num);
    An{nv_idx} = removeOneWayLink(An{nv_idx}); 
    An{nv_idx} = (An{nv_idx} + An{nv_idx}') ./ 2;   
    [Aout{nv_idx},recordA{nv_idx}] = K3_Neig(An{nv_idx}, thK, lyMAX);
    D = diag(sum(Aout{nv_idx}, 2) .^ -0.5);
    normA{nv_idx} = D * Aout{nv_idx} * D;
    eyeA{nv_idx} = eye(n);
end




[Zn] = nelr(X, normA, lambda, eta, tol);
Z = gen_adj_mat_ang(Zn, alpha);
for rep_idx = 1:repNum
    pre_labels = new_spectral_clustering(Z, K);
    [acc, nmi, purity, fmeasure, ri, ari] = mycalc(gt, pre_labels);
    result_tmp(rep_idx,:)  = [acc, nmi, purity, fmeasure, ri, ari];
end
result = [lambda, eta, alpha, nei_num, thK, lyMAX, mean(result_tmp),std(result_tmp), tol];
save(resultname , 'result');
disp([lambda, eta, alpha, nei_num, tol]);
% disp(result_tmp);
disp('ACC        NMI        PURITY        FMS        RI        ARI')
disp([mean(result_tmp)]);
disp([std(result_tmp)]);






