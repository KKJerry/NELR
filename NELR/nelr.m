function [Zn, iter] = nelr(data_views, An, lambda, eta, tol)
    disp('running...');
    rho = 1.1;
    max_mu = 1e10;
    mu = 1e-2;
    maxIter = 200;

    disp(['tol',num2str(tol),'   lambda', num2str(lambda), '   eta', num2str(eta)])

    nv = length(data_views);
    Jn = cell(1, nv);
    Zn = cell(1, nv);
    Yn1 = cell(1, nv);
    Yn2 = cell(1, nv);
    Xtn = cell(1, nv);
    invXn = cell(1, nv);
    min_values = zeros(1, nv);
    min_1_values = zeros(1, nv);
    min_2_values = zeros(1, nv);
    
    % intialize
    for idx = 1 : nv
        [m, n] = size(data_views{idx});
        Jn{idx} = zeros(n, n);
        En{idx} = zeros(m, n);
        Zn{idx} = zeros(n, n);
        Yn1{idx} = zeros(m, n);
        Yn2{idx} = zeros(n, n);
        Xtn{idx} = data_views{idx}' * data_views{idx};
        invXn{idx} = inv(An{idx}' * Xtn{idx} * An{idx} + (1 + (eta / mu) * (nv - 1) * eye(n)));
    end

    iter = 0;
    while iter < maxIter  
        iter = iter + 1;
        for idx = 1 : nv
            data_view = data_views{idx};
            [m, n] = size(data_view);
            % update Jn{idx}
            temp = Zn{idx} + Yn2{idx}/mu;
            temp1 = (temp + temp') / 2;
            [U s V] = svd(temp1);
            s = diag(s);
            r = sum(s > 1/mu); 
            Jn{idx} = U(:, 1 : r) * diag(s(1 : r) - 1/mu) * V(:, 1 : r)';
            
            % update Zn{idx}
            Zn_sum = zeros(size(Zn{idx}));
            for i = 1 : nv
                 if i ~= idx
                     Zn_sum = Zn_sum + Zn{i};
                 end
            end
            Zn{idx} =  invXn{idx} * (An{idx}' * Xtn{idx} - An{idx}' * data_view' * En{idx} + Jn{idx} + (An{idx}' * data_view' * Yn1{idx} - Yn2{idx} + eta * Zn_sum)/mu);

             
            %  update En{idx}           
            xmaz = data_view - data_view * An{idx} * Zn{idx};
            tmp = xmaz + Yn1{idx}/mu;
            deta = lambda / mu;
            for i = 1 : n
                nw = norm(tmp(:, i));
                if deta < nw
                    En{idx}(:,i) = (nw - deta) * tmp(:, i) / nw;
                else
                    En{idx}(:,i)= zeros(length(tmp(:, i)),1);
                end
            end 
            leq1 = xmaz - En{idx};
            leq2 = Zn{idx} - Jn{idx};
            Yn1{idx} = Yn1{idx} + mu * leq1;
            Yn2{idx} = Yn2{idx} + mu * leq2;                       
        end
        mu = min(max_mu, mu * rho);

        for idx = 1 : nv                
            leq1 = data_views{idx} - data_views{idx} * An{idx} * Zn{idx} - En{idx};
            leq2 = Zn{idx} - Jn{idx};
            min_1_values(idx) = norm(leq1,"inf");
            min_2_values(idx) = norm(leq2,"inf");
            min_values(idx) = max(min_1_values(idx), min_2_values(idx));                      
            Zn_sum = 0;
            for i = 1 : nv
                 if i ~= idx
                     Zn_sum = Zn_sum + norm(Zn{idx} - Zn{i}, 'fro')^2;
                 end
            end
        end
        
        isConv = 1;
        for idx = 1 : nv 
            if max(min_values) >= tol 
                isConv = 0;
                break;
            end
        end


        if isConv == 1
            break;
        end
    end
end
