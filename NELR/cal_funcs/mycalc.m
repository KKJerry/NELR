function [acc, nmi, purity, fmeasure, ri, ari] = mycalc(groundLables, currentLabels)
        numClusters = length(unique(groundLables));
        acc = accuracy(groundLables, currentLabels);
%         [sortedLabels] = bestMap(groundLables, currentLabels);
%         acc1 = mean(groundLables==sortedLabels);
                
        class_labels = zeros(1, numClusters);
        for idx =  1 : numClusters
            class_labels(idx) = length(find(groundLables == idx));
        end
        cluster_data = cell(1, numClusters);
        for idx =  1 : numClusters
            cluster_data(1, idx) = { groundLables(currentLabels == idx)' };
        end
        [nmi, purity, ~, ri, ari] = calculate_results(class_labels, cluster_data);

        fmeasure = compute_f(groundLables, currentLabels);
        nmi = nmi*100;
        purity = purity * 100;
        fmeasure = fmeasure * 100;
        ri = ri * 100;
        ari = ari * 100;
end