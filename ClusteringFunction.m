function  T = ClusteringFunction(clustering_data, metric, state_of_interest, diet_of_interest, nr_clusters)
similarity_matrix = zeros(17,17);

    %% State trajectories
    for i = 1:17
        for j = 1:17
            A_1 = squeeze( clustering_data(i,state_of_interest,:,diet_of_interest) );
            A_2 = squeeze( clustering_data(j,state_of_interest,:,diet_of_interest) );
            if metric == 1
                similarity_matrix(i,j) = real( A_1(:)'*A_2(:) )/ ( sqrt(A_1(:)'*A_1(:)) * sqrt(A_2(:)'*A_2(:)) );
            elseif metric ==2
                similarity_matrix(i,j) = sum( ( A_1(:) - A_2(:) ).^2, 'all');
            end
        end
    end

    if metric ==2
        cluster_matrix = squeeze( clustering_data(:,state_of_interest,:,diet_of_interest) );
    end

%% Clustering
if metric == 1
    Z = linkage(similarity_matrix,'complete','cosine');
elseif metric == 2
    Z = linkage(real( cluster_matrix ),'ward');
end
T = cluster(Z,'maxclust',nr_clusters);

%% Plot 
    figure()
    %subplot(121)
    colors = {'b','r','g','k'};
    lines = {'--','-'};
for i = unique(T)'
    plot(1:8, -1*real( squeeze( clustering_data(i==T,state_of_interest,:,diet_of_interest)) ),strcat( lines{i}, colors{i}), 'LineWidth', 2); hold on
end
title('\textbf{State trajectory clustering}', 'Interpreter','latex'); 
xlabel('\textbf{Time points}', 'Interpreter','latex')
ylabel('\textbf{State value}', 'Interpreter','latex')

figure()


cutoff = median([Z(end,3) Z(end,3)]);
h = dendrogram(Z,'ColorThreshold',cutoff);
set(h,'LineWidth',2)


linesColor = cell2mat(get(h,'Color')); % get lines color; 
colorList = unique(linesColor, 'rows');
title('\textbf{Agglomerative clustering of similarity matrix}', 'Interpreter','latex'); 

end