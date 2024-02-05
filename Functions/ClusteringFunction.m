function  T = ClusteringFunction(clustering_data, metric, state_of_interest, diet_of_interest, nr_clusters)

n_inds = size(clustering_data,1);
n_states = size(clustering_data,2);
n_time = size(clustering_data,3);
n_diets = size(clustering_data,4);
similarity_matrix = zeros(n_inds,n_inds);
    %% State trajectories
    for i = 1:n_inds
        for j = 1:n_inds
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

end