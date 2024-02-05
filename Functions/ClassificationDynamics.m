function [basis, classification_dynamics] = ClassificationDynamics(data_tensor, nr_components)

data_matrix = normalize(squeeze( mean( squeeze( mean(data_tensor,1) ), 3) )');

%correlation_matrix = corr(data_matrix);
cosine_similarity_matrix = data_matrix'*data_matrix ./ (vecnorm(data_matrix)'*vecnorm(data_matrix));

Z = linkage(cosine_similarity_matrix,'complete','cosine');
nr_clusters = nr_components;
T = cluster(Z,'maxclust',nr_clusters);
nr_time = size(data_tensor,3);
nr_individual = size(data_tensor,1);
nr_metabolites = size(data_tensor,3);
nr_diets = size(data_tensor,4);
basis = zeros(nr_time,nr_clusters);
list_indices = 1:nr_clusters;%[1 4 2 3];%[4 6 5 3 2 1];
%iter = 0;
for i = list_indices
%     iter = iter + 1;
%     subplot(1,nr_components,iter)
%     plot(1:7,data_matrix(:,T==i),'k');hold on
      basis(:,i) = mean(data_matrix(:,T==i),2)';
%     plot(1:7,basis(:,i),'r','LineWidth',2)
%     xlim([1,7])
%     ylim([min(data_matrix,[],'all'), max(data_matrix,[],'all')])
end
% sgtitle('Dynamical classification')

%basis = basis ./vecnorm(basis,2);
basis_metabolites= zeros(nr_individual,nr_clusters,nr_metabolites,nr_diets);
for i = 1:nr_individual
    for d = 1:nr_diets
        X_tmp = squeeze( data_tensor(i,:,:,d) );
        %basis_metabolites(i,:,:,d) = basis' * X_tmp;
%         basis_metabolites(i,:,:,d) = basis \ X_tmp;
    end
end

% basis_metabolites_pooled_estimate = [];%squeeze( mean(basis_metabolites,[1,4]) )';


% dendrogram(Z)
classification_dynamics = T;
%classification_description_dynamics = {'One peak = 4h',  'Oscil. peak = {3h, 7h}, dip = 5h', 'Fast, peak = 3h','Fast, peak = 2h' ,'Oscil. dip = 4h, peak = 6h', 'Slow, peak > 8h'};
%classification_description_dynamics = {'Fast, peak = 2h', 'Slow, peak > 8h', 'Oscil. dip = 4h, peak = 6h','Fast, peak = 3h'};


end